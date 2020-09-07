function [SensMap,AdditionalOut] = op_CalcSensMaps(SensMapData,AdditionalMaskData,TargetPars,Settings)
%
% op_ReadAndRecoBorjanSpiralData Read and reconstruct data from Borjan Gagoski's Spiral MRSI Sequence
%
% This function was written by Bernhard Strasser, June 2019.
%
%
% The function can read in Spiral MRSI data in the Siemens raw file format ".DAT" and performs
% the reconstruction of the data (Non-Uniform Slow FourierTransform etc.)
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         SpiralDataFile          ...     
% -         SpiralTrajectoryFile    ...     
% -         Settings                ...     
%
% Output:
% -         ?                      ...     
% -         ?                        ...     
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.

% This function expects the input to be of form
% [nCha, nAngInt 

% Input data must be structure with at least fields 'Data' and 'Par'.
% 'Data' must be of size 



%% 0. Preparations

if(~exist('Settings','var'))
    Settings = struct;
end
if(~isfield(Settings,'Debug_ShowSensMaps_flag'))
    Settings.Debug_ShowSensMaps_flag = false;
end
if(~isfield(Settings,'Extrapolate_flag'))
    Settings.Extrapolate_flag = false;
end
if(~isfield(Settings,'RatioToMax4MaskThreshold'))
    Settings.RatioToMax4MaskThreshold = 0.035;
end

if(exist('AdditionalMaskData','var') && isempty(AdditionalMaskData))
    clear AdditionalMaskData;
end
% if(~isfield_recursive(Settings,'Debug.ShowTrajs'))
%     Settings.Debug.ShowTrajs = true;
% end



%% Resize & Cut Image to Fit Target Data
% The image used for the sensitivity map calculation might have different FoV and matrix size than the data which should be coil combined.
% Therefore need to resample the SensMapImage to the target image.
% CAUTION: THE ANGULATION OF THE SensMapImage (AC) AND THE TARGET IMAGE NEED TO BE THE SAME!
% This is done here, although I think a bit complicated. It would be probably easier to interpolate using the original grid and the new grid, and then
% use some interpolation method. This way, even angulations could be handled.

if(TargetPars.FoV_Read ~= SensMapData.Par.FoV_Read(1) || TargetPars.FoV_Phase ~= SensMapData.Par.FoV_Phase(1) || TargetPars.FoV_Partition ~= SensMapData.Par.FoV_Partition(1) || ...
   TargetPars.Pos_Sag(1) ~= SensMapData.Par.Pos_Sag(1) || TargetPars.Pos_Cor(1) ~= SensMapData.Par.Pos_Cor(1) || TargetPars.Pos_Tra(1) ~= SensMapData.Par.Pos_Tra(1) || ...
   TargetPars.nFreqEnc ~= SensMapData.Par.nFreqEnc || TargetPars.nPhasEnc ~= SensMapData.Par.nPhasEnc || TargetPars.nPartEnc ~= SensMapData.Par.nPartEnc || ...
   TargetPars.nSLC ~= SensMapData.Par.nSLC || ...
   TargetPars.SliceNormalVector_x ~= SensMapData.Par.SliceNormalVector_x || TargetPars.SliceNormalVector_y ~= SensMapData.Par.SliceNormalVector_y || ...
   TargetPars.SliceNormalVector_z ~= SensMapData.Par.SliceNormalVector_z || TargetPars.InPlaneRotation ~= SensMapData.Par.InPlaneRotation)

    RoundingFactor = 1;         % Maximum error for calculated FoV: 0.5/RoundingFactor mm
    MatSize_SensMap = size_MultiDims(SensMapData.Data,2:4);
    FoV_D2 = [TargetPars.FoV_Read TargetPars.FoV_Phase TargetPars.FoV_Partition];
    FoV_SensMap = [SensMapData.Par.FoV_Read(1) SensMapData.Par.FoV_Phase(1) SensMapData.Par.FoV_Partition(1)];

    % % Take care of oversampling. THIS IS A HACK! HOW CAN I BETTER FIND OUT WHICH DIRECTION WAS ZEROFILLED?!?!?
    % if(SensMapData.Par.ThreeD_flag && SensMapData.Par.nPartEnc > 1)  
    %     FoV_SensMap(3) = 2*FoV_SensMap(3);
    % else
    %     FoV_SensMap(1) = 2*FoV_SensMap(1);    
    % end

    LeastCommMult = lcm(round(FoV_D2*RoundingFactor),round(FoV_SensMap*RoundingFactor));
    ImresizeSize = LeastCommMult ./ round(FoV_D2*RoundingFactor);
    ImresizeSize = ceil(MatSize_SensMap./ImresizeSize) .* ImresizeSize;        % Dont allow TargetSize < InputSize
    TargetSize = round(FoV_D2 ./ (FoV_SensMap./ImresizeSize));

    dummy3 = zeros([size(SensMapData.Data,1) ImresizeSize]);
    for cha = 1:size(SensMapData.Data,1)
        dummy3(cha,:,:,:) = imresize3(squeeze(SensMapData.Data(cha,:,:,:)),ImresizeSize);        
    end
    SensMapData.Data = dummy3;

    % In case the FoV was too small, zerofill in image-domain
    SensMapData.Data = ZerofillOrCutkSpace(SensMapData.Data,[size(SensMapData.Data,1) TargetSize],0);


    % Imresize in z-dimension
    TargetSize_z = [TargetSize(1:2) TargetPars.DataSize(3)];
    Temp = zeros([size(SensMapData.Data,1) TargetSize_z]);
    for cha = 1:size(SensMapData.Data,1)
        Temp(cha,:,:,:) = imresize3(squeeze(SensMapData.Data(cha,:,:,:)),TargetSize_z);
    end
    SensMapData.Data = Temp;



    % % Resample Data using griddedInterpolant
    % % Assume input points at locations [1,2,...,N_i] for dimension i
    % % Points to which we want to interpolate:
    % 
    % MatSzDiff = - MatSize_SensMap;
    % 
    % xq = (0:5/6:sx)';
    % yq = (0:5/6:sy)';
    % zq = (1:sz)';
    % for cha = 1:size(AC.Data,1)
    %     F = griddedInterpolant(AC.Data(cha,:,:,:));
    %     
    % end

    % Dont know if that is necessary
    SensMapData.Data = circshift(SensMapData.Data,[0 -1 -1 0]);

end
    

%% Calculate Mask

if(exist('AdditionalMaskData','var') && sum(AdditionalMaskData(:) == 0) + sum(AdditionalMaskData(:) == 1) == numel(AdditionalMaskData))% If AdditionalMaskData is mask
    Maskk = imresize(AdditionalMaskData,[size(SensMapData.Data,2) size(SensMapData.Data,3)],'nearest');
else
    ACSoS = squeeze_single_dim(sqrt(sum(abs(SensMapData.Data).^2)),1);
    Maskk = ACSoS > Settings.RatioToMax4MaskThreshold*max(max(max(ACSoS)));     % Kind of arbitrary for now... Better mask creation would be useful!
    for slc = 1:size(SensMapData.Data,4)
        Maskk(:,:,slc) = imfill(Maskk(:,:,slc),'holes');  % Fill holes
    end
end


SensMap.Mask = reshape(Maskk,[1 size(Maskk)]); 


%% Prepare and Perform ESpirit


% ESpirit seems to only handle 2D-data. Thus do the reco slice-by-slice or partition-by-partition

SensMap.Data = zeros(size(SensMapData.Data));
for slc = 1:size(SensMapData.Data,4)

    CurData = SensMapData.Data(:,:,:,slc);
    Data4Espirit.Water_ctkk = reshape(FFTOfMRIData(CurData,0,[2 3],0,1,0),[size(CurData,1) 1 size_MultiDims(CurData,[2 3])]);
    Data4Espirit.NbPtForWaterPhAmp = 1;
    Data4Espirit.Water_kmask = ones(size_MultiDims(CurData,[2 3]));
    Data4Espirit.ImMask2D = squeeze(Maskk(:,:,slc));

    SensMap.Data(:,:,:,slc) = ComputeSENSEProfilesWithESPIRiT_bstr(Data4Espirit);

end
SensMap.Par = SensMapData.Par;
SensMap.Par.DataSize = size(SensMap.Data);


%% Extrapolate a bit

if(Settings.Extrapolate_flag)
    Maskk_NaN = double(SensMap.Mask); Maskk_NaN(Maskk_NaN == 0) = NaN;
    SensMap.Data = SensMap.Data .* Maskk_NaN;
    for cha = 1:size(SensMap.Data,1)
        SensMap.Data(cha,:,:) = single(inpaint_nans(double(squeeze(SensMap.Data(cha,:,:))))); 
    end
    Maskk_Big = reshape(MaskShrinkOrGrow(squeeze(SensMap.Mask),2,1,1),size(SensMap.Mask));
    SensMap.Data = Maskk_Big .* SensMap.Data;
    SensMap.Mask = Maskk_Big;
    
    SensMap.Data(isnan(SensMap.Data)) = 0;
    SensMap.Mask(isnan(SensMap.Mask)) = 0;
   
end


%% Normalize ESpirit Maps

%     for ii=1:32; Normy(ii) = norm(squeeze(SensMap.Data(ii,:,:))); end
%     SensMap.Data = SensMap.Data ./ myrepmat(Normy,size(SensMap.Data)); clear Normy
% blubb = abs(SensMap.Data(end,:)); blubb(blubb == 0) = NaN;
% SensMap.Data = SensMap.Data(1:end-1,:,:,:) / nanmean(blubb)*1;
% SensMap.Data = SensMap.Data(1:end-1,:,:,:);


%% Debug: Show SensMaps

if(Settings.Debug_ShowSensMaps_flag)
    TotChannels = size(SensMap.Data,1);
    Slc = ceil(size(SensMap.Data,4)/2);
    Rows = floor(sqrt(TotChannels));
    Cols = ceil(sqrt(TotChannels));
    Rows(Rows*Cols < TotChannels) = Rows+1;
    figure; 
    maxi = nanmax(abs(SensMap.Data(:)));
    for CurCha = 1:TotChannels
        subplot(Rows,Cols,CurCha)
        imagesc(squeeze(abs(SensMap.Data(CurCha,:,:,Slc))),[0 maxi])
        colorbar
        title(['SensMap Slc ' num2str(Slc) ', Cha ' num2str(CurCha)])
    end

end
    

%% Permute

SensMap = op_PermuteMRData(SensMap,[2 3 4 5 1]);


%% Postparations


SensMap = supp_UpdateRecoSteps(SensMap,Settings);

