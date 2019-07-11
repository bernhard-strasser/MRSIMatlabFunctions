function [Output, AdditionalOut] = op_IterReconstructNonCartMRData(Output,AdditionalIn,ModelFunction,Settings)
%
% read_csi_dat Read in raw data from Siemens
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         file                    ...     Path of MRS(I) file.
% -         DesiredSize             ...     If you want to perform zerofilling or use only a part of the kspace, set
%                                           DesiredSize to your wanted [kx,ky,kz], e.g. [64 64 1].
% -         ReadInDataSets          ...     
%
% Output:
% -         kSpace                      ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x Samples x Averages
% -         Info                        ...     
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


%% 0. Preparations



if(~exist('Settings','var'))
   Settings.Phaseroll_flag = true;
   Settings.DensComp_flag = true;   
    
end
if(~isfield(Settings,'Phaseroll_flag'))
   Settings.Phaseroll_flag = true;    
end
if(~isfield(Settings,'DensComp_flag'))
   Settings.DensComp_flag = true;    
end
if(~isfield(Settings,'DensComp'))
   Settings.DensComp = struct();    
end
if(~isfield(Settings,'ConjInkSpace_flag'))
   Settings.ConjInkSpace_flag = false;    
end
if(~isfield(Settings,'ConjIniSpace_flag'))
   Settings.ConjIniSpace_flag = true;    
end
if(~isfield(Settings,'Correct4SpatialB0_flag'))
   Settings.Correct4SpatialB0_flag = false;    
end
if(~isfield(Settings,'Correct4SpectralB0_flag'))
   Settings.Correct4SpectralB0_flag = false;    
end

if(~isfield(Settings,'CircularSFTFoV_flag'))
   Settings.CircularSFTFoV_flag = false;    
end
if(~isfield(Settings,'DensCompAutoScale_flag'))
   Settings.DensCompAutoScale_flag = false;    
end

% if(~isfield(InData,'RecoPar'))
%     if(~isfield(InData,'Par'))
%         error('InData must have field Par or RecoPar.')
%     end
%     InData.RecoPar = InData.Par;
% end


%% FOV SHIFTs

% NOT YET IMPLEMENTED

% THATS HOW LUKAS HINGERL DOES IT (CODE FROM LUKAS HINGERL):
% %Inplane FOV shift
% kspace_max=sqrt(kx(1,end).^2+ky(1,end).^2);
% yIPFOVshift=-ReadInInfo.General.Ascconv.PosVOI_Sag; %this should be right!
% xIPFOVshift=-ReadInInfo.General.Ascconv.PosVOI_Cor;
% zIPFOVshift=ReadInInfo.General.Ascconv.PosVOI_Tra;
% FOV=ReadInInfo.General.Ascconv.FoV_Phase;
% FOVz=ReadInInfo.General.Ascconv.FoV_Partition;
% 
% 
% for circles=1:nc
%     for samplepoint=1:ns
%         csi_k(:,circles,samplepoint,1,:)=csi_k(:,circles,samplepoint,1,:)*exp(1i*kx(samplepoint,circles)/kspace_max*((REShalbe-0.5)/FOV)*2*pi*xIPFOVshift)*exp(1i*ky(samplepoint,circles)/kspace_max*((REShalbe-0.5)/FOV)*2*pi*yIPFOVshift);
%     end
% end



%% Conj in Beginning

if(Settings.ConjInkSpace_flag)
    Output.Data = conj(Output.Data);
    if(isfield(Output,'NoiseData'))
        Output.NoiseData = conj(Output.NoiseData);
    end
end

%% Try to correct for "tilted trajectory" by doing phaseroll

% Save the data for reconstructing the pseudo-pcg case
if(Settings.Phaseroll_flag)

    nTI = Output.Par.nTempInt;
    vs = Output.Par.vecSize;
    ns = Output.Par.TrajPts;
    nc = Output.Par.nAngInts;
    nrew = Output.Par.RewPts;
    ncha = size(Output.Data,6);

    timeoffset = 0:(ns-1);
    timeoffset = repmat(transpose(timeoffset),[1 vs]);
    Freq = ((0:vs-1)/vs-0.5)/(nrew + ns)*nTI;
    Freq = repmat(Freq,[ns 1]);    

    % This comes from:
    % timeoffset = (0:(ns-1))*Output.Par.ADC_Dt/10^6;
    % sBW = nTI/((nrew + ns)*Output.Par.ADC_Dt/10^6);
    % Freq = -sBW/2 : sBW/vs : (sBW/2 - sBW/vs);
    % Output.Par.ADC_Dt/10^6 cancels out when calculating timeoffset * Freq and so can be omitted
    % the rest is basically the same (-sBW/2:sBW/vs:(sBW/2-sBW/vs) is equivalent to ((0:vs-1)/vs-0.5), and the other constants are
    % the same anyway
    
    phasecorr = exp(-2*1i*pi*timeoffset .* Freq);    
    phasecorr = myrepmat(phasecorr,size(Output.Data));
    

    TiltTrajMat = (phasecorr);   
   
    
else
    
	TiltTrajMat = ones([Output.Par.TrajPts*Output.Par.nAngInts Output.Par.vecSize Output.Par.DataSize(end)]);

end

AdditionalOut.TiltTrajMat = TiltTrajMat;



%% Calculate sft2-Operator

% sft Operator
% Reshape trajectories to expected shape

% Collapse data to a matrix (from [nAngInt x nTrajPoints x nTempInt*vecSize x nCha x nPart*nSlc] to [nAngInt*nTrajPoints x Rest])
SizeData = size(Output.Data); SizeData = cat(2,SizeData,ones([1 5-numel(SizeData)]));



sft2_Oper = sft2_Operator(transpose(squeeze(Output.OutTraj.GM(:,:))*Output.RecoPar.DataSize(1)),transpose(Output.InTraj.GM(:,:)),1);

% Restrict to circular FoV
if(Settings.CircularSFTFoV_flag)
    FoVMask = EllipticalFilter(ones(Output.RecoPar.DataSize(1:2)),[1 2],[1 1 1 Output.RecoPar.DataSize(1)/2-1],1); 
    FoVMask = FoVMask(:);
    sft2_Oper(:,~logical(FoVMask(:))) = 0;
    clear FoVMask;
end


%% Calculate B0-Correction of Spiral Data in Spatial Domain
if(Settings.Correct4SpatialB0_flag)
    t   = (0:Output.RecoPar.TrajPts-1)*Output.RecoPar.ADC_dt/10^9;
    t = repmat(t,[1 1 Output.RecoPar.nAngInts]); t = t(:);
    CurB0 = imresize(AdditionalIn.B0.B0Map,Output.RecoPar.DataSize(1:2));    
    if(isfield(AdditionalIn.B0,'Mask'))
        Mask = imresize(MaskShrinkOrGrow(AdditionalIn.B0.Mask,2,0,1),Output.RecoPar.DataSize(1:2),'nearest');
        CurB0 = CurB0 .* Mask;
    end
    
    AdditionalOut.B0CorrMat_Spatial = exp(transpose(-2*pi*1i*CurB0(:)) .* t);

    sft2_Oper = sft2_Oper .* AdditionalOut.B0CorrMat_Spatial;
    
end

%% Calculate B0CorrMat_Spec

Mask = imresize(AdditionalIn.B0.Mask,Output.RecoPar.DataSize(1:2),'nearest');
AdditionalOut.Mask = Mask;
if(Settings.Correct4SpectralB0_flag)
    % Go to resolution of D1
    CurB0 = imresize(AdditionalIn.B0.B0Map,Output.RecoPar.DataSize(1:2)) ;
    if(exist('Mask','var'))
        CurB0 = CurB0 .* Mask;
    end

    % Round to shift only integer number of points
    HzPerPt = 10^9/Output.Par.Dwelltimes(1) / Output.Par.vecSize;
    % if(Settings.RoundB0ToIntVecPts)
    %     CurB0 = round(CurB0/HzPerPt)*HzPerPt;
    % end

    time   = (0:Output.Par.DataSize(5)-1)*Output.Par.Dwelltimes(1)/10^9;
    time = repmat(time,[Output.Par.DataSize(end) 1]);
    B0CorrMat_Spec = exp(myrepmat(2*pi*1i*CurB0(:),[prod(Output.RecoPar.DataSize(1:2)),Output.RecoPar.DataSize(4:end)]) .* myrepmat(time,[prod(Output.RecoPar.DataSize(1:2)),Output.RecoPar.DataSize(4:end)]));


    AdditionalOut.B0CorrMat_Spec = B0CorrMat_Spec;
else
    AdditionalOut.B0CorrMat_Spec = ones([prod(Output.RecoPar.DataSize(1:2)),Output.RecoPar.DataSize(4:end)]);
end


%% Calculate Density Compensation According to Hoge1997 - Abrupt Changes


if(Settings.DensComp_flag)
    Dummy = Output; Dummy = rmfield(Dummy,'Data'); Dummy = rmfield(Dummy,'NoiseData');
    [Dummy,Dummy2] = op_CalcAndApplyDensComp(Dummy,sft2_Oper,Settings.DensComp);
    AdditionalOut.DCFPreG = Dummy2.DCFPreG;
    clear Dummy;
end
% Output.Data = reshape(Output.Data,[prod(SizeData_k(1:2)) prod(SizeData_k(3:end))]);  

%% Calculate SensitivityMap
% Fake it
Dummy.ResizeMethod = 'Imresize';
Dummy.ScalingMethod = 'UniformSensitivity';   Dummy2.Data = ones(Output.RecoPar.DataSize); Dummy2.RecoPar = Output.RecoPar; %Dummy2.Par = Output.Par; 
[DeleteMe, AddOut] = op_CoilCombineData(Dummy2,AdditionalIn.SensMap,Dummy);
AdditionalOut.SENSE = squeeze_single_dim(AddOut.CoilWeightMap(:,:,:,1,:),3);
clear DeleteMe AddOut Dummy;


Scale = 1 ./ sqrt(sum(abs(AdditionalOut.SENSE(:,:,1,:)).^2,4));
Scale(isinf(Scale) | isnan(Scale) | Scale == 0 ) = 1;
AdditionalOut.SENSE = AdditionalOut.SENSE(:,:,1,:).*Scale.*AdditionalOut.Mask;
AdditionalOut.SENSE(isinf(AdditionalOut.SENSE) | isnan(AdditionalOut.SENSE)) = 0;



AdditionalOut.sft2_Oper = sft2_Oper;
AdditionalOut.TiltTrajMat = reshape(AdditionalOut.TiltTrajMat,[prod(Output.Par.DataSize(1:4)) Output.Par.DataSize(5:end)]);
AdditionalOut.DCFPreG = AdditionalOut.DCFPreG(:);


%% Iterative Reconstruction of CSI Data


Opers = AdditionalOut;
Opers.SamplingOperator = ones(size(Output.Data(:,:,1)));
Opers.B0CorrMat_Spec = reshape(Opers.B0CorrMat_Spec(:,:,1),Output.RecoPar.DataSize([1 2 4]));
Opers.TiltTrajMat = reshape(Opers.TiltTrajMat(:,:,1),[size(Opers.SamplingOperator) size(Opers.TiltTrajMat,2)]);
Opers.DCFPreG = reshape(Opers.DCFPreG,size(Opers.SamplingOperator));


AO = @(x) ModelFunction('NoTransj',x,Opers);
AOT = @(x) ModelFunction('Transj',x,Opers);

[m,n,p] = size(AdditionalOut.B0CorrMat_Spec);
maxiter = 20;
lambda = 0.001;%1e-16;%0.001;
eta =20e-2;

% Output.Data = lr_method(Output.Data,AO,AOT,maxiter,lambda,eta,m,n,Settings, Output.RecoPar);
Output.Data = LS_method(Output.Data,AO,AOT,maxiter,eta);

Output.Data = reshape(Output.Data,[Settings.fov_overgrid*Output.RecoPar.DataSize(1:2) SizeData(3:end-1)]);

bla = ([size(Output.Data,1) size(Output.Data,2)] - Output.RecoPar.DataSize(1:2))/2+1;
Output.Data = Output.Data(bla(1):bla(1)+Output.RecoPar.DataSize(1)-1,bla(1):bla(2)+Output.RecoPar.DataSize(2)-1,:,:,:,:);

% Reconstruct noise-data non-iteratively. Cannot reconstruct noise iteratively, or could we? We actually would want to perform exactly the same reconstruction
% as for the normal data...
if(isfield(Output,'NoiseData'))
    Output.NoiseData = AOT(Output.NoiseData);
end


%% Perform Reconstruction in Slice and z-dimension

% For now just reshape them. We dont have slices or 3D-measurements for now...
Size = size(Output.Data);
Output.Data = reshape(Output.Data, [Size(1:2) prod(Size(3:4)) Size(5:end)]);
if(isfield(Output,'NoiseData'))
    Output.NoiseData = reshape(Output.NoiseData, [Size(1:2) prod(Size(3:4)) Size(5:end)]);
end


%% Conj at End

if(Settings.ConjIniSpace_flag)
    Output.Data = conj(Output.Data);
    if(isfield(Output,'NoiseData'))
        Output.NoiseData = conj(Output.NoiseData);
    end
end

%% Flip left right

% Output.Data = flip(Output.Data,2);


%% Create & Adapt Parameters


Output.RecoPar.DataSize = size(Output.Data);
Output.RecoPar.total_channel_no_reco = 1;      % I actually should make an if-condition to determine if coil combination was really done...



%% Postparations

Output = supp_UpdateRecoSteps(Output,Settings);



