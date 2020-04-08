function [B0Data] = op_ReadRecoGREAndCalcB0Maps(GREData_file,TargetPar,Settings)
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
if(~isfield(Settings,'Volunteer_flag'))
    Settings.Volunteer_flag = true;
end




%%   

 
    test2 = mapVBVD(GREData_file); test2 = test2.image();
    test2 = permute(test2,[2 1 3 5 8 4 6 7]);
    test{1} = double(test2(:,:,:,:,1)); test{2} = double(test2(:,:,:,:,2)); clear test2;
    B0Data.Par = read_ascconv(GREData_file);

    % FFT B0Data-Data
    test{1} = conj(test{1}); test{2} = conj(test{2});
    test{1} = FFTOfMRIData(test{1},0,[2 3],1,1,0)*10^6;
    test{2} = FFTOfMRIData(test{2},0,[2 3],1,1,0)*10^6;

    % Cut FoV if D2GoldStandard and B0Data FoVs are not the same. Generalize that later...
    if(abs(TargetPar.FoV_Read(1) - B0Data.Par.FoV_Read(1)) > 10)
        Sizee = size(test{1});
        test{1} = ZerofillOrCutkSpace(test{1},[Sizee(1) 126 126 Sizee(4)],0);   % The FoV was by mistake 240x240, while the CSI is 220x220
        test{2} = ZerofillOrCutkSpace(test{2},[Sizee(1) 126 126 Sizee(4)],0);   % Therefore need to cut down to 138*220/240 = 126.5 ~ 126

%         % Cut FoV - Method2: Interpolate to a 132x132 matrix
%         Sizee = size(test{1});
%         test2{1} = zeros([Sizee(1) 132*2 132 Sizee(4)]);
%         test2{2} = zeros([Sizee(1) 132*2 132 Sizee(4)]);
%         for ii=1:Sizee(1)
%             for jj=1:Sizee(4)
%                 test2{1}(ii,:,:,jj) = imresize(squeeze(test{1}(ii,:,:,jj)),[132*2 132]);
%                 test2{2}(ii,:,:,jj) = imresize(squeeze(test{2}(ii,:,:,jj)),[132*2 132]);            
%             end
%         end
%         test{1} = ZerofillOrCutkSpace(test2{1},[Sizee(1) 121 121 Sizee(4)],0);   % The FoV was by mistake 240x240, while the CSI is 220x220
%         test{2} = ZerofillOrCutkSpace(test2{2},[Sizee(1) 121 121 Sizee(4)],0);   % Therefore need to cut down to 132*220/240 = 121 ~ 121
%         clear test2

    else
        Sizee = size(test{1});
        test{1} = ZerofillOrCutkSpace(test{1},[Sizee(1) Sizee(2)/2 Sizee(3:end)],0);   % The FoV was by mistake 240x240, while the CSI is 220x220
        test{2} = ZerofillOrCutkSpace(test{2},[Sizee(1) Sizee(2)/2 Sizee(3:end)],0);   % Therefore need to cut down to 138*220/240 = 126.5 ~ 126
    end
    
    
%     % Flip LR
%     test{1} = flip(test{1},3);
%     test{2} = flip(test{2},3);

    % Save data
    B0Data.Data = test;
    clear test

    B0Data = op_ReorderSlices(B0Data);
    
    
    % Average slices
    SlcLimits = [TargetPar.Pos_Tra-TargetPar.VoI_Partition/2, TargetPar.Pos_Tra+TargetPar.VoI_Partition/2]; 
    SlcNo = cell2mat(FindClosestIndex(B0Data.Par.Pos_Tra,SlcLimits));
    SlcNo = [max(SlcNo(:,1)) min(SlcNo(:,2))];
    
    B0Bak = B0Data.Data{1};
    B0Data.Data{1} = mean(B0Data.Data{1}(:,:,:,SlcNo(1):SlcNo(2)),4);
    B0Data.Data{2} = mean(B0Data.Data{2}(:,:,:,SlcNo(1):SlcNo(2)),4);

    % Coil combination of B0map
    B0Data.B0Map = circshift(squeeze(sum(B0Data.Data{1} .* conj(B0Data.Data{2}))),[-1 0]);

    % Calculate Mask
    B0Data.Mag = abs(B0Data.B0Map); 

    % Calculate Mask 
    if(Settings.Volunteer_flag)
        % Using bet
        B0Mag = squeeze(sqrt(sum(abs(B0Bak).^2)));
        PixSiz = [B0Data.Par.FoV_Phase(1)/B0Data.Par.nPhasEnc ...
            B0Data.Par.FoV_Read(1)/B0Data.Par.nFreqEnc B0Data.Par.FoV_Partition(1)/(B0Data.Par.nPartEnc*B0Data.Par.nSLC)];
        BetMask = mrir_brain_mask__BET(B0Mag,PixSiz,Settings.BetPath,Settings.BetOptions);
        LipMask = logical(BetMask.outskin_mask - BetMask.mask);
        Mask = logical(BetMask.mask); clear BetMask
        Mask = Mask(:,:,SlcNo(1):SlcNo(2));
        Mask = imresize3(double(Mask),[size(Mask,1) size(Mask,2) 1],'nearest');
        Mask = MaskShrinkOrGrow(Mask,2,0,1);
        B0Data.BrainMask = logical(Mask);
        LipMask = LipMask(:,:,SlcNo(1):SlcNo(2));
        LipMask = imresize3(double(LipMask),[size(LipMask,1) size(LipMask,2) 1],'nearest');
        LipMask = MaskShrinkOrGrow(LipMask,2,1,1);
        B0Data.LipMask = logical(LipMask);     
        B0Data.Mask = B0Data.LipMask + B0Data.BrainMask;
    else
        % Using threshold
        madd = mad(B0Data.Mag(:));
        medi = median(B0Data.Mag(:));
        B0Data.Mask = B0Data.Mag > medi + 0.1*madd;
        B0Data.Mask = MaskShrinkOrGrow(MaskShrinkOrGrow(B0Data.Mask,3,1,1),2,0,1);
        B0Data.LipMask = zeros(size(B0Data.Mask));
        
    end
    clear B0Bak Mask B0Mag PixSiz madd medi;

    

    % Calculate Frequency Offset
    B0Data.B0Map = angle(B0Data.B0Map)/((B0Data.Par.TEs(2)-B0Data.Par.TEs(1))/10^6*2*pi);

    % Mask B0Map
    B0Data.B0Map = B0Data.B0Map .* B0Data.Mask;

%% Postparations


B0Data = supp_UpdateRecoSteps(B0Data,Settings);



