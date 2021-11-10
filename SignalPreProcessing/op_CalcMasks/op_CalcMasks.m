function [MRStruct] = op_CalcMasks(MRStruct,Settings)
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
if(~isfield(Settings,'BET_flag'))
    Settings.BET_flag = true;
end


%% Calculate Mask 

if(Settings.BET_flag)
    % Using bet
    B0Mag = sqrt(sum(abs(MRStruct.Data).^2,MRStruct.RecoPar.dims.cha)) .* MRStruct.Masks.Mask;
    
    if(isfield(MRStruct.RecoPar,'ContrastNames') && any(strcmpi(MRStruct.RecoPar.ContrastNames,'UNI')))
        TakeContrast = find(strcmpi(MRStruct.RecoPar.ContrastNames,'UNI'));
        S.type = '()'; S.subs = repmat({':'},[1 numel(size(B0Mag))]); S.subs{MRStruct.RecoPar.dims.contrast} = TakeContrast;
        B0Mag = subsref(B0Mag,S); clear S TakeContrast
    end
    
    PixSiz = [MRStruct.Par.FoV_Phase(1)/MRStruct.Par.nPhasEnc ...
        MRStruct.Par.FoV_Read(1)/MRStruct.Par.nFreqEnc MRStruct.Par.FoV_Partition(1)/(MRStruct.Par.nPartEnc*MRStruct.Par.nSLC)];
    BetMask = mrir_brain_mask__BET(B0Mag,PixSiz,Settings.BetPath,Settings.BetOptions);
    LipMask = logical(BetMask.outskin_mask - BetMask.mask);
    Mask = logical(BetMask.mask); clear BetMask

    MRStruct.Masks.BrainMask = logical(Mask);

    MRStruct.Masks.LipMask = logical(LipMask);     
    MRStruct.Masks.Mask = MRStruct.Masks.LipMask + MRStruct.Masks.BrainMask;
    
else
    % Using threshold
    madd = mad(MRStruct.Mag(:));
    medi = median(MRStruct.Mag(:));
    MRStruct.Masks.Mask = MRStruct.Mag > medi + 0.1*madd;
    MRStruct.Masks.Mask = MaskShrinkOrGrow(MaskShrinkOrGrow(MRStruct.Masks.Mask,3,1,1),2,0,1);
    MRStruct.Masks.LipMask = zeros(size(MRStruct.Masks.Mask));

end
clear B0Bak Mask B0Mag PixSiz madd medi;



%% Postparations


MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);



