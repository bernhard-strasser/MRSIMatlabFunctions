function [SensMap,GREData] = op_ReadRecoGREAndCalcSensMaps(GREData_file,Mask,TargetPar,Settings)
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
if(~isfield(Settings,'Imresize_flag'))
    Settings.Imresize_flag = false;
end



%%    
SensMap.Par = read_ascconv(GREData_file);


Temp = mapVBVD(GREData_file);
Temp2 = Temp.image{''};

% Best would be to determine if acquisition was done transversal, coronal, or sagittal. But for some reason that is not given directly in the 
% header of raw data?!?!?!?! Therefore just handle the cases of a B0-map or a CoilAdjScan (first is multi-slice, other is 3D, that's how I
% distinguish here)
if(SensMap.Par.ThreeD_flag)
    permuteVec = [2 4 3 1 5];
    ReadDir = find(permuteVec == 1);
    OverSampVec = ones([1 3]); OverSampVec(ReadDir-1) = 2;
else
    permuteVec = [2 1 3 4 5];
    ReadDir = find(permuteVec == 1);
    OverSampVec = ones([1 3]); OverSampVec(ReadDir-1) = 2;
end

Temp2 = permute(Temp2,permuteVec);
if(~SensMap.Par.ThreeD_flag && SensMap.Par.InterleavedSliceAcquisition)
    SliceOrder = cat(2,1:2:size(Temp2,4),2:2:size(Temp2,4));
    [dum,SliceOrderedInd] = sort(SliceOrder);
    Temp2 = Temp2(:,:,:,SliceOrderedInd,:,:,:,:);
end
SensMap.Data.AC = Temp2(:,:,:,:,1);
SensMap.Data.BC = Temp2(1,:,:,:,2);
    
SensMap.RecoPar = SensMap.Par;
SensMap.InputPar = SensMap.Par;

if(SensMap.Par.ThreeD_flag)
    FFTVec = [2 3 4];
else
    FFTVec = [2 3];
end

dummy1 = FFTOfMRIData(SensMap.Data.AC,0,FFTVec,1,1,0);
dummy2 = FFTOfMRIData(SensMap.Data.BC,0,FFTVec,1,1,0);
if(~SensMap.Par.ThreeD_flag)         % WHY IS THIS NECESSARY?!
    dummy1 = flip(flip(dummy1,2),3);
    dummy2 = flip(flip(dummy2,2),3);    
end

GREData = SensMap;
GREData.Data = dummy1;
TargetPar.DataSize = [TargetPar.nFreqEnc TargetPar.nPhasEnc TargetPar.nPartEnc*TargetPar.nSLC TargetPar.vecSize];


%% Calc Sensmap
SensMap = op_CalcSensMaps(GREData,[],Mask,TargetPar,struct('Extrapolate_flag',true));


%% Imresize SensMaps if Wanted

if(Settings.Imresize_flag)
    DummyStruc.Par = TargetPar;
    DummyStruc.Data = ones(TargetPar.DataSize);
    CurSetts.CoilCombine_D2.ResizeMethod = 'Imresize';
    CurSetts.CoilCombine_D2.ScalingMethod = 'UniformSensitivity';  
    [DummyStruc AdditionalOut] = op_CoilCombineData(DummyStruc,SensMap,CurSetts.CoilCombine_D2);
    
end

%% Postparations


SensMap = supp_UpdateRecoSteps(SensMap,Settings);

