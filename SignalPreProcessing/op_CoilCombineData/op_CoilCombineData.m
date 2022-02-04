function [MRStruct, AdditionalOut] = op_CoilCombineData(MRStruct,CoilWeightMap,Settings)
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
% -         MRStruct                    ...     The Input. The structure containing the MR data. Can have different fields. For this function, it
%                                               must have fields 'Data', and 'Par'.
% -         CoilWeightMap               ...     The Coil-weighting-map. Must have field 'Data'
% -         Settings                    ...     Structure to specify how the coil combination should be performed.
%
% Output:
% -         MRStruct                    ...     The Output (Input is overwritten). Will have field 'RecoPar' additionally.
% -         AdditionalOut               ...     Variable for storing additional output (if asked for by user).   
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: For now, the coil index in CoilWeightMap is the first, but in MRStruct it's the last. Change that later!
% Also, the function cannot handle 3D / multislice data for now. And the channel index is fixed for now.


% This function expects the input to be of form
% ??? 


%% 0. Preparations

if(~exist('Settings','var'))
    Settings = struct;
end
if(isfield(Settings,'ResizeMethod') && strdist(Settings.ResizeMethod,'Imresize') < strdist(Settings.ResizeMethod,'Zerofill'))
    ImResize_flag = true;
else
    ImResize_flag = false;
end
if(isfield(Settings,'ScalingMethod') && strdist(Settings.ScalingMethod,'UniformNoise') < strdist(Settings.ScalingMethod,'UniformSensitivity'))
    UniformSignal_flag = false;
else
    UniformSignal_flag = true;
end
if(~isfield(MRStruct,'Par'))
    MRStruct.Par = struct;
    if(isfield(MRStruct,'Data_file'))
        MRStruct.Par = read_ascconv(MRStruct.Data_file); 
    end
end
if(~isfield(MRStruct.Par,'DataSize'))
    MRStruct.Par.DataSize = size(MRStruct.Data);
end
if(~isfield(MRStruct,'RecoPar'))
    MRStruct.RecoPar = MRStruct.Par;
end

% SoS
if(~exist('CoilWeightMap','var'))
    CoilWeightMap = MRStruct;
    CoilWeightMap.Data = conj(CoilWeightMap.Data(:,:,:,1,:));
    CoilWeightMap.Mask = ones(size_MultiDims(CoilWeightMap.Data,[1 2 3])); 
    UniformSignal_flag = false;
    Settings.ScalingMethod = 'UniformNoise';
end

if(~isfield(CoilWeightMap,'Mask'))
    CoilWeightMap.Mask = ones(size_MultiDims(CoilWeightMap.Data,[1 2 3])); 
end


if(numel(CoilWeightMap.Data) == 1)
    AdditionalOut.CoilWeightMap = CoilWeightMap.Data;
    AdditionalOut.Mask = 1;
    AdditionalOut.Scaling = 1;
    return;
end


%% Resize SensMap

if(ImResize_flag)
    AdditionalOut.CoilWeightMap = zeros([MRStruct.RecoPar.DataSize(1:2) size_MultiDims(CoilWeightMap.Data,3:5)]);
    for cha = 1:prod(size_MultiDims(CoilWeightMap.Data,3:5))
        AdditionalOut.CoilWeightMap(:,:,cha) = imresize(CoilWeightMap.Data(:,:,cha),MRStruct.RecoPar.DataSize(1:2));
    end
else
    AdditionalOut.CoilWeightMap = ZerofillOrCutkSpace(CoilWeightMap.Data,[MRStruct.RecoPar.DataSize(1:3) size_MultiDims(CoilWeightMap.Data,4:5)],1);
end
if(MRStruct.RecoPar.DataSize(3) == 1)      % imresize3 for some reason doesnt work for 2D-input...
    AdditionalOut.Mask = imresize(squeeze(CoilWeightMap.Mask(:,:)),MRStruct.RecoPar.DataSize(1:2),'nearest');       
else 
    AdditionalOut.Mask = logical(imresize3(squeeze(uint8(CoilWeightMap.Mask(:,:,:))),MRStruct.RecoPar.DataSize(1:3),'nearest'));       
end

%% Weight Data

MRStruct.Data = MRStruct.Data .* AdditionalOut.CoilWeightMap;
MRStruct.Data = sum(MRStruct.Data,5);
if(isfield(MRStruct,'NoiseData') && numel(MRStruct.NoiseData) > 1)
    MRStruct.NoiseData = MRStruct.NoiseData .* AdditionalOut.CoilWeightMap;
    MRStruct.NoiseData = sum(MRStruct.NoiseData,5);    
end


%% Scale Data

if(UniformSignal_flag)
    AdditionalOut.Scaling = AdditionalOut.Mask ./ sum(abs(AdditionalOut.CoilWeightMap).^2,5);
else
   AdditionalOut.Scaling = AdditionalOut.Mask ./ sqrt(sum(abs(AdditionalOut.CoilWeightMap).^2,5));   
end
AdditionalOut.Scaling(isinf(AdditionalOut.Scaling) | isnan(AdditionalOut.Scaling)) = 0;

MRStruct.Data = MRStruct.Data .* AdditionalOut.Scaling;
MRStruct.Data(isinf(MRStruct.Data)) = 0;
if(isfield(MRStruct,'NoiseData') && numel(MRStruct.NoiseData) > 1)
    MRStruct.NoiseData = MRStruct.NoiseData .* AdditionalOut.Scaling;
    MRStruct.NoiseData(isinf(MRStruct.NoiseData)) = 0;
end

%% Postparations

MRStruct.RecoPar.DataSize(5) = 1;
MRStruct.RecoPar.total_channel_no_reco = 1;

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);




