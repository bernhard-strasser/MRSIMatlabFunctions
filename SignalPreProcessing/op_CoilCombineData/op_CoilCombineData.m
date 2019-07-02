function [Input, AdditionalOut] = op_CoilCombineData(Input,CoilWeightMap,Settings)
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
% -         Input                    ...     The Input. Must have fields 'Data', and 'Par'.
% -         CoilWeightMap             ...     The Coil-weighting-map. Must have field 'Data'
% -         Settings          ...            Structure to specify how the coil combination should be performed.
%
% Output:
% -         Input                      ...     The Output (Input is overwritten). Will have field 'RecoPar' additionally.
% -         AdditionalOut                        ...  Variable for storing additional output (if asked for by user).   
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: For now, the coil index in CoilWeightMap is the first, but in Input it's the last. Change that later!
% Also, the function cannot handle 3D / multislice data for now. And the channel index is fixed for now.


% This function expects the input to be of form
% ??? 


%% 0. Preparations

if(~exist('Settings','var'))
    Settings = [];
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
if(~isfield(Input,'RecoPar'))
    Input.RecoPar = Input.Par;
end


%% Resize SensMap

if(ImResize_flag)
    AdditionalOut.CoilWeightMap = zeros([Input.RecoPar.DataSize(1:2) Input.RecoPar.DataSize(end)]);
    for cha = 1:size(CoilWeightMap.Data,1)
        AdditionalOut.CoilWeightMap(:,:,cha) = imresize(squeeze(CoilWeightMap.Data(cha,:,:)),Input.RecoPar.DataSize(1:2));
    end
else
    AdditionalOut.CoilWeightMap = ZerofillOrCutkSpace(CoilWeightMap.Data,[size(CoilWeightMap.Data,1) Input.RecoPar.DataSize(1:3)],1);
    AdditionalOut.CoilWeightMap = permute(AdditionalOut.CoilWeightMap,[2 3 1]);
end
AdditionalOut.CoilWeightMap = myrepmat(AdditionalOut.CoilWeightMap,Input.RecoPar.DataSize);
AdditionalOut.Mask(:,:) = imresize(squeeze(CoilWeightMap.Mask(1,:,:)),Input.RecoPar.DataSize(1:2),'nearest');       
   

%% Weight Data

Input.Data = Input.Data .* AdditionalOut.CoilWeightMap;
Input.Data = sum(Input.Data,5);
if(isfield(Input,'NoiseData'))
    Input.NoiseData = Input.NoiseData .* AdditionalOut.CoilWeightMap;
    Input.NoiseData = sum(Input.NoiseData,5);    
end


%% Scale Data

if(UniformSignal_flag)
    AdditionalOut.Scaling = AdditionalOut.Mask ./ sum(abs(AdditionalOut.CoilWeightMap).^2,5);
else
   AdditionalOut.Scaling = AdditionalOut.Mask ./ sqrt(sum(abs(AdditionalOut.CoilWeightMap).^2,5));   
end
AdditionalOut.Scaling(isinf(AdditionalOut.Scaling) | isnan(AdditionalOut.Scaling)) = 0;

Input.Data = Input.Data .* AdditionalOut.Scaling;
Input.Data(isinf(Input.Data)) = 0;
if(isfield(Input,'NoiseData'))
    Input.NoiseData = Input.NoiseData .* AdditionalOut.Scaling;
    Input.NoiseData(isinf(Input.NoiseData)) = 0;
end

%% Adapt Sizes

Input.RecoPar.DataSize = Input.RecoPar.DataSize(1:end-1);
Input.RecoPar.total_channel_no_reco = 1;





