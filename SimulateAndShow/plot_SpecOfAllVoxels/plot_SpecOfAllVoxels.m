function plot_SpecOfAllVoxels(MRStruct, Settings,Mask)
%
% op_PermuteMRData Permute MR Data
%
% This function was written by Bernhard Strasser, Oct 2019.
%
%
% The function just runs the "permute" function on MR data as specified by Settings.PermuteVec.
% It also takes care of things like updating the .Par.DataSize, and permuting the NoiseData, if available.
%
%
% [MRStruct] = op_PermuteMRData(MRStruct,Settings)
%
% Input: 
% -         MRStruct               ...      The structure containing the MR data. Can have different fields. For this function, it
%                                           it has to have field
%                                           * .Data
% -         Settings               ...      Structure with settings how data should be processed. Should have field:
%                                           * .PermuteVec: Data will be permuted to according to this vector.
%
% Output:
% -         MRStruct               ...      The structure containing the MR data. Can have different fields.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: 

%% Preparations

if(~exist('Settings','var'))
    Settings = struct;
end
if(~isfield(Settings,'UseThisInStructMask'))
    Settings.UseThisInStructMask = 'BrainMask';
end
if(isfield(Settings,'UseThisInStructMask') && ~exist('Mask','var') && isfield(MRStruct,(Settings.UseThisInStructMask)))
    Mask = MRStruct.(Settings.UseThisInStructMask);
end
if(~exist('Mask','var'))
    Mask = ones(size_MultiDims(MRStruct.Data,1:3));
end

%% 

MRStruct.Data = MRStruct.Data(:,:,:,:,1);
MRStruct.Data = MRStruct.Data .* Mask;
MRStruct.Data = abs(fftshift(fft(MRStruct.Data,[],4),4));
MRStruct.Data = permute(MRStruct.Data,[4 1 2 3]);
chemy = compute_chemshift_vector(MRStruct.RecoPar);
figure; plot(chemy,MRStruct.Data(:,:))
figure; plot(chemy,sum(MRStruct.Data(:,:),2))



