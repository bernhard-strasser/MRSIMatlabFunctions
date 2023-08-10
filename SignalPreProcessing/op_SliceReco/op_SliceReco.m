function [MRStruct] = op_SliceReco(MRStruct, Settings)
%
% op_PermuteMRData Perform Slice Reco of MR Data
%
% This function was written by Bernhard Strasser, Oct 2019.
%
%
% The function is meant to perform any slice-reconstruction, like reordering slices or Hadamard-decoding.
% For now, it just remomves the slice dimension and does nothing else.
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
% MRStruct:
% -         MRStruct               ...      The structure containing the MR data. Can have different fields.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: 


%% Prep
if(~exist('Settings','var'))
   Settings = struct; 
end
if(~isfield(Settings,'SliceIndex'))
    Settings.SliceIndex = 4;
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


%% Remove Slice or Partition Index

if(size(MRStruct.Data,Settings.SliceIndex) > 1)
    
    RemoveInd = 3;
else
    RemoveInd = Settings.SliceIndex;
end

MRStruct.Data = squeeze_single_dim(MRStruct.Data,RemoveInd); 
if(isfield(MRStruct,'NoiseData'))
    MRStruct.NoiseData = squeeze_single_dim(MRStruct.NoiseData,RemoveInd);
end


%% Update RecoPar
if(isfield(MRStruct.RecoPar,'dimnames'))
    MRStruct.RecoPar.dimnames(RemoveInd) = [];
end
MRStruct.RecoPar.DataSize = cat(2,size(MRStruct.Data),ones([1 6-ndims(MRStruct.Data)]));
MRStruct.RecoPar.nSLC = 1;


%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);



