function [MRStruct] = op_MaskMRData(MRStruct, Settings, Mask)
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


%% Prep
if(~exist('Settings','var'))
   Settings = struct; 
end
if(~isfield(MRStruct,'Par'))
    MRStruct.Par = struct;
    if(isfield(MRStruct,'Data_file'))
        MRStruct.Par = read_ascconv(MRStruct.Data_file); 
    end
end
if(~isfield(MRStruct,'RecoPar'))
    MRStruct.RecoPar = MRStruct.Par;
end
if(~isfield(MRStruct.RecoPar,'DataSize'))
    MRStruct.RecoPar.DataSize = size(MRStruct.Data);
end


%% Find Correct Mask

if(~exist('Mask','var'))
    if(isfield(MRStruct,'BrainMask'))
        Mask = MRStruct.BrainMask;
    elseif(isfield(MRStruct,'Mask'))
        Mask = MRStruct.Mask;
    else
        return
    end
end


%% Mask Data

MRStruct.Data = MRStruct.Data .* Mask;
if(isfield(MRStruct,'NoiseData'))
    MRStruct.NoiseData = MRStruct.NoiseData .* Mask;    
end


%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);



