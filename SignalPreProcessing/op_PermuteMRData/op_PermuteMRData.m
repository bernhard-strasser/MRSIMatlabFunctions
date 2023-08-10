function [MRStruct] = op_PermuteMRData(MRStruct, Settings)
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

if(~isstruct(MRStruct) || numel(fieldnames(MRStruct)) < 1)
    return;
end

if(~exist('Settings','var'))
   Settings = struct; 
else
    if(~isstruct(Settings))
        Bak = Settings; clear Settings;
        Settings.PermuteVec = Bak; clear Bak;
    end
end
if(~isfield(Settings,'PermuteVec'))
    Settings.PermuteVec = 1:numel(size(MRStruct.Data));
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


%% Permute Data

for fieldyy = {'Data','NoiseData','Mask','BrainMask','LipMask'}
    if(isfield(MRStruct,fieldyy{1}))
        if(iscell(MRStruct.(fieldyy{1})))
            for ii = 1:numel(MRStruct.(fieldyy{1}))
                MRStruct.(fieldyy{1}){ii} = permute(MRStruct.(fieldyy{1}){ii},Settings.PermuteVec);
                MRStruct.RecoPar.DataSize{ii} = size(MRStruct.Data{ii});
            end
        else
        MRStruct.(fieldyy{1}) = permute(MRStruct.(fieldyy{1}),Settings.PermuteVec);
            MRStruct.RecoPar.DataSize = size(MRStruct.Data);
    end    
end
end


%% Reorder dimnames
if(isfield(MRStruct.RecoPar,'dimnames'))
    MRStruct.RecoPar.dimnames = MRStruct.RecoPar.dimnames(Settings.PermuteVec);
end
if(isfield(MRStruct.RecoPar,'dimnames_small'))
    MRStruct.RecoPar.dimnames_small = MRStruct.RecoPar.dimnames_small(Settings.PermuteVec);
end




%% Postparations
MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);



