function [MRStructCopyTo] = op_CopyAndMincresampleAllMasks(MRStructCopyFrom,MRStructCopyTo,Settings)
%
% op_CopyAndMincresampleAllMasks Apply minc-resample on all masks of an MRStruct
%
% This function was written by Bernhard Strasser, June 2019.
%
%
% The function 
%
%
% [MRStructCopyTo] = op_CopyAndMincresampleAllMasks(MRStructCopyFrom,MRStructCopyTo,Settings)
%
% Input: 
% -         MRStructCopyFrom        ...     
% -         MRStructCopyTo          ...     
% -         Settings                ...     
%
% Output:
% -         MRStructCopyTo                      ...   MRStruct to which the masks should be copied to,
% and to which the mincresample should be calculated.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: op_Mincresample, supp_UpdateRecoSteps





%% 0. Preparations

if(~exist('Settings','var'))
    Settings = struct;
end
if(~isfield(Settings,'BET_flag'))
    Settings.BET_flag = true;
end


%% Mincresample

MRStructCopyFrom = op_Mincresample(MRStructCopyFrom,MRStructCopyTo);


%% Copy Masks

MaskNames = fieldnames(MRStructCopyFrom.Masks);
for ii = 1:numel(Masks)
    MRStructCopyTo.Masks.(MaskNames{ii}) = MRStructCopyFrom.Masks.(MaskNames{ii});
end


%% Postparations


MRStructCopyTo = supp_UpdateRecoSteps(MRStructCopyTo,Settings);



