function [MRStructResampleSource] = op_Mincresample(MRStructResampleSource,MRStructResampleTarget,Settings)
%
% op_Mincresample Mincresample the data of source from the orientation of source to orientation of target
%
% This function was written by Bernhard Strasser, June 2019.
%
%
% The function can mincresample the data fields .Data, .NoiseData, .Masks.[...] of MRStructResampleSource
% from the orientation&position given by MRStructResampleSource.Par to the orientation&position
% given by MRStructResampleTarget
%
%
% [MRStructResampleSource] = op_Mincresample(MRStructResampleSource,MRStructResampleTarget,Settings)
%
% Input: 
% -         MRStructResampleSource      ...    Source MRStruct from whose data should be resampled 
% -         MRStructResampleTarget      ...    Target MRStruct to whose orientation the source should be resampled to
% -         Settings                    ...     
%
% Output:
% -         MRStructResampleSource      ...     Resampled data
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
if(~isfield(Settings,'tmp_dir'))
    Settings.tmp_dir = '/tmp';
end


%% Create MincFiles


if(~isfield(MRStructResampleSource,'MincPar'))
    MRStructResampleSource = op_CalcMincPars(MRStructResampleSource);
end

% Loop over all fields like .Data, .NoiseData, .Mask, .BrainMask etc. Can this be solved more elegantly instead of just looping over a fixed set of fields?
Fields = {'Data','NoiseData','Masks.Mask','Masks.BrainMask','Masks.LipMask'};
SourceMincFilePath = cell([1 numel(Fields)]); SourceResMincFilePath = SourceMincFilePath;
for ii = 1:numel(Fields)
    temp =  java.util.UUID.randomUUID; myuuid = char(temp.toString);
    SourceMincFilePath{ii} = [Settings.tmp_dir '/Source_' myuuid '.Minc'];
    SourceResMincFilePath{ii} = [Settings.tmp_dir '/Source_' myuuid '_Res.Minc'];
    CurData = MRStructResampleSource;
    CurField = Fields{ii};
    while(~isempty(regexp(CurField,'\.','ONCE')))
        Field1 = regexp(Fields{ii},'\.','split');
        CurData = CurData.(Field1{1});
        CurField = regexprep(CurField,[Field1{1} '.']);
    end
    io_WriteMincFile(SourceMincFilePath{ii},MRStructResampleSource.(Fields{ii}));
end

temp =  java.util.UUID.randomUUID; myuuid = char(temp.toString);
if(~isfield(MRStructResampleTarget,'MincPar'))
    MRStructResampleTarget = op_CalcMincPars(MRStructResampleTarget);
end
TargetMincFilePath = [Settings.tmp_dir '/Target_' myuuid '.Minc'];
io_WriteMincFile(TargetMincFilePath,MRStructResampleTarget);


%% Mincresample

for ii = 1:numel(Fields)
    RunText = sprintf('mincresample %s -like %s %s %s',Settings.MincresampleOptions,TargetMincFilePath,SourceMincFilePath{ii},SourceResMincFilePath{ii});
    unix(RunText);
end


%% Read Minc Files & Delete Them Again


for ii = 1:numel(Fields)
    delete(SourceMincFilePath{ii});
    delete(SourceResMincFilePath{ii});
end
delete(TargetMincFilePath{ii});

    



%% Postparations


MRStructResampleSource = supp_UpdateRecoSteps(MRStructResampleSource,Settings);



