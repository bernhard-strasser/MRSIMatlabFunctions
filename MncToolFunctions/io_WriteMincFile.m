function ErrorOccurred = io_WriteMincFile(MRStruct,MincFile,LikeMincFileOrMincPar)
% io_WriteMincFile Write an Array to minc file
%
% This function was written by Bernhard Strasser, September 2020.
%
%
% The function writes an array InArry to minc file MincFile using as header info either
% another minc-file or some parameters. It basically writes a raw file, and then converts this
% with rawtominc either with the -like option if a path to a minc file is given, or using the
% given parameters.
%
%
% ErrorOccurred = io_WriteMincFile(MRStruct,MincFile,LikeMincFileOrMincPar)
%
% Input: 
% -         MRStruct                       ...     Data to be written.
% -         MincFile                      ...     Path of minc file to which MRStruct should be written.
% -         LikeMincFileOrMincPar         ...     Path of other minc-file with same header, or Parameters containing
%                                                 infor for minc-header.
%
%
% Output:
% -        ErrorOccurred                          ...     1 if error occurred, otherwise 0.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None








%% 0. Preparations

if(~exist('LikeMincFileOrMincPar','var') || isempty(LikeMincFileOrMincPar))
    
    if(~isfield(MRStruct,'MncPars'))
        try
            MRStruct = op_CalcMincPars(MRStruct);
        catch
        end
    end


    if(isfield(MRStruct,'MncPars'))
        LikeMincFileOrMincPar = MRStruct.MncPars;
    else
  
        % Default Pars
        LikeMincFileOrMincPar.DimSize = size_MultiDims(MRStruct.Data,1:3);
        LikeMincFileOrMincPar.StepSize = [-220/LikeMincFileOrMincPar.DimSize(1) -220/LikeMincFileOrMincPar.DimSize(2) 220/LikeMincFileOrMincPar.DimSize(2)];  % Make z same as x and y
        LikeMincFileOrMincPar.ReadNormalVector = [1 0 0];
        LikeMincFileOrMincPar.PhaseNormalVector = [0 1 0];
        LikeMincFileOrMincPar.SliceNormalVector = [0 0 1];
        LikeMincFileOrMincPar.FoV_Read = 220;
        LikeMincFileOrMincPar.FoV_Phase = 220;
        LikeMincFileOrMincPar.FoV_Partition = abs(LikeMincFileOrMincPar.StepSize(1)) * LikeMincFileOrMincPar.DimSize(3);
        
        
        Pos_Minc = [0, 0, 0];      % Omit RotMat
        FoVHalf = [LikeMincFileOrMincPar.FoV_Read(1)/2 LikeMincFileOrMincPar.FoV_Phase(1)/2 -LikeMincFileOrMincPar.FoV_Partition(1)/2];  
        Pos_Minc = Pos_Minc + FoVHalf;
    
        LikeMincFileOrMincPar.POS_X_FirstVoxel = Pos_Minc(1) + abs(LikeMincFileOrMincPar.StepSize(1))/2;         % Be aware that StepRead and StepPhase are reversed and thus the sum is effectively a subtraction.
        LikeMincFileOrMincPar.POS_Y_FirstVoxel = Pos_Minc(2) + abs(LikeMincFileOrMincPar.StepSize(2))/2;
        LikeMincFileOrMincPar.POS_Z_FirstVoxel = Pos_Minc(3);
        if(LikeMincFileOrMincPar.DimSize(3) > 4)
            LikeMincFileOrMincPar.POS_Z_FirstVoxel = LikeMincFileOrMincPar.POS_Z_FirstVoxel + abs(LikeMincFileOrMincPar.StepSize(3))/2;
        end

    end

    
end
   
if(~isstruct(LikeMincFileOrMincPar))
    Info = dir(LikeMincFileOrMincPar);
    if(numel(Info) < 1 || Info.bytes == 0)
        fprintf('\nError in WriteMncFile: Could not find ''LikeMincFileOrMincPar %s''\n',LikeMincFileOrMincPar)
        ErrorOccurred = 1;
        return;
    end
end

% Unfortunately, this doesn't work. Minc needs to be setup before starting matlab! Only then it works!
% if(exist('MncPath','var'))
%     % If it is already loaded, you can still use this function, even if the following command gives an error
%     [tmp, tmpout] = unix(['. ' MncPath]);    
% end



%% Get info from LikeMincFileOrMincPar to reshape/reorder MRStruct

% NEEDS TO BE DONE STILL! THIS IS JUST WORKING FOR ONE CASE NOW (order: z,y,x)
% MRStruct = permute(MRStruct,[3 2 1]);


%% 1. Write Raw File


% Generate random unique id:
temp =  java.util.UUID.randomUUID;
myuuid = char(temp.toString);
RawFile = regexprep(MincFile,'.mnc',['_' myuuid '.raw']);
write_RawFiles(MRStruct.Data,RawFile);



%% 2. Convert Raw File to Mnc

if(~isstruct(LikeMincFileOrMincPar))
    [ErrorOccurred] = unix(['rawtominc -input ' RawFile ' -like ' LikeMincFileOrMincPar ' -float ' MincFile]);
    
else
    bashstring = ['rawtominc -clobber -float ' MincFile ' -input ' RawFile];
    bashstring = sprintf('%s -xstep %8.6f -ystep %8.6f -zstep %8.6f',bashstring, LikeMincFileOrMincPar.StepSize(1), LikeMincFileOrMincPar.StepSize(2),LikeMincFileOrMincPar.StepSize(3));    
    bashstring = sprintf('%s -xstart %8.6f -ystart %8.6f -zstart %8.6f',bashstring, LikeMincFileOrMincPar.POS_X_FirstVoxel,  LikeMincFileOrMincPar.POS_Y_FirstVoxel,  LikeMincFileOrMincPar.POS_Z_FirstVoxel);        
    bashstring = sprintf('%s -xdircos %8.6f %8.6f %8.6f -ydircos %8.6f %8.6f %8.6f -zdircos %8.6f %8.6f %8.6f',bashstring,LikeMincFileOrMincPar.ReadNormalVector,LikeMincFileOrMincPar.PhaseNormalVector,LikeMincFileOrMincPar.SliceNormalVector );            
    bashstring = sprintf('%s %d %d %d',bashstring,LikeMincFileOrMincPar.DimSize(3),LikeMincFileOrMincPar.DimSize(2),LikeMincFileOrMincPar.DimSize(1));     
    [ErrorOccurred] = unix(bashstring);

end






%% 4. Delete Raw File

delete(RawFile);



