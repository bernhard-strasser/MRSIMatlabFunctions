function ErrorOccurred = io_WriteMincFile(InArray,MincFile,LikeMincFileOrMincPar)
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
% ErrorOccurred = io_WriteMincFile(InArray,MincFile,LikeMincFileOrMincPar)
%
% Input: 
% -         InArray                       ...     Data to be written.
% -         MincFile                      ...     Path of minc file to which InArray should be written.
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
    
end
   
if(isstruct(LikeMincFileOrMincPar))
    
else
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



%% Get info from LikeMincFileOrMincPar to reshape/reorder InArray

% NEEDS TO BE DONE STILL! THIS IS JUST WORKING FOR ONE CASE NOW (order: z,y,x)
% InArray = permute(InArray,[3 2 1]);


%% 1. Write Raw File


% Generate random unique id:
temp =  java.util.UUID.randomUUID;
myuuid = char(temp.toString);
RawFile = regexprep(MincFile,'.mnc',['_' myuuid '.raw']);
write_RawFiles(InArray,RawFile);



%% 2. Convert Raw File to Mnc


[ErrorOccurred,tmpout] = unix(['rawtominc -input ' RawFile ' -like ' LikeMincFileOrMincPar ' -float ' MincFile]);




%% 4. Delete Raw File

delete(RawFile);



