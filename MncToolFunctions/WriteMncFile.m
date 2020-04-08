function ErrorOccurred = WriteMncFile(InArray,MncFile,LikeMncFile)
% read_MncFiles Read a simple binary 'raw' file, e.g. created by mnc2raw
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function 
%
%
% image = read_MncFiles(MncFile,ROW,COL,SLC,precision)
%
% Input: 
% -         MncFile                       ...     Path of file.
%
%
% Output:
% -        image                          ...     The data which was read in.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None








%% 0. Preparations

Info = dir(LikeMncFile);
if(numel(Info) < 1 || Info.bytes == 0)
	fprintf('\nError in WriteMncFile: Could not fine ''LikeMncFile %s''\n',LikeMncFile)
	ErrorOccurred = 1;
	return;
end

% Unfortunately, this doesn't work. Minc needs to be setup before starting matlab! Only then it works!
% if(exist('MncPath','var'))
%     % If it is already loaded, you can still use this function, even if the following command gives an error
%     [tmp, tmpout] = unix(['. ' MncPath]);    
% end



%% Get info from LikeMncFile to reshape/reorder InArray

% NEEDS TO BE DONE STILL! THIS IS JUST WORKING FOR ONE CASE NOW (order: z,y,x)
InArray = permute(InArray,[3 2 1]);


%% 1. Write Raw File


% Generate random unique id:
temp =  java.util.UUID.randomUUID;
myuuid = char(temp.toString);
RawFile = regexprep(MncFile,'.mnc',['_' myuuid '.raw']);
write_RawFiles(InArray,RawFile);



%% 2. Convert Raw File to Mnc


[ErrorOccurred,tmpout] = unix(['rawtominc -input ' RawFile ' -like ' LikeMncFile ' -float ' MncFile]);




%% 4. Delete Raw File

delete(RawFile);



