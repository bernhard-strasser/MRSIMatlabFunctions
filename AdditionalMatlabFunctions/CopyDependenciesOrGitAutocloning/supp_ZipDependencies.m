function supp_ZipDependencies(FilesToZip,ZipPath)
%
% MakeAutoCloningFile Create bash script to git clone all necessary repos for certain MATLAB function
%
% This function was written by Bernhard Strasser, March 2014.
%
%
% This function has two functionalities, depending on the input:
% 
% - If GitRepo is the string "Copy", find all the file dependencies for "FilesToZip" and zip
%   those to the zip-file given by ZipPath
% - If GitRepo is a string with the path to a GitRepo folder, find all the dependencies for "FilesToZip",
%   and search if a Git-Repo is available for all those dependencies in the subdirs of "GitRepo". If so,
%   write a git clone command in the bash script file "ZipPath". If it was not found, write a 
%   bash copy (cp) command to that file.
%
%
% hadamard_decoding_x(InArray,ApplyAlongDim)
%
% Input: 
% -     FilesToZip                  ...    The file for which the dependencies should be cloned/copied.
% -     ZipPath                    ...    The path to a .zip file to which all the dependencies are copied.
%
% Output: None (a zip-file or a .sh file is created)
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 




%% 0. Declarations, Preparations, Definitions


% Create full path of FilesToZip if necessary

CurDir = pwd;

% Make to Cell if it is not
if(~iscell(FilesToZip))
    FilesToZip2{1} = FilesToZip; FilesToZip = FilesToZip2; clear FilesToZip2
end

% open ZipPath
if(exist(ZipPath,'file'))
	fprintf('\nProblem: File\n%s\nexists.\nOverwrite?',ZipPath)
	Overwrite = input('[y][n].\n','s');
	if(strcmpi(Overwrite,'n'))
		return;
	end
end

% mkdir if path to ZipPath does not exists
PathToAutoCloningFile = regexprep(ZipPath,'/[^/]*$','');
if(~exist(PathToAutoCloningFile,'dir'))
	mkdir(PathToAutoCloningFile);
end


%% 1.

FileDeps = [];
for ii = 1:numel(FilesToZip)
    
    CurFileToZip = which(FilesToZip{ii});
    if(numel(CurFileToZip) == 0)
        fprintf('\nProblem: Function or Script\n%s\ndoes not exist.\n',FilesToZip)
        continue;
    end

    % Try to cd to folder where file lies in
    CDDir = regexp(CurFileToZip,'.*/','match'); 
    try
        CDDir = CDDir{1};
        cd(CDDir)
    catch ME
    end

    CurFileDeps = getFileDependencies(CurFileToZip);
    FileDeps = unique(union(FileDeps,CurFileDeps));

end


%% Zip

for File = transpose(FileDeps)
    unix(['zip -j -q ' ZipPath ' '  File{:}]);	
end

%% Posparations

cd(CurDir)

