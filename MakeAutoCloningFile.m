function MakeAutoCloningFile(FileToClone,GitRepo,AutoCloningFile)
%
% MakeAutoCloningFile Create bash script to git clone all necessary repos for certain MATLAB function
%
% This function was written by Bernhard Strasser, March 2014.
%
%
% This function has two functionalities, depending on the input:
% 
% - If GitRepo is the string "Copy", find all the file dependencies for "FileToClone" and zip
%   those to the zip-file given by AutoCloningFile
% - If GitRepo is a string with the path to a GitRepo folder, find all the dependencies for "FileToClone",
%   and search if a Git-Repo is available for all those dependencies in the subdirs of "GitRepo". If so,
%   write a git clone command in the bash script file "AutoCloningFile". If it was not found, write a 
%   bash copy (cp) command to that file.
%
%
% hadamard_decoding_x(InArray,ApplyAlongDim)
%
% Input: 
% -     FileToClone                  ...    The file for which the dependencies should be cloned/copied.
% -     GitRepo		                 ...    Either the path to the Git-Repo where all the dependencies are
%                                           searched for and from which they are then cloned
%                                           or the string "Copy" if the files should be directly copied and zipped.
% -     AutoCloningFile              ...    Either the path to a .sh file with all the git clone commands,
%                                           or the path to a .zip file (if GitRepo = 'Copy') to which all the dependencies
%                                           are copied.
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


% Create full path of FileToClone if necessary
FileToClone2 = which(FileToClone);
if(numel(FileToClone2) == 0)
	fprintf('\nProblem: Function or Script\n%s\ndoes not exist.\n',FileToClone)
	return;
end
FileToClone = FileToClone2; clear FileToClone2;

% open AutoCloningFile
if(exist(AutoCloningFile,'file'))
	fprintf('\nProblem: File\n%s\nexists.\nOverwrite?',AutoCloningFile)
	Overwrite = input('[y][n].\n','s');
	if(strcmpi(Overwrite,'n'))
		return;
	end
end

CopyFlag = ~isempty(strfind(GitRepo,'Copy'));

% mkdir if path to AutoCloningFile does not exists
PathToAutoCloningFile = regexprep(AutoCloningFile,'/\w*\.\w*','');
if(~exist(PathToAutoCloningFile,'dir'))
	mkdir(PathToAutoCloningFile);
end


% Open file if it should not be zipped and make a list of all Git Repos
if(~CopyFlag)
	fid = fopen(AutoCloningFile,'w');
	if(fid == -1)
		fprintf('\nProblem: Could not open file\n%s\n',AutoCloningFile)
		return;
	end
	fprintf(fid,'CloneToPath=""\n');	

	GitRepoDirs = genpath(GitRepo); 
	GitRepoDirs = regexp(GitRepoDirs,':','split');
	% Destroy all those entries that do not end with .git
	GitRepoDirs_unnecc = cellfun(@isempty,regexp(GitRepoDirs,'.git(?!/)'));	% Search for those that have .git but no following /. Get those  that are empty (i.e. those that have a following /).
	GitRepoDirs = GitRepoDirs(~GitRepoDirs_unnecc);
	
	
	
end


%% 1.

FileDeps = getFileDependencies(FileToClone);

for File = transpose(FileDeps)
	
	
	if(CopyFlag)
		unix(['zip -j -q ' AutoCloningFile ' '  File{:}]);
	else
		% Search for such a function in the GitRepo or copy it
		LastForwardSlash = strfind(File{:},'/'); 
		PreviousLastForwardSlash = LastForwardSlash(end-1);
		LastForwardSlash = LastForwardSlash(end);
		SearchForGitRepo = File{:}(LastForwardSlash+1:end-2);
		FoundGitRepo_Log = ~cellfun(@isempty,regexp(GitRepoDirs,SearchForGitRepo));
		
		% Write a clone command or a copy command, depending on if the GitRepo found
		if(sum(FoundGitRepo_Log) > 0)
			fprintf(fid,'git clone %s $CloneToPath/%s\n',GitRepoDirs{FoundGitRepo_Log},SearchForGitRepo);
		else
			Subdir = File{:}(PreviousLastForwardSlash+1:LastForwardSlash);
			if(isempty(regexpi(SearchForGitRepo,Subdir(1:end-1))))
				Subdir = '';
			end
			fprintf(fid,'mkdir $CloneToPath/%s\n',Subdir);						
			fprintf(fid,'cp %s $CloneToPath/%s%s.m\n',File{:},Subdir,SearchForGitRepo);			
		end
	end
	
end








%% Posparations

if(~CopyFlag)
	fclose(fid);
end


