function FoundStringInFiles = supp_SearchDependencies(FilesToSearch,SearchString)
%
% MakeAutoCloningFile Create bash script to git clone all necessary repos for certain MATLAB function
%
% This function was written by Bernhard Strasser, March 2014.
%
%
% This function searches a matlab script/function and all its dependencies for the keyword SearchString:
% 
%
%
% FoundStringInFiles = supp_SearchDependencies(FilesToSearch,SearchString)
%
% Input: 
% -     FilesToSearch                   ...   The file for which the dependencies should be searched.
% -     SearchString                    ...   The keyword for which all dependencies should be searched.
%
% Output: 
% -		FoundStringInFiles				...	  List of files where SearchString was found.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 




%% 0. Declarations, Preparations, Definitions


% Create full path of FilesToSearch if necessary

CurDir = pwd;

% Make to Cell if it is not
if(~iscell(FilesToSearch))
    FilesToSearch2{1} = FilesToSearch; FilesToSearch = FilesToSearch2; clear FilesToSearch2
end


%% 1.

FileDeps = [];
for ii = 1:numel(FilesToSearch)
    
    CurFileToZip = which(FilesToSearch{ii});
    if(numel(CurFileToZip) == 0)
        fprintf('\nProblem: Function or Script\n%s\ndoes not exist.\n',FilesToSearch)
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


%% Search Files

fprintf('\nFound string %s in the following dependencies:',SearchString)
FoundLogical = false([1 numel(FileDeps)]);
for CurFileNo = 1:numel(FileDeps)
    CurText = fileread(FileDeps{CurFileNo});
    parts = regexp(CurText, SearchString,'ONCE');
    if(~isempty(parts))
        fprintf('\n%s',FileDeps{CurFileNo})
        FoundLogical(CurFileNo) = true;
    end
end
fprintf('\n')
FoundStringInFiles = FileDeps(FoundLogical);

%% Posparations

cd(CurDir)

