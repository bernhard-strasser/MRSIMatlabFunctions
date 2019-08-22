function [OrigPaths,NewPaths] = RemovePath(PathToRemove,Verbosity)
%
% ReplacePath Replace matlab-functions by matlab-functions found in PathToAdd 
%
% This function was written by Bernhard Strasser, June 2019.
%
%
% The functions searches the PathToAdd for matlab functions and then searches the current functions for matches.
% If there is a match in the function names, the current one will be replaced by the new one. If it doesn't find a match, it will just add the new function.
%
% 
%
%
% [OrigPaths,NewPaths] = ReplacePath(PathToAdd,Verbosity)
%
% Input: 
% -         PathToAdd                   ...    The Path where matlab-functions are searched, and which will replace the old matlab-functions
% -         Verbosity                   ...    If 0, no text output will be printed, if 1 only warning-like outputs are printed, if 2 all outputs are printed
%
% Output:
% -         OrigPaths                   ...     All original matlab-functions that were replaced
% -         NewPaths                    ...     All new matlab-functions that replaced the old ones
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations

if(~exist('Verbosity','var'))
	Verbosity = 0;
end


% 0.3 Definitions
    
OrigPaths = path;



%% Get Existing Paths

OrigPaths = regexp(OrigPaths,':','split');
test = contains(OrigPaths,'.git');
OrigPaths = OrigPaths(~test);
% CurPath_LastFold = regexprep(OrigPaths,'.*/','');


%% Get Paths that should be removed

PathToRemove2 = genpath(PathToRemove);
PathToRemove2 = regexp(PathToRemove2,':','split');

% Remove git folders
test = contains(PathToRemove2,'.git');
PathToRemove2 = PathToRemove2(~test); PathToRemove2(cellfun(@isempty,PathToRemove2)) = [];



%% Compare Paths, and Remove


% Remove paths
RmPaths = [];
for ii = 1:numel(PathToRemove2)
    
    CurRm = strcmp(OrigPaths,PathToRemove2{ii});
    CurRm = OrigPaths(CurRm);
    
    if(~isempty(CurRm))
        RmPaths = strcat(RmPaths,CurRm{1},':');  
        
		if(Verbosity > 1)
        	fprintf('\nRemoving path: %s.',CurRm{1})
		end        
    else
		if(Verbosity > 0)
        	fprintf('\nWarning: Path %s not found on current path.',PathToRemove2{ii})
		end
    end

end

if(Verbosity > 0)
    fprintf('\n')
end
RmPaths = RmPaths(1:end-1);
rmpath(RmPaths)



%%
NewPaths = path;



%% 2. Postparations







