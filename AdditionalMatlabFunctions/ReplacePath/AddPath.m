function [OrigPaths,NewPaths] = AddPath(PathToAdd,Verbosity)
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
	Verbosity = 1;
end


% 0.3 Definitions
    
OrigPaths = path;


%% Get Paths that should be added

PathToAdd2 = genpath(PathToAdd);
PathToAdd2 = regexp(PathToAdd2,':','split');

% Remove git folders
test = contains(PathToAdd2,'.git');
PathToAdd2 = PathToAdd2(~test); PathToAdd2(cellfun(@isempty,PathToAdd2)) = [];



%% Add Path

blubb = strcat(PathToAdd2,':'); blubb = strcat(blubb{:}); blubb(end) = [];
addpath(blubb)


%% Create 2nd Output
NewPaths = path;



%% 2. Postparations







