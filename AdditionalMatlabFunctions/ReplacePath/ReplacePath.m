function [OrigPaths,NewPaths] = ReplacePath(PathToAdd,PathToRemove,Verbosity)
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


%% If PathToRemove is given, short-circuit

% In that case, just remove old path, and add new path.
if(exist('PathToRemove','var') && ~isempty(PathToRemove))
    RemovePath(PathToRemove,Verbosity);
    AddPath(PathToAdd,Verbosity);
    NewPaths = path;
    return;
end




%% Get Existing Functions

OrigPaths2 = regexp(OrigPaths,':','split');

% Remove git folders
test = contains(OrigPaths2,'.git');
OrigPaths2 = OrigPaths2(~test); OrigPaths2(cellfun(@isempty,OrigPaths2)) = [];

% Remove all matlabroot subdirectories. We don't want to touch them!!!
OrigPaths2 = OrigPaths2(cellfun(@isempty,regexp(OrigPaths2,matlabroot)));

% Find all .m files in OrigPaths2
LoadedFun = [];
LoadedFunPaths = [];
for ii = 1:numel(OrigPaths2)
    Dummy = dir(OrigPaths2{ii}); Dummy = {Dummy.name};
    Dummy = regexp(Dummy(3:end),'.*\.m','match','once');
    Dummy(cellfun(@isempty,Dummy)) = [];
    Dummy = regexprep(Dummy,'\.m','');
    LoadedFun = cat(2,LoadedFun,Dummy);
    LoadedFunPaths = cat(2,LoadedFunPaths,repmat(OrigPaths2(ii),size(Dummy)));
end


%% Get Functions that should be added

PathToAdd2 = genpath(PathToAdd);
PathToAdd2 = regexp(PathToAdd2,':','split');

% Remove git folders
test = contains(PathToAdd2,'.git');
PathToAdd2 = PathToAdd2(~test); PathToAdd2(cellfun(@isempty,PathToAdd2)) = [];

% Find all .m files in FunsToAddPaths
FunsToAdd = [];
FunsToAddPaths = [];
for ii = 1:numel(PathToAdd2)
    Dummy = dir(PathToAdd2{ii}); Dummy = {Dummy.name};
    Dummy = regexp(Dummy(3:end),'.*\.m','match','once');
    Dummy(cellfun(@isempty,Dummy)) = [];
    Dummy = regexprep(Dummy,'\.m','');
    FunsToAdd = cat(2,FunsToAdd,Dummy);
    FunsToAddPaths = cat(2,FunsToAddPaths,repmat(PathToAdd2(ii),size(Dummy)));
end




%% Compare Functions
% If identical, remove current function
% In any case, add new functions

RmPath = [];
AddPathlis = [];
for ii = 1:numel(FunsToAdd)
        
    FoundFun = find(strcmp(LoadedFun,FunsToAdd{ii}));
    if(numel(FoundFun) > 0)
        
%         if(strcmp([FunsToAddPaths{ii} '/' FunsToAdd{ii}],[LoadedFunPaths{ii} '/' FunsToAdd{ii}]))
%             continue
%         end
        
        if(Verbosity > 1)
            fprintf('\nReplacing function %s with %s/%s.',FunsToAdd{ii},FunsToAddPaths{ii},FunsToAdd{ii})
        end
        for jj = FoundFun
            RmPath = strcat(RmPath,LoadedFunPaths{jj},':');
        end
    else
        if(Verbosity > 0)
            fprintf('\nFunction %s not found on current search-paths. Adding path %s/%s without replacement.',FunsToAdd{ii},FunsToAddPaths{ii},FunsToAdd{ii})
        end
    end
    AddPathlis = strcat(AddPathlis,FunsToAddPaths{ii},':');
end
RmPath = RmPath(1:end-1);
AddPathlis = AddPathlis(1:end-1);


%% Remove and add Functions

if(~isempty(RmPath))
    rmpath(RmPath);
end
if(~isempty(AddPathlis))
    addpath(AddPathlis);
end


%%
NewPaths = path;



%% 2. Postparations







