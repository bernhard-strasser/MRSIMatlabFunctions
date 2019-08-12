function [OrigPaths,NewPaths] = ReplacePath(PathToAdd,CellSearchFor)
%
% ReplacePath search for occurances of certain strings in a cell string
%
% This function was written by Bernhard Strasser, [month] [year].
%
%
% You can give two Cell Arrays of strings, CellToSearch and CellSearchFor to the function.
% The function will then tell you, where ANY of the strings in CellSearchFor matches
% CellToSearch.
%
% E.g. CellToSearch = {'A','B','BCDEFA'}; CellSearchFor = {'A','G'}
% Output: [1] [] [6], because A occurs in first and third cell of CellToSearch.
% 
%
%
% [A,B] = read_csi_dat_1_10(inputvar1,inputvar2)
%
% Input: 
% -         CellToSearch                   ...    This is the first input
% -         CellSearchFor                   ...    And this the second
%
% Output:
% -         Matches                           ...     This is output A
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations


% 0.3 Definitions
    
OrigPaths = path;



%% Get Existing Paths

OrigPaths = regexp(OrigPaths,':','split');
test = contains(OrigPaths,'.git');
OrigPaths = OrigPaths(~test);
% CurPath_LastFold = regexprep(OrigPaths,'.*/','');


%% Get Paths that should be added

PathToAdd2 = genpath(PathToAdd);
PathToAdd2 = regexp(PathToAdd2,':','split');

% Remove git folders
test = contains(PathToAdd2,'.git');
PathToAdd2 = PathToAdd2(~test); PathToAdd2(cellfun(@isempty,PathToAdd2)) = [];

% % Remove exact matches (folders that should be replaced
% ExactMatchesInd = regexp_SearchCellInCell(PathToAdd2,strcat(OrigPaths,'$'));
% PathToAdd2(cell2mat(ExactMatchesInd)) = [];

% Remove the path-part PathToAdd (e.g. if PathToAdd is '/sub1/sub2' and it has the subfolder '/sub1/sub2/sub3/sub4', only search '/sub3/sub4' in the original paths
% that should be replaced. If we only used /sub4 we would find probably too many replacements)
PathToAdd2_LastFold = PathToAdd2;
PathToAdd2_LastFold(2:end) = regexprep(PathToAdd2(2:end),PathToAdd,''); PathToAdd2_LastFold(cellfun(@isempty,PathToAdd2_LastFold)) = [];

% If we passed only one folder for replacing, always add this
if(isempty(PathToAdd2_LastFold))
    PathToAdd2_LastFold = regexprep(PathToAdd,'.*/','/');
end

% Add $, so that regexp is forced to match only end of path
% (otherwise e.g. /sub1/sub2 will be found in /sub1/sub2/sub3)
PathToAdd2_LastFold = strcat(PathToAdd2_LastFold,'$');



%% Compare Paths, and Replace

ReplaceInd = regexp_SearchCellInCell(OrigPaths,PathToAdd2_LastFold);

% Remove paths
for ii = 1:numel(ReplaceInd)
    
    CurRm = [];
    for jj = 1:numel(ReplaceInd{ii})
        CurRm = strcat(CurRm,OrigPaths{ReplaceInd{ii}(jj)},':');
    end
    CurRm = CurRm(1:end-1);
        
    
    if(~isempty(CurRm))
        if(strcmp(CurRm,PathToAdd2{ii}))       % dont text things like 'Replacing /sub1/sub2 with /sub1/sub2'
            continue
        end
        fprintf('\nReplacing path: %s with %s.',CurRm,PathToAdd2{ii})
        rmpath(CurRm)
    else
        fprintf('\nPath %s not found on current search-paths. Adding path %s without replacement.',PathToAdd2_LastFold{ii}(1:end-1),PathToAdd2{ii})
    end
    addpath(PathToAdd2{ii})
    
end

fprintf('\n')









%% old
% Pathlis2 = cell([1 2*numel(RmPaths)]);
% Pathlis2(1:2:end) = RmPaths;
% Pathlis2(2:2:end) = {':'};
% RmPaths2 = [Pathlis2{:}];
% if(~isempty(RmPaths2))
%     rmpath(RmPaths2);
% end
% 
% 
% % Add Paths
% Pathlis2 = cell([1 2*numel(PathToAdd2)]);
% Pathlis2(1:2:end) = PathToAdd2;
% Pathlis2(2:2:end) = {':'};
% PathsToAdd3 = [Pathlis2{:}];
% if(~isempty(PathsToAdd3))
%     for ii = 1:numel(PathToAdd2)
%        fprintf('\nReplacing paths: %s with %s.',RmPaths{ii},PathToAdd2{ii}) 
%     end
%     
%     
%     addpath(PathsToAdd3);
% end




%%
NewPaths = path;



%% 2. Postparations







