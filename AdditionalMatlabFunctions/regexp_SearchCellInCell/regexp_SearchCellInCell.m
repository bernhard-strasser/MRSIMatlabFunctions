function [MatchInd,Matches2,varargout] = regexp_SearchCellInCell(CellToSearch,CellSearchFor,regexpi_flag,varargin)
%
% regexp_SearchCellInCell search for occurances of certain strings in a cell string
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

Matches = 0;

if(nargin < 2)
	fprintf('\nError: Too little input. Aborting.\n')
	return
end
if(~exist('regexpi_flag','var'))
	regexpi_flag = false;
end
if(regexpi_flag)
	regexp_fhandle = @regexpi;
else
	regexp_fhandle = @regexp;
end

% 0.3 Definitions
    


%% 1. To be on the safe side, replace all special characters with \[special character]

%CellSearchFor2 = regexprep(CellSearchFor,'\W','\\$0');
CellSearchFor2 = CellSearchFor;










%% 2. Run regexp on that.


MatchInd = cell([1 numel(CellSearchFor2)]);
Matches2 = zeros([1 numel(CellToSearch)]);
for i = 1:numel(CellSearchFor2)
	
	if(~isempty(varargin))
		if(nargout > 1)
			[Matches,varargout] = feval(regexp_fhandle,CellToSearch,CellSearchFor2{i},varargin);
		else
			[Matches] = feval(regexp_fhandle,CellToSearch,CellSearchFor2{i},varargin{:});		
		end
	else
		if(nargout > 1)
			[Matches,varargout] = feval(regexp_fhandle,CellToSearch,CellSearchFor2{i});
		else
			[Matches] = feval(regexp_fhandle,CellToSearch,CellSearchFor2{i});		
		end
    end
    Matches2 = double(~cellfun(@isempty,Matches)) + Matches2;
	MatchInd{i} = find(~cellfun(@isempty,Matches));
	
end


%% 2. Postparations

% fclose(fid)






