function SearchHistory(SearchString,CaseInsensitive_flag)
%
% SearchHistory Searches the command history. 
%
% This function was written by Bernhard Strasser, July 2014.
%
%
% The function searches the command history (file prefdir/history.m) for SearchString.
%
%
% SearchHistory(SearchString,CaseInsensitive_flag)
%
% Input: 
% -     SearchString                   ...    String you want to search for.
% -     CaseInsensitive_flag           ...    If you want to search case insensitive, set this flag to true or 1.
%
% Output:
% None
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat_1_0

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


% 0.1 Preparations
if(nargin < 1)
	fprintf('\nWarning: Thou shalt specify a string to search for!\n')
	return
end
if(~exist('CaseInsensitive_flag','var'))
	CaseInsensitive_flag = false;
end
if(CaseInsensitive_flag)
	SearchFunction = 'regexpi';
else
	SearchFunction = 'regexp';
end
	

% 0.2 Declarations


% 0.3 Definitions
    



%% 1. Search

preffy = prefdir;

history_fid = fopen([preffy '/history.m'],'r');
HistoryString = fscanf(history_fid,'%c');

Daty = regexp(HistoryString,'%{1,2}--.{12,33}--%{1,2}','match');
HistoryString_DateSplit = regexp(HistoryString,'%{1,2}--.{15,33}--%{1,2}','split');
search_match = feval(SearchFunction,HistoryString_DateSplit,['\n[^\n]*' SearchString '[^\n]*\n'],'match');





%% 2. Print Out

if( ~iscell(search_match) || ~ismember(0,cellfun(@isempty,search_match)) )
	fprintf('\n%s not found in history.\n',SearchString)
	return;
end

for DateIndex = 1:numel(search_match)
	
	if(ismember(0,cellfun(@isempty,search_match{DateIndex})))
		fprintf('\n%s:\n', Daty{DateIndex})
		for LineIndex = 1:numel(search_match{DateIndex})
			fprintf('%s',search_match{DateIndex}{LineIndex})
		end
		fprintf('\n')
	end
	
end



%% 5. Postparations







