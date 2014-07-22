function SearchHistory(SearchString,CaseInsensitive_flag)
%
% SearchHistory Searches the file prefdir/history.m. 
%
% This function was written by Bernhard Strasser, July 2014.
%
%
% The function masks the data in k-space, so that k-space values outside of an ellipsoid-shaped mask are set to zero. The mask can be a
% 3d-ellipsoid, or an 2d-ellipse. The equation for the mask is
% mask = {(x,y,z) E R³ | (x/a)² + (y/b)² + (z/c)² <= R²}
% a, b, c, and R can be chosen by the user.
%
%
% [OutArray,mask] = EllipticalFilter_x_y(OutArray,ApplyAlongDims,EllipsoidCoefficients,PerformFFT_flag)
%
% Input: 
% -     OutArray                     ...    Input array to which the filter should be applied
% -     ApplyAlongDims               ...    Along these dimensions the filter is applied. If this vector has two elements, a two dimensional 
%                                          Filter is applied. Otherwise, a 3d filter is used.
% -     EllipsoidCoefficients        ...    The values for [a b c R], which determine the shape and size of the ellipsoid. For two dimensional
%                                          Filter, set c = 1;
% -     PerformFFT_flag              ...    If it is 0, the image gets Fourier transformed to k-space before applying the filter, 
%                                          and transformed back to image domain afterwards
%
% Output:
% -     IND                     ...     c
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

% 0.2 Declarations


% 0.3 Definitions
    



%% 1. Search

preffy = prefdir;

history_fid = fopen([preffy '/history.m'],'r');
HistoryString = fscanf(history_fid,'%c');

Daty = regexp(HistoryString,'%{1,2}--.{12,33}--%{1,2}','match');
HistoryString_DateSplit = regexp(HistoryString,'%{1,2}--.{15,33}--%{1,2}','split');
search_match = regexp(HistoryString_DateSplit,['\n[^\n]*' SearchString '[^\n]*\n'],'match');





%% 2. Print Out

if( ~iscell(search_match) || ~ismember(0,cellfun(@isempty,search_match)) )
	fprintf('\n%s not found in history.\n',SearchString)
	return;
end

for DateIndex = 1:numel(search_match)
	
	if(ismember(0,cellfun(@isempty,search_match{DateIndex})))
		fprintf('\n\n%s:\n', Daty{DateIndex})
		for LineIndex = 1:numel(search_match{DateIndex})
			fprintf('%s',search_match{DateIndex}{LineIndex})
		end
	end
	
end



%% 5. Postparations







