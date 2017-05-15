function out_array = squeeze_single_dim(in_array,squeeze_dim)
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
if(numel(size(in_array)) < squeeze_dim)
    display([char(10) sprintf('Error: Your input array has less than %d dimensions', squeeze_dim)])
    out_array = in_array;
    return
end

if(size(in_array,squeeze_dim) > 1)
    display([char(10) 'Warning in squeeze_single_dim: The dimension to squeeze is non-singleton. Take First element of dimension ' num2str(squeeze_dim) '.'])
	regy = cell([1 numel(size(in_array))]);
	regy(:) = {':,'};
	regy{squeeze_dim} = '1,';
	regy{end} = regy{end}(1);
	
	% if it is the last index, its done.
	if(squeeze_dim == numel(size(in_array)))
		out_array = eval(['in_array(' [regy{:}] ');']);
		return;
	else
		in_array = eval(['in_array(' [regy{:}] ');']);		
	end
	
end
     

% 0.2 Declarations


% 0.3 Definitions
    




%% 1. Squeeze

new_dimensions = setxor(1:numel(size(in_array)), squeeze_dim);  % numel(size(..)) gives the number of dimensions of in_array; so 1:numel... gives [1 2 3 ... dimensionality(in_array)]; setxor: removes squeeze_dim out of this vec
size_in_array = size(in_array);                                 % size of the array
new_size = size_in_array(new_dimensions);                       % gives a vector with the same size as size(in_array) but where the size_in_array(squeeze_dim) does not occur

out_array = reshape(in_array, new_size);


