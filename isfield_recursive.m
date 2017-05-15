function FieldExists_flag = isfield_recursive(InStruct,RecursiveFieldsToCheck)
%
% ExponentialFilter Apply an Exponential filter to time-domain signals (e.g. FIDs)
%
% This function was written by Wolfgang Bogner 2013, revised by Bernhard Strasser, October 2013.
%
%
% The function computes an exponential filter in Hertz
%
%
% [OutArray,exp_filter_funct] = ExponentialFilter(InArray,dwelltime,ApplyAlongDim,exp_filter_Hz)
%
% Input: 
% -         InArray                     ...    Input array to which the filter should be applied
% -         dwelltime                   ...    The dwelltime in [us], i.e. the time between two consecutive time points.
% -         ApplyAlongDim               ...    Along this dimension the filter is applied. 
% -         exp_filter_Hz               ...    The filter strength in Hz
%
% Output:
% -         OutArray                    ...     The filtered output array
% -         exp_filter_funct            ...     The filter
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy:

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


% Define variables if not defined






%% 1. Compute exponential Time-Domain Filter

FieldExists_flag = true;

RecursiveFieldsToCheck_cell = regexp(RecursiveFieldsToCheck,'\.','split');

CurStruct = InStruct;
for CurFieldNo = 1:numel(RecursiveFieldsToCheck_cell)
	CurField = RecursiveFieldsToCheck_cell{CurFieldNo};
	
	if(isfield(CurStruct,CurField))
		CurStruct = CurStruct.(CurField);
	else
		FieldExists_flag = false;
		break
	end
	
end


