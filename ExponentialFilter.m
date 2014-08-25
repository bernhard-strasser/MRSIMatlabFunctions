function [OutArray,exp_filter_mat] = ExponentialFilter(InArray,dwelltime,exp_filter_Hz,ApplyAlongDim)
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
if(~exist('dwelltime','var'))
    dwelltime = 1;
end
if(~exist('exp_filter_Hz','var'))
    exp_filter_Hz = 3;
end
if(~exist('ApplyAlongDim','var'))
    ApplyAlongDim = numel(size(InArray));
end


% Define vecSize
vecSize = size(InArray,ApplyAlongDim);





%% 1. Compute exponential Time-Domain Filter

dwelltime_in_s = dwelltime/1000000000;

t= 0:dwelltime_in_s:dwelltime_in_s*(vecSize-1);
exp_filter_funct = exp(-exp_filter_Hz*t);     %exp(-t/a) wobei "1/a" = "exp_filter" Linebroadening in Hz




%% 2. Replicate Filter to size of InArray

% Find out along which dimensions the filter has to be replicated in order to match the size of InArray
ReplicateAlongDims = setdiff(1:numel(size(InArray)),ApplyAlongDim);

% Compute the size to which the filter has to be replicated. The first element is always 1, because this is the dimension of the filter (e.g. [1024 1])
ReplicateToSize = ones([1 numel(size(InArray))]);
loopy2 = 1;
for loopy = ReplicateAlongDims
    loopy2 = loopy2 + 1;
    ReplicateToSize(loopy2) = size(InArray,loopy);
end


% Replicate the filter, e.g. [1024 1] -> [1024 64 64]
exp_filter_mat = repmat(transpose(exp_filter_funct),ReplicateToSize);


% shiftdim: [1024 64 64] -> [64 64 1024], reshape: [64 64 1024] -> [64 64 1 1024]
exp_filter_mat = shiftdim(exp_filter_mat, numel(size(InArray)) - ApplyAlongDim + 1);
exp_filter_mat = reshape(exp_filter_mat, size(InArray));







%% 3. Apply Hamming Filter


OutArray = exp_filter_mat.*InArray;


