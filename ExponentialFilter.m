function [OutArray,exp_filter_funct] = ExponentialFilter(InArray,ApplyAlongDims,exp_filter_Hz, dwelltime)
%
% EllipticalFilter_1_0 Apply an elliptical filter to k-space data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function masks the data in k-space, so that k-space values outside of an ellipsoid-shaped mask are set to zero. The mask can be a
% 3d-ellipsoid, or an 2d-ellipse. The equation for the mask is
% mask = {(x,y,z) E R³ | (x/a)² + (y/b)² + (z/c)² <= R²}
% a, b, c, and R can be chosen by the user.
%
%
% [A,B] = read_csi_dat_1_10(inputvar1,inputvar2)
%
% Input: 
% -         InArray                     ...    Input array to which the filter should be applied
% -         ApplyAlongDims              ...    Along these dimensions the filter is applied. If this vector has two elements, a two dimensional 
%                                              Filter is applied. Otherwise, a 3d filter is used.
% -         EllipsoidCoefficients       ...    The values for [a b c R], which determine the shape and size of the ellipsoid. For two dimensional
%                                              Filter, set c = 1;
% -         InputIskSpace_flag          ...    If it is 0, the image gets Fourier transformed to k-space before applying the filter, 
%                                              and transformed back to image domain afterwards
%
% Output:
% -         OutArray                    ...     The filtered/masked output array
% -         mask                        ...     The mask of the filter
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat_1_0

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


OutArray = InArray; 

%% 1. run Parameter.m

run ./tmp/Parameter.m;


%% 2. Compute exponential Time-Domain Filter

dwelltime_in_s = dwelltime/1000000000;

t= 0:dwelltime_in_s:dwelltime_in_s*(vecSize-1);
exp_filter_funct = exp(-exp_filter_Hz*t);     %exp(-t/a) wobei "1/a" = "exp_filter" Linebroadening in Hz


%% 3. Apply Hamming Filter

OutArray = (permute(repmat(exp_filter_funct', [1 ROW COL SLC]),[2 3 4 1])).*OutArray;         %apply exponential filter


