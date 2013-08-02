function [csi_out, Skip_kPoints_Spatial] = Skip_kPoints_1_1(csi_in, Keep_kPoints_ElementaryCell)
%
% EllipticalFilter_x_y Apply an elliptical filter to k-space data
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
% [OutArray,mask] = EllipticalFilter_x_y(InArray,ApplyAlongDims,EllipsoidCoefficients,InputIskSpace_flag)
%
% Input: 
% -     InArray                     ...    Input array to which the filter should be applied
% -     ApplyAlongDims              ...    Along these dimensions the filter is applied. If this vector has two elements, a two dimensional 
%                                          Filter is applied. Otherwise, a 3d filter is used.
% -     EllipsoidCoefficients       ...    The values for [a b c R], which determine the shape and size of the ellipsoid. For two dimensional
%                                          Filter, set c = 1;
% -     InputIskSpace_flag          ...    If it is 0, the image gets Fourier transformed to k-space before applying the filter, 
%                                          and transformed back to image domain afterwards
%
% Output:
% -     OutArray                    ...     The filtered/masked output array
% -     Skip_kPoints_Spatial        ...     The mask of the filter
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat_1_0

% Further remarks: 




%% 0. Preparations


% Compute the Ratios of the size(Input)/size(ElementaryCell)
SizeRatio = size(squeeze(csi_in(1,:,:,:,1))) ./ size(Keep_kPoints_ElementaryCell);


% Compute Elementary Cell of the k_points to skip
Skip_kPoints_ElementaryCell = ~Keep_kPoints_ElementaryCell;




%% 1. Replicate Elementary Cell

% Replicate in spatial dimensions.
Skip_kPoints_Spatial = repmat(Skip_kPoints_ElementaryCell, SizeRatio);

% Replicate in channel and time dimensions
Skip_kPoints = myrepmat_1_0(Skip_kPoints_Spatial, [size(csi_in,1) 1 1 1 size(csi_in,5)], 2);




%% 2. Set kPoints to Zero 

csi_out = csi_in;
csi_out(Skip_kPoints) = [];



