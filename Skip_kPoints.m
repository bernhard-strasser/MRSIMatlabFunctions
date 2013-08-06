function [OutData, Skip_kPoints_Spatial] = Skip_kPoints(InData, Keep_kPoints_ElementaryCell)
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
SizeRatio = size(squeeze(InData(1,:,:,:,1))) ./ size(Keep_kPoints_ElementaryCell);


% Compute Elementary Cell of the k_points to skip
Skip_kPoints_ElementaryCell = ~Keep_kPoints_ElementaryCell;




%% 1. Replicate Elementary Cell

% Replicate in spatial dimensions.
Skip_kPoints_Spatial = repmat(Skip_kPoints_ElementaryCell, floor(SizeRatio));

% Append beginning of pattern if elementary cell cannot be replicated to the size of InData
Skip_kPoints_Spatial = cat(1,Skip_kPoints_Spatial,Skip_kPoints_Spatial(1:size(InData,2)-size(Skip_kPoints_Spatial,1),:,:));
Skip_kPoints_Spatial = cat(2,Skip_kPoints_Spatial,Skip_kPoints_Spatial(:,1:size(InData,3)-size(Skip_kPoints_Spatial,2),:));
Skip_kPoints_Spatial = cat(3,Skip_kPoints_Spatial,Skip_kPoints_Spatial(:,:,1:size(InData,4)-size(Skip_kPoints_Spatial,3)));

%figure; imagesc(Skip_kPoints_Spatial)
% Append zeros if elementary cell cannot be replicated to the size of InData
% if(mod(SizeRatio(1),1) ~= 0)
%     Append_RowNo = size(InData,2) - size(Skip_kPoints_Spatial,1);
% 	Append_RowAtBeginning = ones([floor(Append_RowNo/2) size(Skip_kPoints_Spatial,2)]);
% 	Append_RowAtEnd = ones([ceil(Append_RowNo/2) size(Skip_kPoints_Spatial,2)]);
%     Skip_kPoints_Spatial = cat(1,Append_RowAtBeginning,Skip_kPoints_Spatial,Append_RowAtEnd);
% end
% if(mod(SizeRatio(2),1) ~= 0)
%     Append_ColNo = size(InData,3) - size(Skip_kPoints_Spatial,2);
% 	Append_ColAtBeginning = ones([size(Skip_kPoints_Spatial,1) floor(Append_ColNo/2)]);
% 	Append_ColAtEnd = ones([size(Skip_kPoints_Spatial,1) ceil(Append_ColNo/2)]);
%     Skip_kPoints_Spatial = cat(2,Append_ColAtBeginning,Skip_kPoints_Spatial,Append_ColAtEnd);
% end
Skip_kPoints_Spatial = logical(Skip_kPoints_Spatial);
%figure; imagesc(Skip_kPoints_Spatial)

% Replicate in channel and time dimensions
Skip_kPoints = repmat(Skip_kPoints_Spatial, [1 1 1 size(InData,5)]);



%% 2. Set kPoints to Zero 

OutData = zeros(size(InData));
for cha = 1:size(InData,1)
    csi_out_dummy = InData(cha,:,:,:,:);
    csi_out_dummy(Skip_kPoints) = 0;
    OutData(cha,:,:,:,:) = csi_out_dummy;
end


