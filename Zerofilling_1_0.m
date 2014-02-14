function [OutArray,mask] = Zerofilling_1_0(InArray,Zerofill_To,InputIskSpace_flag)
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
% -     mask                        ...     The mask of the filter
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat_1_0

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


% 0.1 Preparations

if(nargin < 2)
    OutArray = InArray;
    return
end
if(~exist('InputIskSpace_flag','var'))
   InputIskSpace_flag = true; 
end 


% 0.2 Declarations


% 0.3 Definitions
    
OutArray = InArray;

size_OutArray = size(OutArray);
AppendZeros = round(Zerofill_To - size_OutArray);
AppendZeros(AppendZeros < 0) = 0;
ApplyAlongDims = find(AppendZeros > 0);



 





%% 1. FFT to k-space

if(~InputIskSpace_flag)

    for filter_dim = ApplyAlongDims
        OutArray = ifftshift(OutArray,filter_dim);
        OutArray = ifft(OutArray,[],filter_dim);
        OutArray = fftshift(OutArray,filter_dim);
    end  

end




%% 2. Compute Mask


for dummy_dim = ApplyAlongDims

    
	AppendZeros_AtBeginning = zeros([size_OutArray(1:dummy_dim-1) ceil(AppendZeros(dummy_dim)/2) size_OutArray(dummy_dim+1:end)]);
	AppendZeros_AtEnd = zeros([size_OutArray(1:dummy_dim-1) floor(AppendZeros(dummy_dim)/2) size_OutArray(dummy_dim+1:end)]);
    OutArray = cat(dummy_dim,AppendZeros_AtBeginning,OutArray,AppendZeros_AtEnd);
    size_OutArray = size(OutArray);
    
end


%% 4. FFT to ImageSpace


if(~InputIskSpace_flag)
    
    for filter_dim = ApplyAlongDims
        OutArray = ifftshift(OutArray,filter_dim);
        OutArray = fft(OutArray,[],filter_dim);
        OutArray = fftshift(OutArray,filter_dim);
    end  
    
end




%% 5. Postparations







