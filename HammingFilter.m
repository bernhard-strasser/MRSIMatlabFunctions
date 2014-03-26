function [OutArray,HammingFilter] = HammingFilter(OutArray,ApplyAlongDims,InputIskSpace_flag)
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
% -         OutArray                     ...   Input array to which the filter should be applied. For memory reasons InArray = OutArray.
% -         ApplyAlongDims              ...    Along these dimensions the filter is applied. If this vector has two elements, a two dimensional 
%                                              Filter is applied. Otherwise, a 3d filter is used.
% -         InputIskSpace_flag          ...    If it is 0, the image gets Fourier transformed to k-space before applying the filter, 
%                                              and transformed back to image domain afterwards
%
% Output:
% -         OutArray                    ...     The filtered/masked output array
% -         HammingFilter               ...     The values of the Hamming filter in k-Space.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat_1_0

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


%OutArray = InArray; 


%% 1. FFT to k-space

if(~InputIskSpace_flag)
    for hamming_dim = ApplyAlongDims
        OutArray = ifftshift(OutArray,hamming_dim);
        OutArray = ifft(OutArray,[],hamming_dim);
        OutArray = fftshift(OutArray,hamming_dim);
    end
end





%% 2. Compute Hamming Filter


HammingFilter = ones(size(OutArray));

for hamming_dim = ApplyAlongDims                                            % Compute Hamming filter in each dimension seperately
    
    
    RepmatToSizeOfMatrixIn = size(OutArray);                                % The Hamming-filter must have the same size as the OutArray
    RepmatToSizeOfMatrixIn(hamming_dim) = 1;                                % Do not repmat in that dimension, in which hamming filtering is performed.
    


    %calculate Hamming filter
    hamming_1D = (hamming(size(OutArray,hamming_dim)));

    if(hamming_dim == 1)
        hamming_1D_LeadingOnes = hamming_1D;
    else
        reshape_to = horzcat(ones([1 hamming_dim-1]), numel(hamming_1D));       % This creats e.g. a vector [1 1 1 1 64]
        hamming_1D_LeadingOnes = reshape(hamming_1D,reshape_to);                % Reshapes hamming filter to above size
    end
    
    HammingFilter = repmat(hamming_1D_LeadingOnes, RepmatToSizeOfMatrixIn) ...  % Replicate 1d-Hamming to the matrix size, 
                    .* HammingFilter;                                           % Multiply this array with the previous calculated array
                                                                                % but now the hamming-variation is in another dimension
end





%% 3. Apply Hamming Filter

OutArray = OutArray .* HammingFilter;






%% 4. FFT to Image Space

if(~InputIskSpace_flag)
    
    for hamming_dim = ApplyAlongDims
        OutArray = ifftshift(OutArray,hamming_dim);
        OutArray = fft(OutArray,[],hamming_dim);
        OutArray = fftshift(OutArray,hamming_dim);
    end
    
end


%% 5. Postparations


