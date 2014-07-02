function [OutArray,HammingFilter] = HammingFilter(OutArray,ApplyAlongDims,FilterWidth,RadialOrOuterProduct,InputIskSpace_flag)
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
% -         OutArray                    ...		Input array to which the filter should be applied. For memory reasons InArray = OutArray.
% -         ApplyAlongDims              ...		Along these dimensions the filter is applied. If this vector has two elements, a two dimensional 
%												Filter is applied. Otherwise, a 3d filter is used.
% -			FilterWidth					...		Same as in spectroscopy sequences. Filter Width of 100 (%) means normal hamming filter,
%												filter width of n % means filter is only applied on n % (n/2 % on left, n/2 % on right) of
%												the data, the rest of the data is untouched (filter is set to 1 there).
% -         RadialOrOuterProduct        ...		Input: 'Radial' or 'OuterProduct'. Default: 'Radial'. There are two different methods to
%												create an n-dimensional filter from an 1-d filter: The OuterProduct is nothing more
%												than applying the 1d-filter in each dimension consecutively, 
%												i.e. Filter_nD(x1,x2,...,xn) = Filter_1D(x1)*Filter_1D(x2)*...*Filter_1D(xn)
%												Radial in contrary does Filter(x1,x2,...,xn) = Filter_1D(||x - M||), where x = (x1,x2,...,xn), 
%												and M is something like a center (around which point the filter should be applied, e.g.
%												the k-Space center). (Note: In fact it is a little more complicated, because if we want to
%												filter a 64x32 matrix, we have actually a elliptical filter, not radial...)
% -         InputIskSpace_flag          ...		If it is 0, the image gets Fourier transformed to k-space before applying the filter, 
%												and transformed back to image domain afterwards
%
% Output:
% -         OutArray                    ...		The filtered/masked output array
% -         HammingFilter               ...		The values of the Hamming filter in k-Space.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat_1_0

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


%OutArray = InArray; 
if(~exist('OutArray','var'))
	fprintf('\nHamming Filtering could not be performed.\nMore input needed.')
	return
end
if(~exist('ApplyAlongDims','var'))
	fprintf('\nWARNING: No dimensions specified along which Hamming Filter should be applied.\nFiltering along dim 1 with size %d',size(OutArray,1))
	ApplyAlongDims = 1;
end
if(~exist('FilterWidth','var'))
	FilterWidth = 100;
end
if( ~exist('RadialOrOuterProduct','var') || (~strcmpi(RadialOrOuterProduct,'Radial') && ~strcmpi(RadialOrOuterProduct,'OuterProduct')) )
	RadialOrOuterProduct = 'Radial';
end
if(~exist('InputIskSpace_flag','var'))
	InputIskSpace_flag = true;
end

Size_OutArray = size(OutArray);



%% 1. FFT to k-space

if(~InputIskSpace_flag)
    for hamming_dim = ApplyAlongDims
        OutArray = ifftshift(OutArray,hamming_dim);
        OutArray = ifft(OutArray,[],hamming_dim);
        OutArray = fftshift(OutArray,hamming_dim);
    end
end




%% 2. Compute Hamming Filter




if(strcmpi(RadialOrOuterProduct,'Radial'))

	HammingFilter = NaN(Size_OutArray(ApplyAlongDims));
	FilterSize = ceil(FilterWidth/100 * Size_OutArray(ApplyAlongDims));
	OnesSize = Size_OutArray(ApplyAlongDims) - FilterSize;
	OnesRadius = OnesSize/2;
	NormFactor = sqrt(sum(  (OnesRadius./(sqrt(numel(ApplyAlongDims))*FilterSize)).^2  ));
	
	
	
	% Loop over all Array entries... This is slow, but was easiest to implement now...
	% Loop Preps
	LinIndex = 1:prod(Size_OutArray(ApplyAlongDims));
	MatInd = cell( 1, numel(ApplyAlongDims) );
	for LinLoopy = LinIndex
		% Get from the linear index to the x,y,z,... in order to calculate the distance to k-Space center
		[MatInd{:}] = ind2sub(size(HammingFilter),LinLoopy);
		MatInd2 = [MatInd{:}] - (size(HammingFilter)+1)/2;
		
		% Dont process anything that is outside of the ellipsoid (consider ellipse equation!)
		if(sum((MatInd2 ./ (size(HammingFilter)/2)).^2) > 1 )
			continue;
		end
		
		if(	sum((MatInd2 ./ OnesRadius).^2) <= 1 )
			HammingFilter(LinLoopy) = 1;
		else
			HammingFilter(LinLoopy) = sum((MatInd2./FilterSize).^2);
		end
		
	end

	HammingFilter(HammingFilter ~= 1) = 0.54 + 0.46*cos(2*pi*(sqrt(HammingFilter(HammingFilter ~= 1)) - NormFactor));		

	
	% Extrapolate the filter to the edges
	HammingSize = size(HammingFilter);
	for hamming_dim = 1:numel(ApplyAlongDims)
		ZerosSize = size(HammingFilter); ZerosSize(hamming_dim) = 1;
		HammingFilter = cat(hamming_dim,zeros(ZerosSize),HammingFilter,zeros(ZerosSize));
	end
	HammingFilter = inpaint_nans(HammingFilter,4);
	HammingFilter = reshape(HammingFilter,HammingSize+2); %#ok
	% Crop out the good data
	CropString = '2:end-1,'; CropString = repmat(CropString,[1 numel(ApplyAlongDims)]); CropString(end) = [];
	HammingFilter = eval(['HammingFilter(' CropString ');']);
	
	
	% Replicate To size of OutArray
	HammingFilter = myrepmat(HammingFilter,Size_OutArray);
	
	
	
else
	
	HammingFilter = zeros(size(OutArray));	
	for hamming_dim = ApplyAlongDims                                            % Compute Hamming filter in each dimension seperately


		RepmatToSizeOfMatrixIn = size(OutArray);                                % The Hamming-filter must have the same size as the OutArray
		RepmatToSizeOfMatrixIn(hamming_dim) = 1;                                % Do not repmat in that dimension, in which hamming filtering is performed.



		%calculate Hamming filter
		n = size(OutArray,hamming_dim);
		hamming_1D = hamming(  ceil(FilterWidth/100*n)  );
		hamming_1D = cat(1, hamming_1D(1:ceil(ceil(FilterWidth/100*n)/2)), ones([n-numel(hamming_1D) 1]), hamming_1D(ceil(ceil(FilterWidth/100*n)/2)+1:end));


		if(hamming_dim == 1)
			hamming_1D_LeadingOnes = hamming_1D;
		else
			reshape_to = horzcat(ones([1 hamming_dim-1]), numel(hamming_1D));       % This creates e.g. a vector [1 1 1 1 64]
			hamming_1D_LeadingOnes = reshape(hamming_1D,reshape_to);                % Reshapes hamming filter to above size
		end

		HammingFilter = repmat(hamming_1D_LeadingOnes, RepmatToSizeOfMatrixIn) ...  % Replicate 1d-Hamming to the matrix size, 
						.* HammingFilter;                                           % Multiply this array with the previous calculated array
																					% but now the hamming-variation is in another dimension
	end

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


