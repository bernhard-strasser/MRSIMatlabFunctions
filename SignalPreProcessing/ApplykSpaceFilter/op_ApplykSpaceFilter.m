function [MRStruct,kSpaceFilter] = op_ApplykSpaceFilter(MRStruct,Settings)
%
% kSpaceFilter Apply a filter to (k-space) data
%
% This function was written by Bernhard Strasser, July 2012 - July 2014.
%
%
% The function applies a filter to the input-data. You can specify along which dimensions this should be done (e.g. [2 4 5 6]). You can
% also specify the filter width in %, e.g. a Settings.FilterOption of 70 % leaves the inner 30% of the k-Space untouched and applies the filter only
% along the 70 % of the outer data. You can also specify how the multi-dimensional filter is created, i.e. by the OuterProduct, 
% kSpaceFilter2D(x,y) = kSpaceFilter1D(x) * kSpaceFilter1D(y), or by a Radial 'rotation' (e.g. creating a 1D-hamming filter along the x-axis and then rotating
% that around the z-axis to get a 2D-filter. If you tell the function that the Input is already in kSpace, no fft is performed before applying
% the filter.
%
%
% [MRStruct,kSpaceFilter] = op_ApplykSpaceFilter(MRStruct,Settings)
%
% Input: 
% -         MRStruct                    ...   Input/Output array to which the filter should be applied. For memory reasons InArray = MRStruct.
% -         Settings.ApplyAlongDims              ...   Along these dimensions the filter is applied. If this vector has two elements, a two dimensional 
%                                             Filter is applied. Otherwise, a 3d filter is used.
% -         Settings.FilterTypeName              ...   Name of filter that should be applied. Examples: 'hamming', 'hann', 'hanning', 'chebwin', 'gausswin', see 'help window'
% -         Settings.FilterOption                ...   In case of hamming filter: Same as in spectroscopy sequences. Filter Width of 100 (%) means normal hamming filter,
%                                             filter width of n % means filter is only applied on n % (n/2 % on left, n/2 % on right) of
%                                             the data, the rest of the data is untouched (filter is set to 1 there). 
%                                             See p74 (A.2-26) Spectroscopy Manual.
%                                             Non-Hamming filter: option that can be specified by chosen window.
% -         Settings.RadialOrOuterProduct       ...    Input: 'Radial' or 'OuterProduct'. Default: 'Radial'. There are two different methods to
%                                             create an n-dimensional filter from an 1-d filter: The OuterProduct is nothing more
%                                             than applying the 1d-filter in each dimension consecutively, 
%                                             i.e. Filter_nD(x1,x2,...,xn) = Filter_1D(x1)*Filter_1D(x2)*...*Filter_1D(xn)
%                                             Radial in contrary does Filter(x1,x2,...,xn) = Filter_1D(||x - M||), where x = (x1,x2,...,xn), 
%                                             and M is something like a center (around which point the filter should be applied, e.g.
%                                             the k-Space center). (Note: In fact it is a little more complicated, because if we want to
%                                             filter a 64x32 matrix, we have actually a elliptical filter, not radial...)
% -         Settings.InputIskSpace_flag         ...    If it is 0, the image gets Fourier transformed to k-space before applying the filter, 
%                                             and transformed back to image domain afterwards
%
% Output:
% -         MRStruct                   ...    The filtered/masked output array
% -         kSpaceFilter              ...    The values of the filter in k-Space.
%
%
% EXAMPLE: 
% - Test = op_ApplykSpaceFilter(MRStructIn,struct('FilterTypeName','chebwin','FilterOption',50,'InputIskSpace_flag',0,'ApplyAlongDims',[1 2]));
% - Test = op_ApplykSpaceFilter(MRStructIn,struct('FilterTypeName','hamming','FilterOption',100,'InputIskSpace_flag',0,'ApplyAlongDims',[1 2]));
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat_1_0

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


%MRStruct = InArray; 
if(~exist('MRStruct','var'))
	fprintf('\nFiltering could not be performed.\nMore input needed.')
	return
end
if(~isfield(Settings,'ApplyAlongDims'))
	fprintf('\nWARNING: No dimensions specified along which Filter should be applied.\nFiltering along dim 1 with size %d',size(MRStruct.Data,1))
	Settings.ApplyAlongDims = 1;
end
if(~isfield(Settings,'FilterTypeName'))
    Settings.FilterTypeName = 'hamming';
end
if(numel(Settings.FilterOption)==0)
    Settings.FilterOption = [];
end
if( ~isfield(Settings,'RadialOrOuterProduct') || (~strcmpi(Settings.RadialOrOuterProduct,'Radial') && ~strcmpi(Settings.RadialOrOuterProduct,'OuterProduct')) )
	Settings.RadialOrOuterProduct = 'OuterProduct';
end
if(~isfield(Settings,'InputIskSpace_flag'))
	Settings.InputIskSpace_flag = true;
end

Size_MRStruct = size(MRStruct.Data);



%% 1. FFT to k-space

if(~Settings.InputIskSpace_flag)
    for filter_dim = Settings.ApplyAlongDims
        MRStruct.Data = ifftshift(MRStruct.Data,filter_dim);
        MRStruct.Data = ifft(MRStruct.Data,[],filter_dim);
        MRStruct.Data = fftshift(MRStruct.Data,filter_dim);
    end
end




%% 2. Compute  Filter




if(strcmpi(Settings.RadialOrOuterProduct,'Radial'))
    if(strcmp(Settings.FilterTypeName,'hamming'))

        kSpaceFilter = NaN(Size_MRStruct(Settings.ApplyAlongDims));
        FilterSize = ceil(Settings.FilterOption/100 * Size_MRStruct(Settings.ApplyAlongDims));
        OnesSize = Size_MRStruct(Settings.ApplyAlongDims) - FilterSize;
        OnesRadius = OnesSize/2;
        NormFactor = sqrt(sum(  (OnesRadius./(sqrt(numel(Settings.ApplyAlongDims))*FilterSize)).^2  ));



        % Loop over all Array entries... This is slow, but was easiest to implement now...
        % Loop Preps
        LinIndex = 1:prod(Size_MRStruct(Settings.ApplyAlongDims));
        MatInd = cell( 1, numel(Settings.ApplyAlongDims) );
        for LinLoopy = LinIndex
            % Get from the linear index to the x,y,z,... in order to calculate the distance to k-Space center
            [MatInd{:}] = ind2sub(size(kSpaceFilter),LinLoopy);
            MatInd2 = [MatInd{:}] - (size(kSpaceFilter)+1)/2;

            % Dont process anything that is outside of the ellipsoid (consider ellipse equation!)
            if(sum((MatInd2 ./ (size(kSpaceFilter)/2)).^2) > 1 )
                continue;
            end

            if(	sum((MatInd2 ./ OnesRadius).^2) <= 1 )
                kSpaceFilter(LinLoopy) = 1;
            else
                kSpaceFilter(LinLoopy) = sum((MatInd2./FilterSize).^2);
            end

        end

        kSpaceFilter(kSpaceFilter ~= 1) = 0.54 + 0.46*cos(2*pi*(sqrt(kSpaceFilter(kSpaceFilter ~= 1)) - NormFactor));		


        % Extrapolate the filter to the edges
        kSpaceSize = size(kSpaceFilter);
        for filter_dim = 1:numel(Settings.ApplyAlongDims)
            ZerosSize = size(kSpaceFilter); ZerosSize(filter_dim) = 1;
            kSpaceFilter = cat(filter_dim,zeros(ZerosSize),kSpaceFilter,zeros(ZerosSize));
        end
        kSpaceFilter = inpaint_nans(kSpaceFilter,4);
        kSpaceFilter = reshape(kSpaceFilter,kSpaceSize+2); %#ok
        % Crop out the good data
        CropString = '2:end-1,'; CropString = repmat(CropString,[1 numel(Settings.ApplyAlongDims)]); CropString(end) = [];
        kSpaceFilter = eval(['kSpaceFilter(' CropString ');']);


        % Replicate To size of MRStruct
        kSpaceFilter = myrepmat(kSpaceFilter,Size_MRStruct);
    else
       	error('Option ''Radial'' currently only supported for Hamming filter.') 
    end
	
else
	
	kSpaceFilter = ones(size(MRStruct.Data));	
	for filter_dim = Settings.ApplyAlongDims                                            % Compute filter in each dimension seperately


		%calculate filter
		n = size(MRStruct.Data,filter_dim);
        
        if(strcmp(Settings.FilterTypeName,'hamming'))
            if(isfield(Settings,'FilterOption'))
                filter_1D = hamming(  ceil(Settings.FilterOption/100*n)  );
                filter_1D = cat(1, filter_1D(1:ceil(ceil(Settings.FilterOption/100*n)/2)), ones([n-numel(filter_1D) 1]), filter_1D(ceil(ceil(Settings.FilterOption/100*n)/2)+1:end));
            else
                filter_1D = hamming( n );
                filter_1D = cat(1, filter_1D(1:ceil(n/2)), ones([n-numel(filter_1D) 1]), filter_1D(ceil(n/2)+1:end));
            end
        else
            if(isfield(Settings,'FilterOption'))
                filter_1D = window(Settings.FilterTypeName,n,Settings.FilterOption);
            else
                 filter_1D = window(Settings.FilterTypeName,n);               
            end
        end



        % Define the dimension of Size_MRStruct which corresponds to the first (and only) dimension of filter_1D
        DimensionCorrespondence = zeros([1 numel(Size_MRStruct)]); DimensionCorrespondence(filter_dim) = 1;
		kSpaceFilter = myrepmat(filter_1D, Size_MRStruct,DimensionCorrespondence) .* kSpaceFilter;  % Replicate 1d-filter to the matrix size, 
                                                                                                    % Multiply this array with the previous calculated array
                                                                                                    % but now the filter-variation is in another dimension
                                                                                                    
	end

end



%% 3. Apply Filter

MRStruct.Data = MRStruct.Data .* kSpaceFilter;






%% 4. FFT to Image Space

if(~Settings.InputIskSpace_flag)
    
    for filter_dim = Settings.ApplyAlongDims
        MRStruct.Data = ifftshift(MRStruct.Data,filter_dim);
        MRStruct.Data = fft(MRStruct.Data,[],filter_dim);
        MRStruct.Data = fftshift(MRStruct.Data,filter_dim);
    end
    
end


%% 5. Postparations

