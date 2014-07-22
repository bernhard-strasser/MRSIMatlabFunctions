function IND = sub2ind_extended(ArraySize,AdditionalReplication,varargin)
%
% sub2ind_extended 
%
% This function was written by Bernhard Strasser, July 2014.
%
%
% The function masks the data in k-space, so that k-space values outside of an ellipsoid-shaped mask are set to zero. The mask can be a
% 3d-ellipsoid, or an 2d-ellipse. The equation for the mask is
% mask = {(x,y,z) E R³ | (x/a)² + (y/b)² + (z/c)² <= R²}
% a, b, c, and R can be chosen by the user.
%
%
% [OutArray,mask] = EllipticalFilter_x_y(OutArray,ApplyAlongDims,EllipsoidCoefficients,PerformFFT_flag)
%
% Input: 
% -     OutArray                     ...    Input array to which the filter should be applied
% -     ApplyAlongDims               ...    Along these dimensions the filter is applied. If this vector has two elements, a two dimensional 
%                                          Filter is applied. Otherwise, a 3d filter is used.
% -     EllipsoidCoefficients        ...    The values for [a b c R], which determine the shape and size of the ellipsoid. For two dimensional
%                                          Filter, set c = 1;
% -     PerformFFT_flag              ...    If it is 0, the image gets Fourier transformed to k-space before applying the filter, 
%                                          and transformed back to image domain afterwards
%
% Output:
% -     IND                     ...     c
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat_1_0

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


% 0.1 Preparations
if(nargin < 3)
	IND = 0;
	return
end
if(~exist('varargin','var') || numel(varargin) == 0)
	IND = 0;
	return
end

% 0.2 Declarations


% 0.3 Definitions
    



% Consistency Checks
if(isstruct(AdditionalReplication))
	AddRepNotWorking = ~isfield(AdditionalReplication,'DimPos') || ~isfield(AdditionalReplication,'Mat') || ~isfield(AdditionalReplication,'AddToDims') ...
	|| size(AdditionalReplication.Mat,2) ~= numel(AdditionalReplication.AddToDims);

	% If this is ok test if there would be non-positive indices
	if(~AddRepNotWorking)
		AddRepMin = min(AdditionalReplication.Mat,[],1);
		AddRepMax = max(AdditionalReplication.Mat,[],1);
		for dim = 1:numel(AdditionalReplication.AddToDims)
			if(min(varargin{AdditionalReplication.AddToDims(dim)}) + AddRepMin(dim) < 0 ||  ...
			max(varargin{AdditionalReplication.AddToDims(dim)}) + AddRepMax(dim) > ArraySize(AdditionalReplication.AddToDims(dim)))
				AddRepNotWorking = true;
				break;
			end
		end
	end
	if(AddRepNotWorking)
		fprintf('\n\nWARNING: AdditionalReplication has wrong properties, see help sub2ind_extended.\nIgnore AdditionalReplication variable.\n\n')
		AdditionalReplication = 0;
	end
end


% Make Input as wanted
% Make all input vectors to size N x 1
for cellno = 1:numel(varargin)
	varargin{cellno} = int16(reshape(varargin{cellno},[numel(varargin{cellno}) 1]));
end









%% 1.Compute linear indices


% What does the code do here? The computation of the linear indices looks quite complicated, and in fact -- it is :)
% Well in principle the target points are just all points all x- and y-values of all channels. The relative distance of the
% source points to a target point is computed within the kernel, and this relative information is applied to the target points
% in order to get the source points. Then the linear indices can be computed.

[sizes_x, sizes_y] = cellfun(@size,varargin); clear sizes_y;


if( isstruct(AdditionalReplication) )
	
end

sub2indString = 'sub2ind(ArraySize';
for dim = 1:numel(ArraySize)
	
	if( isstruct(AdditionalReplication) && sum(dim == AdditionalReplication.AddToDims) > 0 )
		AddDim = 0;
		if(AdditionalReplication.DimPos <= dim)
			AddDim = 1;
		end
		RepmatTo = ones([1 size(ArraySize,2)+numel(AdditionalReplication.DimPos)]); RepmatTo(AdditionalReplication.DimPos) = size(AdditionalReplication.Mat,1); RepmatTo(dim+AddDim) = numel(varargin{dim});
		DimCorr = zeros([1 AdditionalReplication.DimPos]); DimCorr2 = DimCorr; 
		DimCorr(AdditionalReplication.DimPos) = 1; DimCorr(dim+AddDim) = 0; DimCorr2(AdditionalReplication.DimPos) = 0; DimCorr2(dim+AddDim) = 1;
		AddPts = int16(myrepmat(AdditionalReplication.Mat(:,dim == AdditionalReplication.AddToDims),RepmatTo,DimCorr));
		varargin{dim} = myrepmat(varargin{dim},RepmatTo,DimCorr2) + AddPts;
	else
	
		% Replicate each of the vectors to an array of the outer products of all vectors
		DimCorr = zeros([1 numel(sizes_x)]); DimCorr(dim) = 1;
		varargin{dim} = myrepmat(varargin{dim},sizes_x,DimCorr);
		varargin{dim} = reshape(varargin{dim},[1 numel(varargin{dim})]);

	end
	
	sub2indString = [sub2indString ', uint32(varargin{' num2str(dim) '})'];
end
sub2indString = [sub2indString ');'];

IND = eval(sub2indString);



% Source_Channels = int16(transpose(1:nChannel));
% Source_Channels = repmat(Source_Channels, [1 no_SourcePoints nx_ACS_wo_border ny_ACS_wo_border]);
% 
% % Create spatial info
% Source_x = int16(kernelsize{KernelIndex}(1)+1:nx_ACS-kernelsize{KernelIndex}(2));
% Source_y = int16(kernelsize{KernelIndex}(3)+1:ny_ACS-kernelsize{KernelIndex}(4));
% 
% % Copy Source to Target Points
% Target_x = Source_x;
% Target_y = Source_y;
% 
% % Apply Relative Info
% Source_x = repmat(transpose(Source_x), [1 no_SourcePoints]);
% Source_y = repmat(transpose(Source_y), [1 no_SourcePoints]);
% Source_x = Source_x + repmat(reshape(int16(SrcRelativeTarg{KernelIndex}(:,1)),[1 no_SourcePoints]),[size(Source_x,1) 1]);
% Source_y = Source_y + repmat(reshape(int16(SrcRelativeTarg{KernelIndex}(:,2)),[1 no_SourcePoints]),[size(Source_y,1) 1]);
% 
% % Replicate spatial info
% Source_x = repmat(Source_x, [1 1 nChannel ny_ACS_wo_border]);
% Source_y = repmat(reshape(Source_y, [ny_ACS_wo_border no_SourcePoints]), [1 1 nChannel nx_ACS_wo_border]);
% 
% % Reorder Source Points
% Source_x = permute(Source_x, [3 2 1 4]);
% Source_y = permute(Source_y, [3 2 4 1]);
% 
% 
% % Target Points
% Target_Channels = reshape(Source_Channels(:,1,:,:), [nChannel nx_ACS_wo_border ny_ACS_wo_border]);
% Target_x = repmat(Target_x,[nChannel 1 ny_ACS_wo_border]);
% Target_y = repmat(reshape(Target_y,[1 1 numel(Target_y)]),[nChannel nx_ACS_wo_border 1]);
% 
% 
% % Linear Indices
% Target_linear = sub2ind( ...
% [nChannel nx_ACS ny_ACS], uint32(reshape(Target_Channels, [1 numel(Target_Channels)])), uint32(reshape(Target_x, [1 numel(Target_x)])), uint32(reshape(Target_y, [1 numel(Target_y)])));    
% Source_linear = sub2ind( ...
% [nChannel nx_ACS ny_ACS], uint32(reshape(Source_Channels, [1 numel(Source_Channels)])), uint32(reshape(Source_x, [1 numel(Source_x)])), uint32(reshape(Source_y, [1 numel(Source_y)])));
% % For Reconstructing MRSI data, uint32 has to be changed to uint64
% 
% if(max(Source_linear) > 2^31)
% 	max(Source_linear)
% 	display('Change to uint64 in code. Aborting.')
% 	OutData = InData; weights = 0;
% 	return
% end










%% 5. Postparations







