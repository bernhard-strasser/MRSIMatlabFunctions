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

[AccessPointSize, sizes_y] = cellfun(@size,varargin); clear sizes_y;


% We have to replicate also
if( isstruct(AdditionalReplication) )
	AccessPointSize = [AccessPointSize(1:AdditionalReplication.DimPos-1) size(AdditionalReplication.Mat,1) AccessPointSize(AdditionalReplication.DimPos:end)];
end

sub2indString = 'sub2ind(ArraySize';
for dim = 1:numel(ArraySize)
	DimCorr = zeros([1 numel(AccessPointSize)]); DimCorr(dim) = 1;
	
	
	% Replicate the part of AdditionalReplication and add to varargin.
	if( isstruct(AdditionalReplication) && sum(dim == AdditionalReplication.AddToDims) > 0 )
		RepmatTo = [numel(varargin{dim}) size(AdditionalReplication.Mat,1)];
		AddPts = int16(myrepmat(AdditionalReplication.Mat(:,dim == AdditionalReplication.AddToDims),RepmatTo,[0 1]));
		varargin{dim} = myrepmat(varargin{dim},RepmatTo,[1 0]) + AddPts;
		if(AdditionalReplication.DimPos <= dim)
			varargin{dim} = transpose(varargin{dim});
			DimCorr(dim) = 0; DimCorr(dim+1) = 2; DimCorr(AdditionalReplication.DimPos) = 1;
		else
			DimCorr(AdditionalReplication.DimPos) = 2;
		end
	end
	
	% Replicate each of the vectors to an array of the outer products of all vectors

	varargin{dim} = myrepmat(varargin{dim},AccessPointSize,DimCorr);
	varargin{dim} = reshape(varargin{dim},[1 numel(varargin{dim})]);

	
	sub2indString = [sub2indString ', uint32(varargin{' num2str(dim) '})'];
end
sub2indString = [sub2indString ');'];

IND = eval(sub2indString);








%% 5. Postparations







