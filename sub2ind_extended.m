function IND = sub2ind_extended(ArraySize,AdditionalReplication,varargin)
%
% sub2ind_extended Compute the linear indices of an array with extended features w.r.t. sub2ind.
%
% This function was written by Bernhard Strasser, July 2014.
%
%
% The function does the same as the internal sub2ind function, but with 2 additional features:
% 1) If you specify a range for I1,I2 etc. with min(In) <= size(ArraySize,n) <= max(In), you get a whole list of
%    indices. In the original sub2ind you have to can do this only indirectly, by giving the whole list already,
%    e.g. IND = sub2ind([2 2 2],[1 1 1 1],[1 1 2 2],[1 2 1 2]);
%    in comparison, this can achieved by this function simply by
%    IND = sub2ind_extended([2 2 2],0,1,1:2,1:2);
%    which is less flexible (the order of the output indices is determined) but much easier to use.
%
% 2) You can also compute the linear indices when additionally replicating the SIZ along one dimension
%    and add a relative position matrix along this dimension.
%    Example: If you have an 3x4x5-array, and for each of these points you want to have a set of 6 points
%    which can be computed by taking the y- and z-point of the array and adding
%    [-1 -1 -1 0 0 0] to the y-points and [-1 0 1 -1 0 1] to the z-points, then you can do this by
%    setting AdditionalReplication.Mat = [-1 -1; -1 0; -1 1; 0 -1; 0 0; 0 1] and 
%    AdditionalReplication.AddToDims = [2 3]. AdditionalReplication.DimPos can be set to any value,
%    the output will differ.
%
%
% IND = sub2ind_extended(SIZ,AdditionalReplication,I1,I2,...,IN)
%
% Input: 
% -     ArraySize              ...    Size of Array for which the linear indices should be computed.
% -     AdditionalReplication  ...    Structure with fields (all mandatory):
%                                     * .Mat
%                                     * .AddToDims
%                                     * .DimPos
%                                     If there should be performed an additional replication along dimension n,
%                                     so if you wantfor each point a set of additional points which have a
%                                     relative position of Mat to the other original points of dimensions AddToDims,
%                                     you can achieve this by setting AdditionalReplication.DimPos = n; 
%                                     AdditionalReplication.AddToDims = AddToDims; AdditionalReplication.Mat = Mat;
% -     I1,I2,...,IN           ...    The ranges for which the linear indices should be computed for all the dimensions
%                                     of ArraySize. Thus, numel(ArraySize) = N must hold.
%
% Output:
% -     IND                    ...    The linear indices of the array which is of size 
%                                     [1 numel(I1)*numel(I2)*...*numel(IN)*size(AdditionalReplication.Mat,1)
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


% We have to replicate also due to the AdditionalReplication
IncreasedDim = 0;
if( isstruct(AdditionalReplication) )
	AccessPointSize = [AccessPointSize(1:AdditionalReplication.DimPos-1) size(AdditionalReplication.Mat,1) AccessPointSize(AdditionalReplication.DimPos:end)];
	ArraySize = [ArraySize(1:AdditionalReplication.DimPos-1) size(AdditionalReplication.Mat,1) ArraySize(AdditionalReplication.DimPos:end)];
	IncreasedDim = 1;
end

sub2indString = 'sub2ind(ArraySize';
for dim = 1:numel(ArraySize)-IncreasedDim
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







