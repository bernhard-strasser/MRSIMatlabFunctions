function IND = sub2ind_extended(ArraySize,AdditionalReplication,EntangledDims,varargin)
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
%                                     If this is a little mystic 4 you, write edit sub2ind_extended and go to the
%                                     bottom, where there is an example from opencaipirinha.
% -      EntangledDims         ...    The array of the entangled dimensions.
%                                     If this option is used, e.g.
%                                     EntangledDims = [2 3], the inputs of varargin{2} and varargin{3}
%                                     will be considered 'entangled'. This means that those inputs will
%                                     be treated as encoding already the positions in the matrix of
%                                     ArraySize(2:3), e.g. if          
%                                     varargin{2} = [1 3 4] and varargin{3} = [6 7 7], then it will be
%                                     considered that the matrix elements (1,6), (3,7), (4,7) should be
%                                     encoded by varargin{2:3}. If the option was not used, the    
%                                     matrix elements
%                                     (1,6), (1,7), (4,7), (3,6), (3,7), (3,7), (4,6),(4,7),(4,7)
%                                     would be encoded instead.
% -     I1,I2,...,IN (varargin)       ...    The ranges for which the linear indices should be computed for all the dimensions
%                                     of ArraySize. Thus, numel(ArraySize) = N must hold.
%
% Output:
% -     IND                    ...    The linear indices of the array which is of size 
%                                     [1 numel(I1)*numel(I2)*...*numel(IN)*size(AdditionalReplication.Mat,1)]
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat

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
    



% % Consistency Checks
% if(isstruct(AdditionalReplication))
% 	AddRepNotWorking = ~isfield(AdditionalReplication,'DimPos') || ~isfield(AdditionalReplication,'Mat') || ~isfield(AdditionalReplication,'AddToDims') ...
% 	|| size(AdditionalReplication.Mat,2) ~= numel(AdditionalReplication.AddToDims);
% 
% 	% If this is ok test if there would be non-positive indices
% 	if(~AddRepNotWorking)
% 		AddRepMin = min(AdditionalReplication.Mat,[],1);
% 		AddRepMax = max(AdditionalReplication.Mat,[],1);
% 		for dim = 1:numel(AdditionalReplication.AddToDims)
% 			if(min(varargin{AdditionalReplication.AddToDims(dim)}) + AddRepMin(dim) < 0 ||  ...
% 			max(varargin{AdditionalReplication.AddToDims(dim)}) + AddRepMax(dim) > ArraySize(AdditionalReplication.AddToDims(dim)))
% 				AddRepNotWorking = true;
% 				break;
% 			end
% 		end
% 	end
% 	if(AddRepNotWorking)
% 		fprintf('\n\nWARNING: AdditionalReplication has wrong properties, see help sub2ind_extended.\nIgnore AdditionalReplication variable.\n\n')
% 		AdditionalReplication = 0;
% 	end
% end


% Make Input as wanted
% Make all input vectors to size N x 1
for cellno = 1:numel(varargin)
	varargin{cellno} = int16(reshape(varargin{cellno},[max(size(varargin{cellno})) min(size(varargin{cellno}))]));
end









%% 1.Compute linear indices


% What does the code do here? The computation of the linear indices looks quite complicated, and in fact -- it is :)
% Well in principle the target points are just all points all x- and y-values of all channels. The relative distance of the
% source points to a target point is computed within the kernel, and this relative information is applied to the target points
% in order to get the source points. Then the linear indices can be computed.

[AccessPointSize, sizes_y] = cellfun(@size,varargin);

% % Find out where varargin has a negative value
% NegValuePos = double(and(AccessPointSize == 1,varargin{AccessPointSize == 1} < 0)); NegValuePos(NegValuePos > 0) = -varargin{NegValuePos > 0};
% % Copy that varargin-element that is indicated by the absolute value of the negative, to the position of the negative values
% varargin{NegValuePos > 0} = varargin{NegValuePos(NegValuePos > 0)};
% % 
% AccessPartOfvararginMat = NegValuePos; AccessPartOfvararginMat(AccessPartOfvararginMat > 0) = 2:sum(AccessPartOfvararginMat > 0); AccessPartOfvararginMat(

if(numel(EntangledDims) > 1)
	AccessPointSize(EntangledDims(2:end)) = 1;
end


% We have to replicate also due to the AdditionalReplication
if( isstruct(AdditionalReplication) )
	AccessPointSize = [AccessPointSize(1:AdditionalReplication.DimPos-1) size(AdditionalReplication.Mat,1) AccessPointSize(AdditionalReplication.DimPos:end)];
end

AddDim = 0;
sub2indString = 'sub2ind(ArraySize';
for dim = 1:numel(ArraySize)
	
	DimCorr = zeros([1 numel(AccessPointSize)]); 
	if(numel(EntangledDims) > 1 && sum(dim == EntangledDims(2:end)) > 0)
		DimCorr(EntangledDims(1)+AddDim) = 1;	
	else
		DimCorr(dim+AddDim) = 1;
	end
	
	% Replicate the part of AdditionalReplication and add to varargin.
	if( isstruct(AdditionalReplication) && sum(dim == AdditionalReplication.AddToDims) > 0 )
		RepmatTo = [numel(varargin{dim}) size(AdditionalReplication.Mat,1)];
		AddPts = int16(myrepmat(AdditionalReplication.Mat(:,dim == AdditionalReplication.AddToDims),RepmatTo,[0 1]));
		varargin{dim} = myrepmat(varargin{dim},RepmatTo,[1 0]) + AddPts;
		if(AdditionalReplication.DimPos <= dim)
			varargin{dim} = transpose(varargin{dim});
			AddDim = 1;
			if(numel(EntangledDims) > 1 && sum(dim == EntangledDims(2:end)) > 0)
				DimCorr(EntangledDims(1)+AddDim) = 2;				
			else
				DimCorr(dim+AddDim) = 2;
			end
			DimCorr(AdditionalReplication.DimPos) = 1;
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






%% Example

% Since the code and its description is complicated and may be obscure for you, I added an example.
% The function call at the top results in the same as the bulk down (i.e. find(Source_linear2 ~= Source_linear) = empty-matrix, 
% i.e. find(Target_linear2 ~= Target_linear) = empty-matrix).


% AdditionalReplication.DimPos = 2; AdditionalReplication.AddToDims = [2 3]; AdditionalReplication.Mat = int16(SrcRelativeTarg{KernelIndex});
% Source_linear2 = sub2ind_extended([nChannel nx_ACS ny_ACS], ...
% 				 AdditionalReplication, ...
% 				 1:nChannel, kernelsize{KernelIndex}(1)+1:nx_ACS-kernelsize{KernelIndex}(2), kernelsize{KernelIndex}(3)+1:ny_ACS-kernelsize{KernelIndex}(4));
% Target_linear2 = sub2ind_extended([nChannel nx_ACS ny_ACS], ...
% 				 0, ...
% 				 1:nChannel, kernelsize{KernelIndex}(1)+1:nx_ACS-kernelsize{KernelIndex}(2), kernelsize{KernelIndex}(3)+1:ny_ACS-kernelsize{KernelIndex}(4));

% Compute the Source and Target Points as linear indices
% What does the code do here? Let's assume ... (ohh... did I fall asleep? Hm.)
% No, really: The source and the target points are computed for the ACS data. This is done by computing the linear indices of both
% (to avoid slow loops). The computation of the linear indices looks quite complicated, and in fact -- it is :)
% Well in principle the target points are just all points all x- and y-values of all channels. The relative distance of the
% source points to a target point is computed within the kernel, and this relative information is applied to the target points
% in order to get the source points. Then the linear indices can be computed.
%
% nx_ACS_wo_border = nx_ACS - sum(kernelsize{KernelIndex}(1:2));
% ny_ACS_wo_border = ny_ACS - sum(kernelsize{KernelIndex}(3:4));
% 
% % Source Points
% % Create channel info
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




%% 5. Postparations







