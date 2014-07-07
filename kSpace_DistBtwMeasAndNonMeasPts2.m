function QualityMeasure = kSpace_DistBtwMeasAndNonMeasPts2(PtsMeas, CellSize,EachNonMeasPtHasOwnKernel_flag)
%
% kSpace_Distance Compute quality measure of kSpace Pattern.
%
% This function was written by Bernhard Strasser, July 2013.
%
%
% This function 1) calculates the distance of a non-measured kSpace point to all the measured ones.
%               2) Does this for all non-measured points.
%               3) Calculates a quality measure for those distances (like maximum, mean or combination).
%               4) Does this for all non-measured points.
%               5) Computes the same quality measure using all non-measured points (maximum, mean or combination).
%
%
% [csi,csi_kspace] = read_csi_1_4(csi_path, zerofill_to_nextpow2_flag, zerofilling_fact, x_shift,y_shift)
%
% Input: 
% -         PtsMeas                     ...     The measured points as linear index (e.g. if A=[1 0;0 1], 1...measured --> PtsMeas = [1 4])
% -         CellSize                    ...     The size of the elementary caipirinha cell. in the above example: CellSize = [2 2]
%                                               Use at least 2*[Rx,Ry], where Rx and Ry would be the GRAPPA acceleration factors in x and y dir.
%                                               Otherwise you get only 1 possible undersampling cell.
%
% Output:
% -         QualityMeasure              ...     The quality measure of the kSpace Pattern, e.g. the mean distance of non-measured to measured points,
%                                               or the maximum of the distances between measured and non-measured ones.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: none yet.



%% 0. PREPARATIONS

% Assign standard values to variables if nothing is passed to function.
if(nargin < 2)
    display([ char(10) 'Gimme more input, Ma''am!' char(10) ])
    return;
end 
if(nargin < 3)
	EachNonMeasPtHasOwnKernel_flag = true;
end

MinKernelPts = 20;




%% 1. Create & Replicate Cell

ElCell = zeros(CellSize);
ElCell(PtsMeas) = 1;

KernelPtsPerElCell = sum(sum(ElCell));
MinReps = MinKernelPts / KernelPtsPerElCell;
MinReps(MinReps < 9) = 9;

PossibleReps = [1 3 5 7 9 11].^2;
PossibleReps(PossibleReps < MinReps) = NaN;
BestRep = nanmin(PossibleReps);
BestRep = sqrt(BestRep) + 2;


RepCell = repmat(ElCell,[BestRep BestRep]);





%% 2. Define Kernel Size(s)


CentralCell = [ floor(BestRep/2) floor(BestRep/2) ];

if(EachNonMeasPtHasOwnKernel_flag)
	NoKernels = sum(sum(~ElCell));
	kernelsize_std = [1 1 1 1];	
	
	CentralCellMatIndexFirstVoxel = CentralCell .* size(ElCell) + 1;
	CentralCellLinearIndex = zeros(size(RepCell));
	CentralCellLinearIndex(CentralCellMatIndexFirstVoxel(1) : CentralCellMatIndexFirstVoxel(1) + size(ElCell,1) - 1, ...
							   CentralCellMatIndexFirstVoxel(2) : CentralCellMatIndexFirstVoxel(2) + size(ElCell,2) - 1) = ~ElCell;
	CentralCellLinearIndex = find(CentralCellLinearIndex);
	[CentralCellMatIndex_x, CentralCellMatIndex_y] = ind2sub(size(RepCell), CentralCellLinearIndex);	
	
else
	NoKernels = 1;
	kernelsize_std = [0, size(ElCell,1)-1, 0, size(ElCell,2)-1];
	
	CentralCellMatIndex_x = CentralCell(1) .* size(ElCell,1) + 1;
	CentralCellMatIndex_y = CentralCell(2) .* size(ElCell,2) + 1;
	
end



CutOutCell = cell([1 NoKernels]);
kernelsize = cell([1 NoKernels]);
for KernelIndex = 1:NoKernels
	% Compute kernel size for processed Kernel
	kernelsize{KernelIndex} = kernelsize_std;
	
	no_SourcePoints = 0;
	while(no_SourcePoints < MinKernelPts)
		kernelsize{KernelIndex} = kernelsize{KernelIndex} + 1;
		kernel_dummy = RepCell(					CentralCellMatIndex_x(KernelIndex) - kernelsize{KernelIndex}(1): ...
												CentralCellMatIndex_x(KernelIndex) + kernelsize{KernelIndex}(2), ...
												CentralCellMatIndex_y(KernelIndex) - kernelsize{KernelIndex}(3): ...
												CentralCellMatIndex_y(KernelIndex) + kernelsize{KernelIndex}(4)		);

		no_SourcePoints = sum(sum(kernel_dummy));
	end

	% Test if the kernel size can be reduced in one direction without losing any Source points
	for kernelcrop = {[1 0 0 0], [0 1 0 0], [0 0 1 0], [0 0 0 1]}

		for crop_repeat = 1:max(size(ElCell))    % In GRAPPA-like kernels, one has to remove several lines/rows, if R_x > 2 | R_y > 2
			kernel_dummy_cropped = kernel_dummy(1 + kernelcrop{1}(1) : end - kernelcrop{1}(2), 1 + kernelcrop{1}(3) : end - kernelcrop{1}(4));

			no_SourcePoints_cropped = sum(sum(kernel_dummy_cropped));

			if(no_SourcePoints_cropped == no_SourcePoints)
				kernel_dummy = kernel_dummy_cropped;
				kernelsize{KernelIndex} = kernelsize{KernelIndex} - kernelcrop{1};
			end
		end

	end
	
	
	
	% Cut out the relevant matrix
	RepCell(CentralCellMatIndex_x(KernelIndex), CentralCellMatIndex_y(KernelIndex)) = 2;
	CutOutCell{KernelIndex} = RepCell(CentralCellMatIndex_x(KernelIndex) - kernelsize{KernelIndex}(1) : CentralCellMatIndex_x(KernelIndex) + kernelsize{KernelIndex}(2), ...
									  CentralCellMatIndex_y(KernelIndex) - kernelsize{KernelIndex}(3) : CentralCellMatIndex_y(KernelIndex) + kernelsize{KernelIndex}(4) );
	RepCell(CentralCellMatIndex_x(KernelIndex), CentralCellMatIndex_y(KernelIndex)) = 0;								  


end





%% 3. Compute the distance between all measured points for all non-measured ones


if(EachNonMeasPtHasOwnKernel_flag)
	MinDist = zeros([1 numel(CutOutCell)]);
	MinDist2 = zeros([1 numel(CutOutCell)]);
	for KernelIndex = 1:numel(CutOutCell)

		
		[PtsNotMeas_x, PtsNotMeas_y] = find(CutOutCell{KernelIndex} == 2);
		CutOutCell{KernelIndex}(CutOutCell{KernelIndex} == 2) = 0;
		PtsMeasKernel = find(CutOutCell{KernelIndex});

		

		% Find the x-y coordinates of the measured and not measured points within the cell.
		[PtsMeas_x, PtsMeas_y] = ind2sub(size(CutOutCell{KernelIndex}),PtsMeasKernel);
		

		% Compute the x and y distances
		Dist_x = abs(PtsMeas_x - PtsNotMeas_x);
		Dist_y = abs(PtsMeas_y - PtsNotMeas_y);	

		% Compute the Euclidean distance
		distance_mat = sqrt(Dist_x.^2 + Dist_y.^2);

		% Get min distance for each non-meas pt and weight that inversely with the number of occurences
		MinDist_dum = min(distance_mat,[],1);
		MinDist_dum2 = distance_mat; MinDist_dum2(distance_mat == repmat(MinDist_dum,[size(distance_mat,1) 1])) = Inf;
		MinDist_dum2 = min(MinDist_dum2,[],1); 

		MinDist_dum = MinDist_dum ./ sqrt(sum(distance_mat == repmat(MinDist_dum,[size(distance_mat,1) 1]),1));
		MinDist_dum2 = MinDist_dum2 ./ sqrt(sum(distance_mat == repmat(MinDist_dum2,[size(distance_mat,1) 1]),1));
		MinDist_dum2 = 0;
		
		MinDist(KernelIndex) = MinDist_dum;
		MinDist2(KernelIndex) = MinDist_dum2;
		

		%MinDist1 = min(distance_mat,[],1) ./ 

		% sum distance of closest 3 pts
		distance_sum3_dum = sort(distance_mat);
		distance_sum3(KernelIndex) = sum(distance_sum3_dum(1:3));

		% And the mean distance (mena over all measured points) 
		distance_mean = sum(distance_mat,1)/(numel(PtsMeas_x)^2);

		% And the min distance (min over all measured points) 
		distance_min = min(distance_mat,[],2);	






	end
	
else
	
	PtsMeasKernel = find(CutOutCell{1});
	PtsNotMeasKernel = find(~CutOutCell{1});
	
	% Replicate the Measured and Not Measured Points in such a way to get all combinations
	% (for every non-measured point we want the distance to all measured ones)
	PtsMeas_rep = repmat(PtsMeasKernel,[1 numel(PtsNotMeasKernel)]);
	PtsNotMeas_rep = repmat(transpose(PtsNotMeasKernel),[numel(PtsMeasKernel) 1]);	
	
	% Find the x-y coordinates of the measured and not measured points within the cell.
	[PtsMeas_x, PtsMeas_y] = ind2sub(size(CutOutCell{1}),PtsMeas_rep);
	[PtsNotMeas_x, PtsNotMeas_y] = ind2sub(size(CutOutCell{1}),PtsNotMeas_rep);
	
	% Compute the x and y distances
	Dist_x = abs(PtsMeas_x - PtsNotMeas_x);
	Dist_y = abs(PtsMeas_y - PtsNotMeas_y);	
	
	% Compute the Euclidean distance
	distance_mat = sqrt(Dist_x.^2 + Dist_y.^2);
	
	% Get min distance for each non-meas pt and weight that inversely with the number of occurences
	MinDist = min(distance_mat,[],1);
	MinDist2 = distance_mat; MinDist2(distance_mat == repmat(MinDist,[size(distance_mat,1) 1])) = Inf;
	MinDist2 = min(MinDist2,[],1); 
	
	MinDist = MinDist ./ sqrt(sum(distance_mat == repmat(MinDist,[size(distance_mat,1) 1]),1));
	MinDist2 = MinDist2 ./ sqrt(sum(distance_mat == repmat(MinDist2,[size(distance_mat,1) 1]),1));
	MinDist2 = 0;
	
	%MinDist1 = min(distance_mat,[],1) ./ 
	
	% sum distance of closest 3 pts
	distance_sum3 = sort(distance_mat);
	distance_sum3 = sum(distance_sum3(1:3));	
	
	% And the mean distance (mena over all measured points) 
	distance_mean = sum(distance_mat,1)/(numel(PtsMeas_x)^2);
	
	% And the min distance (min over all measured points) 
	distance_min = min(distance_mat,[],2);	
end











%% 1. Compute the quality measure

QualityMeasure(1) = max(distance_min);
QualityMeasure(2) = max(distance_mean);
QualityMeasure(3) = mean(MinDist+MinDist2);
QualityMeasure(4) = max(distance_sum3);







