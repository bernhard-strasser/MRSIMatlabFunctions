function [OutData,weights,kernelsize,SrcRelativeUpLeCorner,TrgRelativeUpLeCorner] = ...
opencaipirinha_MRSI_Fast(OutData, ACS, UndersamplingCell, quiet_flag, MinKernelSrcPts,weights,kernelsize,SrcRelativeUpLeCorner,TrgRelativeUpLeCorner) 
% 
% opencaipirinha_MRSI Reconstruct MRSI and MRI Data Undersampled With caipirinha Patterns
% 
%  [OutData,weights]=opencaipirinha_MRSI(InData, ACS, UndersamplingCell, MinKernelSrcPts) 
%       
%   Input:      
% -     InData                 Undersampled Data            (size: [#coils, nx, ny, nSlice, nTime]) (nTime = 1 for MRI). For memory reasons: InData = OutData;
%                              The undersampled points must be set to zero, NOT MISSING!
% -     ACS                    AutoCalibration Signal       (size: [#coils, nx_ACS, ny_ACS, nSlice_ACS])     
% -     UndersamplingCell      Elementary Cell, logical array (or array containing only 1's and 0's) which tells you the measured points (1's) and omitted points (0's). 
%                              Gets replicated to spatial size of InData to define the undersampling pattern. 
%                              Thus spatial size of InData and ACS must be integer multiple of UndersamplingCell.
% -     MinKernelSrcPts        The minimum source points in the kernel that are fitted to the target points. 20 is a good value, as that is standard in GRAPPA.
% -     weights                See Output. Used if the weights are already known and should not be computed from the ACS data.
% -     kernelsize             See Output. Used if the weights are already known and should not be computed from the ACS data. 
% -     SrcRelativeUpLeCorner  See Output. Used if the weights are already known and should not be computed from the ACS data.
% 
%   Output:
% -     OutData                Reconstructed Output Data    (size: [#coils, k_x, k_y, slc, nTime]) (nTime = 1 for MRI)
% -     weights                Weights for Reconstruction  
% -     kernelsize             The kernelsizes for each kPoint within the elementary cell  
% -     SrcRelativeUpLeCorner  For each target point (= non-measured kPoint) within the elementary cell this is a set of source points (= measured kPoints) around it
%                              that are used to reconstruct the target point.
%              
%




% November 2012 - November 2013 Bernhard Strasser
% September (?) 2013 Lucia Navarro (using precomputed GRAPPA weights, kernelsizes and SrcrelativeTarg to only reconstruct data, and not computing that from ACS data)




%   Some things to think about when using this code:
%
%           - The ACS lines used for reconstruction are NOT included in the final reconstructed
%             data sets. If you want to do this, you can do it after reconstruction. Please check
%             that they are lined up with the reconstructed data. I have not checked this.
%
%           - Since the ACS lines are not included, feel free to use a different imaging sequence
%             for acquisition of the ACS lines. We have seen some advantage in doing this when
%             the sequence used for the reduced acquisition is very flow sensitive, for example.
%             This can also be faster in many cases.
%
%           - The matrix problems here are normally so overdetermined that numerical conditioning 
%             is normally not necessary. You could try something out on your problem and see what
%             you get. I have never seen any real advantage.
%




% TODO: 
%       1)   If it is not too complicated, compute different kernels for the border voxels. 
%            Because the border voxels are computed with a lot of zeros, but with the "normal" reconstruction weights.
%            Yet, these "normal" weights assume that there is still data in all source points, and not zeros.
%       2)   A ball-shaped kernel instead of a cube might be a little better. Only little advantages expected.



%% 0. Preparation

% Assign standard values to variables if nothing is passed to function.

if(~exist('OutData','var'))
    display([char(10) 'You should consider inputting data which can be reconstructed. Aborting . . .'])
    weights = 0;    
    return
end

if(~exist('ACS','var'))
    display([char(10) 'I need an Auto Calibration Signal (ACS) for drinking GRAPPA/CAIPIRINHA! Aborting . . .'])
    weights = 0;    
    return
end

if(~exist('UndersamplingCell','var'))
    display([char(10) 'Please tell me how your data was undersampled by inputting ''UndersamplingCell''. Aborting . . .'])
    weights = 0;    
    return
end

if(~exist('MinKernelSrcPts','var'))
    MinKernelSrcPts = 20;
end

if(sum(sum(sum(UndersamplingCell))) == numel(UndersamplingCell))
    fprintf('\nNothing to do, as the UndersamplingCell only contains ones')
    weights = 0;
    return
end


if(~exist('quiet_flag','var'))
	quiet_flag = false;
end


% Further Preparations

% Get the size of both the input data and the autocalibration data
[nChannel,nx,ny, nSlice, nTime] = size(OutData);
[nChannel_ACS,nx_ACS,ny_ACS, nSlice_ACS]=size(ACS);

% Check for dimension size
if(nChannel_ACS~=nChannel)
    disp('Error! The number of coils has to be the same for both inputs! Aborting . . .')
    return;
end

% Convert UndersamplingCell to logical (Make us whole again)
UndersamplingCell = logical(UndersamplingCell);

% Define the Reduction Factor R
%R = numel(UndersamplingCell)/sum(sum(sum(UndersamplingCell)));

% Create Array with all Target (linear) indices, which are in the elementary cell
ElementaryCellLinearIndex = find(~UndersamplingCell);



%% 1. Calculate weights

if(~quiet_flag)
	tic;
end
if(~exist('weights','var'))

    % Fancy Text Message
	if(~quiet_flag)
		fprintf('\nFilling the glasses (Calculating weights)')
	end

    % What does the code do here? Let's assume an UndersamplingCell like this:
    % UndersmaplingCell = 1 0 1 0
    %                     0 1 1 0
    %                     1 0 1 1
    %                     0 1 0 1
    % 
    % with 1 = sampled point (source point), 0 = non-sampled point (target point)
    %
    % The zeros have to be reconstructed, thus one has to find weights which tells us how to combine the source points to get the target points.
    % This is achieved by fitting the source points to one target point at a time with apropriate kernels within the ACS data.
    % SO FOR EVERY TARGET POINT, A DIFFERENT KERNEL IS DEFINED, AND A DIFFERENT WEIGHTING SET IS COMPUTED!
    %
    % The devil in the detail: 
    % 1) Replicates UndersamplingCell to the size of the ACS data
    % 2) Loop over all target points within the central lower right cell
    % 3) Define for each of these target points the kernelsize, by steadily increasing the kernelsize 
    %    until the number of source points within it is higher than MinKernelSrcPts
    % 4) Try to reduce the kernelsize at all four sides (left, right, up, down with respect to the target point) without losing source points
    % 5) Test if the such defined kernel is equal to a kernel previously defined. If this is the case, just copy the weights for this kernel and skip the rest.
    % 6) For each kernel, loop over all possibilities to find that kernel within the ACS data 
    %    (the ACS data is fully sampled, so all points can be defined as source or target points, however it is necessary for the fitting).
    % 7) Gather all these possibilities in one big matrix.
    % 8) Compute the weights for the processed kernel
    % 9) Process next kernel.
    %
    % For visualization, run the script Caipirinha_Visualizing_SameAsFunction.m.




    % Create the undersampling pattern of the ACS data, to determine the source and target points for computing the weights
    UndersamplingPattern_ACS = repmat(UndersamplingCell, [floor(nx_ACS/size(UndersamplingCell,1)) floor(ny_ACS/size(UndersamplingCell,2))]);

    % Find the central lower right elementary cell replication within the UndersamplingPattern_ACS, and the indices of the first (upper left) voxel of this elementary cell.
    % In this cell, the kernelsize is computed for all kernels.
    CentralCell_ACS = [ floor(nx_ACS/(size(UndersamplingCell,1)*2)) floor(ny_ACS/(size(UndersamplingCell,2)*2)) ];
    CentralCellMatIndexFirstVoxel_ACS = CentralCell_ACS .* size(UndersamplingCell) + 1;
    %CentralCellLinearIndexFirstVoxel_ACS = sub2ind(size(UndersamplingPattern_ACS), CentralCellMatIndexFirstVoxel_ACS(1), CentralCellMatIndexFirstVoxel_ACS(2));

	if(~quiet_flag)
		fprintf('\nFilling Glass')
	end

	% Compute kernel size
	kernelsize = [0 size(UndersamplingCell,1)-1 0 size(UndersamplingCell,2)-1] - 1; % Always include the whole undersampling cell
	no_SourcePoints = 0;
	while(no_SourcePoints < MinKernelSrcPts)
		kernelsize = kernelsize + 1;

		kernel_dummy = UndersamplingPattern_ACS(CentralCellMatIndexFirstVoxel_ACS(1) - kernelsize(1): ...
												CentralCellMatIndexFirstVoxel_ACS(1) + kernelsize(2), ...
												CentralCellMatIndexFirstVoxel_ACS(2) - kernelsize(3): ...
												CentralCellMatIndexFirstVoxel_ACS(2) + kernelsize(4));

		no_SourcePoints = sum(sum(kernel_dummy));
	end


	% Test if the kernel size can be reduced in one direction without losing any Source points or making it smaller than the Undersampling cell
	kernel_dummy_sumrow = sum(kernel_dummy,1);
	kernel_dummy_sumcol = sum(kernel_dummy,2);
	KernelCropCellLog = ~[kernel_dummy_sumcol(1),kernel_dummy_sumcol(end),kernel_dummy_sumrow(1),kernel_dummy_sumrow(end)];
	KernelCropCellDummy = {[1 0 0 0],[0 1 0 0],[0 0 1 0],[0 0 0 1]};
	KernelCropCell = KernelCropCellDummy(KernelCropCellLog);

	for kernelcrop = KernelCropCell

		for crop_repeat = 1:max(size(UndersamplingCell))    % In GRAPPA-like kernels, one has to remove several lines/rows, if R_x > 2 | R_y > 2
			kernel_dummy_cropped = kernel_dummy(1 + kernelcrop{1}(1) : end - kernelcrop{1}(2), 1 + kernelcrop{1}(3) : end - kernelcrop{1}(4));

			no_SourcePoints_cropped = sum(sum(kernel_dummy_cropped));

			if(no_SourcePoints_cropped == no_SourcePoints && ~ismember(0,kernelsize(2:2:end) - kernelcrop{1}(2:2:end) >= size(UndersamplingCell)) )
				kernel_dummy = kernel_dummy_cropped;
				kernelsize = kernelsize - kernelcrop{1};
			else
				break;
			end
		end

	end


	% Compute the Target Points
	kernel = kernel_dummy; 
	clearvars kernel_dummy* no_SourcePoints_cropped kernelcrop ElementaryCellLinearIndex crop_repeat KernelCropCell*;
	kernel_TargetPoint = false(size(kernel));
	kernel_TargetPoint(kernelsize(1)+1:kernelsize(1)+size(UndersamplingCell,1),kernelsize(3)+1:kernelsize(3)+size(UndersamplingCell,2)) = ~UndersamplingCell;
	[TrgRelativeUpLeCorner_x,TrgRelativeUpLeCorner_y] = find(kernel_TargetPoint);
	TrgRelativeUpLeCorner = [TrgRelativeUpLeCorner_x-(kernelsize(1)+1),TrgRelativeUpLeCorner_y-(kernelsize(3)+1)];	% Here the reference point is simply defined as the left upper point of the
	no_TargetPoints = size(TrgRelativeUpLeCorner,1);
	
	
	% And the Source Points
	kernel_SourcePoints = kernel;


	[SrcRelativeUpLeCorner_x,SrcRelativeUpLeCorner_y] = find(kernel_SourcePoints);
	SrcRelativeUpLeCorner = [SrcRelativeUpLeCorner_x-(kernelsize(1)+1),SrcRelativeUpLeCorner_y-(kernelsize(3)+1)];	% Here the reference point is simply defined as the left upper point of the
	clear SrcRelativeUpLeCorner_x SrcRelativeUpLeCorner_y TargPos_x TargPos_y TrgRelativeUpLeCorner_x TrgRelativeUpLeCorner_y






	% Compute the Source and Target Points as linear indices
	% What does the code do here? Let's assume ... (ohh... did I fall asleep? Hm.)
	% No, really: The source and the target points are computed for the ACS data. This is done by computing the linear indices of both
	% (to avoid slow loops). The computation of the linear indices looks quite complicated, and in fact -- it is :)
	% Well in principle the target points are just all points all x- and y-values of all channels. The relative distance of the
	% source points to a target point is computed within the kernel, and this relative information is applied to the target points
	% in order to get the source points. Then the linear indices can be computed.


	AdditionalReplication.DimPos = 2; AdditionalReplication.AddToDims = [2 3]; AdditionalReplication.Mat = int16(SrcRelativeUpLeCorner);
	Source_linear = sub2ind_extended([nChannel nx_ACS ny_ACS], ...
					 AdditionalReplication,0, ...
					 1:nChannel, kernelsize(1)+1:nx_ACS-kernelsize(2), kernelsize(3)+1:ny_ACS-kernelsize(4));
	
	AdditionalReplication.Mat = int16(TrgRelativeUpLeCorner);
	Target_linear = sub2ind_extended([nChannel nx_ACS ny_ACS], ...
					 AdditionalReplication,0, ...
					 1:nChannel, kernelsize(1)+1:nx_ACS-kernelsize(2), kernelsize(3)+1:ny_ACS-kernelsize(4));

	if(max(Source_linear) > 2^31)
		max(Source_linear)
		display('Change to uint64 in code. Aborting.')
		weights = 0;
		return
	end

	weights = zeros([nChannel*no_TargetPoints nChannel*no_SourcePoints nSlice_ACS]);
	for SliceIndex = 1:nSlice_ACS

		% Slice ACS data
		ACS_sliced = ACS(:,:,:,SliceIndex);

		% The Target Points, size: nChannel x kernel repetitions in ACS data
		TargetPoints_ACS = ACS_sliced(Target_linear);
		TargetPoints_ACS = reshape(TargetPoints_ACS, [nChannel*no_TargetPoints numel(TargetPoints_ACS)/(nChannel*no_TargetPoints)]);

		% The Source Points, size: nChannel*no_SourcePoints x kernel repetitions in ACS data
		SourcePoints_ACS = ACS_sliced(Source_linear);    
		SourcePoints_ACS = reshape(SourcePoints_ACS, [nChannel*no_SourcePoints numel(SourcePoints_ACS)/(nChannel*no_SourcePoints)]);



		% find weights by fitting the source data to target data.
		% The pinv averages over all kernel repetitions in a weighted way ('least square' solution)
		% size: nChannel x nChannel*no_SourcePoints
		weights(:,:,SliceIndex) = TargetPoints_ACS * pinv(SourcePoints_ACS); 

	end




	if(~quiet_flag)
		fprintf('... %f sec \n',toc)                                                          
		fprintf('Aaaaahhh! This SMELL! \n') 
	end

else
	if(~quiet_flag)
		fprintf('... %f sec \n',toc)                                                          
		fprintf('Aaaaahhh! This SMELLS much BETTER! \n') 
	end
end

clearvars *ACS*;





%% 2. Apply weights

% Fancy Text Message
if(~quiet_flag)
	tic; fprintf('\nDrinking CAIPIRINHA (Applyig weights) ')
end


% Create an extended matrix, extended by the kernelsize. This extension is necessary, so that also the target points at the border of the matrix can be reconstructed, using the kernel.
% Write undersampled data, and zero pad this data


% Prepare the enlarged/extended Reconstruction Matrix
Reco_dummy = zeros([nChannel, nx+sum(kernelsize(1:2))-(size(UndersamplingCell,1)-1), ny+sum(kernelsize(3:4))-(size(UndersamplingCell,2)-1), nSlice, nTime]); % This has to be changed for "TODO 2)"
Reco_dummy(:,kernelsize(1)+1:kernelsize(1)+size(OutData,2),kernelsize(3)+1:kernelsize(3)+size(OutData,3),:,:) = OutData;
Reco_dummy_CheckIfRecoNecess = squeeze(abs(Reco_dummy(1,:,:,1,1)) > 0);

if(~quiet_flag)
	fprintf('\nDrinking Glass\tGUUULLPP . . . UAHHHH . . . ')
end

% Compute all the linear indices of the target points within Reco_dummy of the processed Kernel
% UndersamplingCell_CurrentTargetPoint = false(size(UndersamplingCell));											% Mark the current target point in the undersampling cell
% UndersamplingCell_CurrentTargetPoint(ElementaryCellLinearIndex) = true; 
% TargetPoints = false([size(Reco_dummy,2) size(Reco_dummy,3)]);


% This has to be changed for "TODO 2)"
% What is done here?
% 1) The inner part of the TargetPoints are called, because the target points can only be in that part of Reco_dummy which implements the InData, not the zerofilled part 
%    (we dont want to reconstruct the zeros, do we?)
% 2) The UndersamplinCell with the current target point gets replicated to the size of the InData. If the InData is not an integer multiple of the UndersamplingCell then
%    this will lead to an error here.
% 3) Then the x- and y-coordinates of the target points are computed.
% UC_dummy = repmat(UndersamplingCell_CurrentTargetPoint, [floor(nx/size(UndersamplingCell,1)) floor(ny/size(UndersamplingCell,2))]);
% UC_dummy = cat(  1,UC_dummy, UC_dummy(1:nx-size(UC_dummy,1),:,:)  );
% UC_dummy = cat(  2,UC_dummy, UC_dummy(:,1:ny-size(UC_dummy,2),:)  );
% %UC_dummy = cat(  3,UC_dummy, UC_dummy(:,:,1:nSlice-size(UC_dummy,3))  );  
% TargetPoints(kernelsize(1)+1:end-kernelsize(2),kernelsize(3)+1:end-kernelsize(4),:) = UC_dummy;
% [TargetPoints_x, TargetPoints_y] = find(TargetPoints);


% Loop over all target points for the processed kernel
SourcePoints = zeros([nChannel*size(SrcRelativeUpLeCorner,1) nTime]);
for xLoopy = kernelsize(1)+1:size(UndersamplingCell,1):size(Reco_dummy,2)-kernelsize(2)
	for yLoopy = kernelsize(3)+1:size(UndersamplingCell,2):size(Reco_dummy,3)-kernelsize(4)
	
		SourceIndices = [xLoopy + SrcRelativeUpLeCorner(:,1), yLoopy + SrcRelativeUpLeCorner(:,2)];
		
		if(sum(Reco_dummy_CheckIfRecoNecess(sub2ind(size(Reco_dummy_CheckIfRecoNecess),SourceIndices(:,1),SourceIndices(:,2)))) < 2)
			continue
		end
		
		
		
		TargetIndices = [xLoopy + TrgRelativeUpLeCorner(:,1), yLoopy + TrgRelativeUpLeCorner(:,2)];
		
		for SliceIndex = 1:nSlice
			% Compute the Source Points for that specific TargetIndex
			for Sourceloopy = 1:size(SourceIndices,1)
				SourcePoints((Sourceloopy-1)*nChannel+1:Sourceloopy*nChannel,:) = Reco_dummy(:,SourceIndices(Sourceloopy,1),SourceIndices(Sourceloopy,2),SliceIndex,:);
			end		
			% Reconstruct data by applying weights to Sourcepoints.			
			Reco_dummy2 = reshape(weights(:,:,SliceIndex)*SourcePoints, [nChannel size(TargetIndices,1) 1 1 nTime]);
			for Targetloopy = 1:size(TargetIndices,1)
				Reco_dummy(:,TargetIndices(Targetloopy,1),TargetIndices(Targetloopy,2),SliceIndex,:)=Reco_dummy2(:,Targetloopy,:);
			end
		end
		
	end
end




OutData = Reco_dummy(:,kernelsize(1)+1:end-kernelsize(2),kernelsize(3)+1:end-kernelsize(4),:,:);                                           %Crop out the good data.

if(~quiet_flag)
	fprintf('... %f sec \n',toc)
end






%% ~. Postparations

% Fancy Text Message
if(~quiet_flag)
	fprintf('Mmmmhhh, I love CAIPIRINHA! \n')    
end


