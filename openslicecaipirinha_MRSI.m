function [OutData,weights]=openslicecaipirinha_MRSI(InData, ACS_or_weights, FoV_shifts, quiet_flag,precision, AliasedSlices, kernelsize) 
% 
% openslicegrappa_MRSI Reconstruct the Slices of MRSI and MRI Data Ehen Only the Sum of Those Slices Was Measured.
% 
%  [OutData,weights] = openslicegrappa_MRSI(InData, ACS_or_weights, FoV_shifts, AliasedSlices, kernelsize);
%       
%   Input:      
% -     InData             Undersampled Data            (size: [#coils, nx, ny, nSlice, nTime]) (nTime = 1 for MRI)
%                          nSlice must be nSlice = nSlice_ACS/R_Slice, where R_Slice is the acceleration factor in slice direction.
% -     ACS_or_weights     AutoCalibration Signal       (size: [#coils, nx_ACS, ny_ACS, nSlice_ACS])  
%                          Must have as many slices as the OutData should have. Use similar sequence parameters as for the InData.
%                          If ACS_or_weights is a cell, ACS_or_weights are treated as weights, and the weights are not computed then.
% -     FoV_shifts         The multiples of the FoV with which each slice was shifted 
%                          (against the original slice position) in x- and y direction.
%                          For each slice and for both directions (x and y) this must be given.
%                          E.g.: [0 0; 0.2 0.3; 0.5 0.5] would shift slice 0 by nothing,
%                          slice 1 by 0.2*FoV_x in x and 0.3*FoV_y in x direction and
%                          slice 3 by 0.5*FoV in both directions.
%                          size: [nSlice_ACS 2]
% -     kernelsize         The kernelsize. Input should be either one of those formats (examples produce exactly the same results):
%                               *   [1 4]-vector: These give the maximum distance of the kernel from the target point
%                                   to the left and right (x direction), and to down and up (y direction). E.g.: [2 2 2 2]
%                               *   [1 2]-vector: The kernelsize in x- and y-direction. E.g. [5 5]
%                               *   scalar: The kernel size in both directions. E.g. 5.
% -     AliasedSlices      A matrix telling the function which slices are aliased.
%                          E.g.: [1 3 5; 2 4 6]. This tells the program:
%                          AliasedSlices(1,:) = 1 3 5 --> Slice 1 of InData contains the aliased Slices of 1,3 and 5
%                          AliasedSlices(2,:) = 2 4 6 --> Slice 2 of Indata contains the aliased Slices of 2,4 and 6
% 
%   Output:
% -     OutData            Reconstructed Output Data    (size: [#coils, nx, ny, nSlice_ACS, nTime]) (nTime = 1 for MRI)
% -     weights            Weights for Reconstruction              
%              
%

% August 2013 - Bernhard Strasser


% TODO: 
%       1)   If it is not too complicated, compute different kernels for the border voxels. 
%            Because the border voxels are computed with a lot of zeros, but with the "normal" reconstruction weights.
%            Yet, these "normal" weights assume that there is still data in all source points, and not zeros.
%       2)   A ball-shaped kernel instead of a cube might be a little better. Only little advantages expected.





%% 0. Preparation




% Assign standard values to variables if nothing is passed to function.

if(~exist('InData','var'))
    display([char(10) 'I need Cachaca ( = data which should be processed (InData)) for preparing CAIPIRINHA! Aborting . . .'])
    return
end

if(~exist('ACS_or_weights','var'))
    display([char(10) 'I need Glasses ( = Auto Calibration Signal (ACS)) for preparing CAIPIRINHA! Aborting . . .'])
    return
end
if(~exist('quiet_flag','var'))
	quiet_flag = false;
end
if(~exist('precision','var') || ~ischar(precision))
	precision = 'double';
end


% Initialize kernelsize
if(~exist('kernelsize','var'))
    kernelsize = [2 2 2 2];
end
if(numel(kernelsize)==1)
    kernelsize = repmat(kernelsize(1),[1 4]);
elseif(numel(kernelsize) < 4)
    kernelsize = [floor(kernelsize(1)/2) floor(kernelsize(1)/2) floor(kernelsize(2)/2) floor(kernelsize(2)/2)];
end


% Fancy Text Message
if(~quiet_flag)
	fprintf('\n\nLet the lemon-slice-caipirinha party start!')
end


% Further Preparations


% kernelsize
kernelsize_x = sum(kernelsize(1:2))+1;
kernelsize_y = sum(kernelsize(3:4))+1;

% Get the size of both the input data and the autocalibration data
[nChannel,nx,ny, nSlice, nTime] = size(InData);

if(iscell(ACS_or_weights))
    weights = ACS_or_weights; clear ACS_or_weights
    
    if(exist('FoV_shifts','var'))
        nSlice_ACS = size(FoV_shifts,1);
    elseif(exist('AliasedSlices','var'))
        nSlice_ACS = numel(AliasedSlices);
    else
        display([char(10) 'I need either FoV-Shifts or AliasedSlices! Aborting . . .'])
        return
    end    
    
else
    
    ACS = ACS_or_weights; clear ACS_or_weights
    [nChannel_ACS,nx_ACS,ny_ACS, nSlice_ACS]=size(ACS);

    % Check for dimension size
    if(nChannel_ACS~=nChannel)
        disp([char(10) 'Error! The number of coils has to be the same for both inputs! Aborting . . .'])
        OutData = InData;
        weights = 0;
        return;
    end

end


if(nSlice == nSlice_ACS)
    disp([char(10) 'Nothing to do, since the InData has the same amount of slices as the ACS data . . .'])
    OutData = InData;
    weights = 0;
    return;
end

% Initialize FoV-shifts
if(~exist('FoV_shifts','var'))
    FoV_shifts = zeros([nSlice_ACS 2]);
end

% Initialize AliasedSlices
if(~exist('AliasedSlices','var'))                           % Produce patterns like [1 3; 2 4; 1 3; 2 4] or [1 4; 2 5; 3 6; 1 4; 2 5; 3 6].
    R_Slice = nSlice_ACS/size(InData,4);
    nSlice = size(InData,4);    
    AliasedSlices = zeros([nSlice R_Slice]);
    for dummy = 1:nSlice
       AliasedSlices(dummy,:) = dummy:nSlice:nSlice_ACS;
    end 
    %AliasedSlices = repmat(AliasedSlices,[R_Slice 1]);
    clear R_Slice dummy
end



%% 1. Shift the ACS FoV

if(~exist('weights','var'))
    ACS = kSpace_FoVShift(ACS,FoV_shifts);
end




%% 2. Calculate weights


if(~exist('weights','var'))
    
    % Fancy Text Message
	if(~quiet_flag)
		tic;
		fprintf('\nFilling the glasses and cutting lemon slices (Calculating weights)')
	end

    % What does the code do here? Let's assume ... (ohh... did I fall asleep? Hm.)
    % No, really: The source and the target points are computed for the ACS data. This is done by computing the linear indices of both
    % (to avoid slow loops). The computation of the linear indices looks quite complicated, and in fact -- it is :)
    % Well in principle the target points are just all points all x- and y-values of all channels. The relative distance of the
    % source points to a target point is computed within the kernel, and this relative information is applied to the target points
    % in order to get the source points. Then the linear indices can be computed.
    % For visualization, run the script SliceGrappa_Visualizing.m


    % Compute the Source and Target Points as linear indices
    nx_ACS_wo_border = nx_ACS - sum(kernelsize(1:2));
    ny_ACS_wo_border = ny_ACS - sum(kernelsize(3:4));

    kernel_Src = ones([sum(kernelsize(1:2))+1, sum(kernelsize(3:4))+1]);
    kernel_Targ = zeros([sum(kernelsize(1:2))+1, sum(kernelsize(3:4))+1]);
    kernel_Targ(kernelsize(1)+1,kernelsize(2)+1) = 1;

    no_SrcPts = sum(sum(sum(kernel_Src)));

    [SrcRelativeTarg_x,SrcRelativeTarg_y] = find(kernel_Src);
    [TargPos_x, TargPos_y] = find(kernel_Targ);
    SrcRelativeTarg = [SrcRelativeTarg_x-TargPos_x,SrcRelativeTarg_y-TargPos_y];
    clear SrcRelativeTarg_x SrcRelativeTarg_y TargPos_x TargPos_y

    % Source Points
    % Create channel info
    Source_Channels = int16(transpose(1:nChannel));
    Source_Channels = repmat(Source_Channels, [1 no_SrcPts nx_ACS_wo_border ny_ACS_wo_border]);

    % Create spatial info
    Source_x = int16(kernelsize(1)+1:nx_ACS-kernelsize(2));
    Source_y = int16(kernelsize(3)+1:ny_ACS-kernelsize(4));

    % Copy Source to Target Points
    Target_x = Source_x;
    Target_y = Source_y;

    % Apply Relative Info
    Source_x = repmat(transpose(Source_x), [1 no_SrcPts]);
    Source_y = repmat(transpose(Source_y), [1 no_SrcPts]);
    Source_x = Source_x + repmat(reshape(int16(SrcRelativeTarg(:,1)),[1 no_SrcPts]),[size(Source_x,1) 1]);
    Source_y = Source_y + repmat(reshape(int16(SrcRelativeTarg(:,2)),[1 no_SrcPts]),[size(Source_y,1) 1]);

    % Replicate spatial info
    Source_x = repmat(Source_x, [1 1 nChannel ny_ACS_wo_border]);
    Source_y = repmat(reshape(Source_y, [ny_ACS_wo_border no_SrcPts]), [1 1 nChannel nx_ACS_wo_border]);

    % Reorder Source Points
    Source_x = permute(Source_x, [3 2 1 4]);
    Source_y = permute(Source_y, [3 2 4 1]);


    % Target Points
    Target_Channels = reshape(Source_Channels(:,1,:,:), [nChannel nx_ACS_wo_border ny_ACS_wo_border]);
    Target_x = repmat(Target_x,[nChannel 1 ny_ACS_wo_border]);
    Target_y = repmat(reshape(Target_y,[1 1 numel(Target_y)]),[nChannel nx_ACS_wo_border 1]);


    % Linear Indices
    Target_linear = sub2ind( ...
    [nChannel nx_ACS ny_ACS], uint32(reshape(Target_Channels, [1 numel(Target_Channels)])), uint32(reshape(Target_x, [1 numel(Target_x)])), uint32(reshape(Target_y, [1 numel(Target_y)])));    
    Source_linear = sub2ind( ...
    [nChannel nx_ACS ny_ACS], uint32(reshape(Source_Channels, [1 numel(Source_Channels)])), uint32(reshape(Source_x, [1 numel(Source_x)])), uint32(reshape(Source_y, [1 numel(Source_y)])));
    % For Reconstructing MRSI data, uint32 has to be changed to uint64

    if(max(Source_linear) > 2^31)
        max(Source_linear)
        display('Change to uint64 in code. Aborting.')
        OutData = InData; weights = 0;
        return
    end




    % Iterate over all Slices
    weights = cell([1 nSlice_ACS]);
    for SourceSliceIndex = 1:nSlice                                 % Iterate over the slices of the InData (so the aliased or source slices)

        % The Source Points, size: nChannel*no_SrcPts x kernel repetitions in ACS data
        ACS_sum = sum(ACS(:,:,:,AliasedSlices(SourceSliceIndex,:)),4);
        SourcePoints_ACS = ACS_sum(Source_linear);
        SourcePoints_ACS = reshape(SourcePoints_ACS, [nChannel*no_SrcPts numel(SourcePoints_ACS)/(nChannel*no_SrcPts)]);     

        for TargetSliceIndex = AliasedSlices(SourceSliceIndex,:)                    % Iterate over the corresponding (corr to source slices) Target Slices.
                                                                                    % So if Slice 1 of the InData contains the aliased slices [1 3 5], then compute the weights
            % The Target Points, size: nChannel x kernel repetitions in ACS data    % For Slices 1,3 and 5 to be able to reconstruct those out of Slice 1 of the InData.
            ACS_sliced = ACS(:,:,:,TargetSliceIndex);    
            TargetPoints_ACS = ACS_sliced(Target_linear);
            TargetPoints_ACS = reshape(TargetPoints_ACS, [nChannel numel(TargetPoints_ACS)/nChannel]);        

            % find weights by fitting the source data to target data.
            % The pinv averages over all kernel repetitions in a weighted way ('least square' solution)
            % size: nChannel x nChannel*no_SrcPts
            weights{TargetSliceIndex} = TargetPoints_ACS * pinv(SourcePoints_ACS); 

        end
	end

	if(~quiet_flag)
		fprintf('\nCutting took\t... %f sec \n',toc)                                                          
		fprintf('Aaaaahhh! This SMELL! \n')
	end

end





%% 3. Apply weights

% Fancy Text Message
if(~quiet_flag)
	tic; fprintf('\nDrinking Caipirinha with slices of lemon (Applyig weights).')
end


% Create an extended matrix, extended by the kernelsize. This extension is necessary, so that also the target points at the border of the matrix can be reconstructed, using the kernel.
Reco_dummy = feval(precision,zeros([nChannel nx+sum(kernelsize(1:2)) ny+sum(kernelsize(3:4)) nSlice nTime]));	% Convert to singles if necessary
Reco_dummy(:,kernelsize(1)+1:end-kernelsize(2),kernelsize(3)+1:end-kernelsize(4),:,:) = InData;
InPlane_UndersamplingPattern = abs(squeeze(Reco_dummy(1,:,:,1,1))) > 0;

% Create the OutData which has to have as many slices as the ACS data because we want to reconstruct all those slices.
OutData = feval(precision,zeros([nChannel nx ny nSlice_ACS nTime]));
for SourceSliceIndex = 1:nSlice
	
	if(~quiet_flag)
		fprintf('\nDrinking Caipirinha with lemon slice %d\tGUUULLPP . . . UAHHHH . . . ',SourceSliceIndex)   
	end
	weights_temp = cat(1,weights{AliasedSlices(SourceSliceIndex,:)});

    % Loop over all target points for the processed kernel
	for yloopy = kernelsize(1)+1:kernelsize(1)+size(OutData,2)
		for xloopy = kernelsize(3)+1:kernelsize(3)+size(OutData,3)
			if(InPlane_UndersamplingPattern(xloopy,yloopy) == 0)
				continue;
			end
			% Sourcepoints are the points of all channels and all x and y points within the kernel around the target point.
			SourcePoints = reshape(Reco_dummy(:,xloopy-kernelsize(1):xloopy+kernelsize(2), ...
			yloopy-kernelsize(3):yloopy+kernelsize(4),SourceSliceIndex,:), [nChannel*kernelsize_x*kernelsize_y nTime]);        	

			OutData_dummy = reshape(weights_temp*SourcePoints, [nChannel 1 1 size(AliasedSlices,2) nTime]);
			OutData(:,xloopy-kernelsize(1),yloopy-kernelsize(3),AliasedSlices(SourceSliceIndex,:),:) = OutData_dummy;	
		end
	end

end


if(~quiet_flag)
	fprintf('\nDrinking took\t...\t%f sec.',toc)
end
clearvars -except quiet_flag OutData FoV_shifts weights



%% 4. Undo FoV-Shifts

OutData = kSpace_FoVShift(OutData,-FoV_shifts);




%% ~. Postparations

% Fancy Text Message
if(~quiet_flag)
	fprintf('\nMmmmhhh, I love Caipirinha! \n')    
end


