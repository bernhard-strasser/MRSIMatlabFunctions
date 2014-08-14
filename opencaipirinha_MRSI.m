function [OutData,weights,kernelsize,SrcRelativeTarg]=opencaipirinha_MRSI(OutData, ACS, UndersamplingCell, quiet_flag, MinKernelSrcPts,weights,kernelsize,SrcRelativeTarg) 
% 
% opencaipirinha_MRSI Reconstruct MRSI and MRI Data Undersampled With caipirinha Patterns
% 
%  [OutData,weights]=opencaipirinha_MRSI(InData, ACS, UndersamplingCell, MinKernelSrcPts) 
%       
%   Input:      
% -     InData             Undersampled Data            (size: [#coils, nx, ny, nSlice, nTime]) (nTime = 1 for MRI). For memory reasons: InData = OutData;
%                          The undersampled points must be set to zero, NOT MISSING!
% -     ACS                AutoCalibration Signal       (size: [#coils, nx_ACS, ny_ACS, nSlice_ACS])     
% -     UndersamplingCell  Elementary Cell, logical array (or array containing only 1's and 0's) which tells you the measured points (1's) and omitted points (0's). 
%                          Gets replicated to spatial size of InData to define the undersampling pattern. 
%                          Thus spatial size of InData and ACS must be integer multiple of UndersamplingCell.
% -     MinKernelSrcPts    The minimum source points in the kernel that are fitted to the target points. 20 is a good value, as that is standard in GRAPPA.
% -     weights            See Output. Used if the weights are already known and should not be computed from the ACS data.
% -     kernelsize         See Output. Used if the weights are already known and should not be computed from the ACS data. 
% -     SrcRelativeTarg    See Output. Used if the weights are already known and should not be computed from the ACS data.
% 
%   Output:
% -     OutData            Reconstructed Output Data    (size: [#coils, k_x, k_y, slc, nTime]) (nTime = 1 for MRI)
% -     weights            Weights for Reconstruction  
% -     kernelsize         The kernelsizes for each kPoint within the elementary cell  
% -     SrcRelativeTarg    For each target point (= non-measured kPoint) within the elementary cell this is a set of source points (= measured kPoints) around it
%                          that are used to reconstruct the target point.
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

% Number of different kernels. For each kernel, the reconstruction weights are computed.
nKernels = sum(sum(sum(~UndersamplingCell)));


% Create Array with all Target (linear) indices, which are in the elementary cell
ElementaryCellLinearIndex = find(~UndersamplingCell);

% Create a Vector with the Kernel Correspondance (which kernels are the same)
KernelCorrespondence = zeros([1 nKernels]);


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


    % Find indices of target points, which are in the above computed central cell
    % Around every of these points, a kernel will be defined, and for all these kernels, the weights will be computed.
    CentralCellLinearIndex_ACS = zeros(size(UndersamplingPattern_ACS));
    CentralCellLinearIndex_ACS(CentralCellMatIndexFirstVoxel_ACS(1) : CentralCellMatIndexFirstVoxel_ACS(1) + size(UndersamplingCell,1) - 1, ...
                               CentralCellMatIndexFirstVoxel_ACS(2) : CentralCellMatIndexFirstVoxel_ACS(2) + size(UndersamplingCell,2) - 1) = ~UndersamplingCell;
    CentralCellLinearIndex_ACS = find(CentralCellLinearIndex_ACS);
    [CentralCellMatIndex_ACS_x, CentralCellMatIndex_ACS_y] = ind2sub(size(UndersamplingPattern_ACS), CentralCellLinearIndex_ACS);



    % Preallocate kernelsize, kernel and the kernel_TargetPoint (indicating the TargetPoint of the corresponding kernel)
    kernelsize = cell([1 nKernels]);               % each cell element is 1x4-matrix, with elements [x_up x_down y_left y_right]
    kernel = cell([1 nKernels]);
    kernel_TargetPoint = cell([1 nKernels]);
    SrcRelativeTarg = cell([1 nKernels]);
    weights = cell([1 nKernels]);


    % Iterate over all Kernels
	for KernelIndex = 1:nKernels
		
		if(~quiet_flag)
			fprintf('\nFilling Glass %d',KernelIndex)
		end
		
        % Compute kernel size for processed Kernel
        kernelsize{KernelIndex} = [1 1 1 1];
        no_SourcePoints = 0;
		while(no_SourcePoints < MinKernelSrcPts)
            kernelsize{KernelIndex} = kernelsize{KernelIndex} + 1;

            kernel_dummy = UndersamplingPattern_ACS(CentralCellMatIndex_ACS_x(KernelIndex) - kernelsize{KernelIndex}(1): ...
                                                    CentralCellMatIndex_ACS_x(KernelIndex) + kernelsize{KernelIndex}(2), ...
                                                    CentralCellMatIndex_ACS_y(KernelIndex) - kernelsize{KernelIndex}(3): ...
                                                    CentralCellMatIndex_ACS_y(KernelIndex) + kernelsize{KernelIndex}(4));

            no_SourcePoints = sum(sum(kernel_dummy));
		end


		% Test if the kernel size can be reduced in one direction without losing any Source points
		kernel_dummy_sumrow = sum(kernel_dummy,1);
		kernel_dummy_sumcol = sum(kernel_dummy,2);
		KernelCropCellLog = ~[kernel_dummy_sumcol(1),kernel_dummy_sumcol(end),kernel_dummy_sumrow(1),kernel_dummy_sumrow(end)];
		KernelCropCellDummy = {[1 0 0 0],[0 1 0 0],[0 0 1 0],[0 0 0 1]};
		KernelCropCell = KernelCropCellDummy(KernelCropCellLog);
	
		for kernelcrop = KernelCropCell

			for crop_repeat = 1:max(size(UndersamplingCell))    % In GRAPPA-like kernels, one has to remove several lines/rows, if R_x > 2 | R_y > 2
                kernel_dummy_cropped = kernel_dummy(1 + kernelcrop{1}(1) : end - kernelcrop{1}(2), 1 + kernelcrop{1}(3) : end - kernelcrop{1}(4));

                no_SourcePoints_cropped = sum(sum(kernel_dummy_cropped));

				if(no_SourcePoints_cropped == no_SourcePoints)
                    kernel_dummy = kernel_dummy_cropped;
                    kernelsize{KernelIndex} = kernelsize{KernelIndex} - kernelcrop{1};
				else
					break;
				end
			end

		end
        clear kernel_dummy_cropped no_SourcePoints_cropped



        % Test for kernel equality. If the now defined kernel is equal to a formerly processed one, don't do the rest, but copy the weights of this former kernel
        kernel_found = false;
        for test_against_kernel_no = 1:KernelIndex-1
            if(isequal(kernelsize{test_against_kernel_no},kernelsize{KernelIndex}))
                if(isequal(squeeze(kernel{test_against_kernel_no}(1,:,:)),kernel_dummy))        % kernel gets replicated with channels. --> Must index into kernel.
                    % Copy Weights of kernel{test_against_kernel_no}
                    [weights{KernelIndex}] = deal(weights{test_against_kernel_no});
                    SrcRelativeTarg{KernelIndex} = SrcRelativeTarg{test_against_kernel_no};
					KernelCorrespondence(KernelIndex) = test_against_kernel_no;
                    kernel_found = true;
                    break
                end
            end
        end
        %kernel{KernelIndex} = kernel_dummy;    
        kernel{KernelIndex} = myrepmat_1_0(logical(kernel_dummy), [nChannel size(kernel_dummy)], 2);
        kernel_TargetPoint{KernelIndex} = false(size(kernel_dummy));
        kernel_TargetPoint{KernelIndex}(kernelsize{KernelIndex}(1)+1,kernelsize{KernelIndex}(3)+1) = true;
        %kernel_TargetPoint{KernelIndex} = myrepmat_1_0(kernel_TargetPoint{KernelIndex}, [nChannel size(kernel_TargetPoint{KernelIndex})], 2);
        if(kernel_found)
            continue
        end


        [SrcRelativeTarg_x,SrcRelativeTarg_y] = find(kernel_dummy);
        [TargPos_x, TargPos_y] = find(kernel_TargetPoint{KernelIndex});
        SrcRelativeTarg{KernelIndex} = [SrcRelativeTarg_x-TargPos_x,SrcRelativeTarg_y-TargPos_y];
        clear SrcRelativeTarg_x SrcRelativeTarg_y TargPos_x TargPos_y






        % Compute the Source and Target Points as linear indices
        % What does the code do here? Let's assume ... (ohh... did I fall asleep? Hm.)
        % No, really: The source and the target points are computed for the ACS data. This is done by computing the linear indices of both
        % (to avoid slow loops). The computation of the linear indices looks quite complicated, and in fact -- it is :)
        % Well in principle the target points are just all points all x- and y-values of all channels. The relative distance of the
        % source points to a target point is computed within the kernel, and this relative information is applied to the target points
        % in order to get the source points. Then the linear indices can be computed.


		AdditionalReplication.DimPos = 2; AdditionalReplication.AddToDims = [2 3]; AdditionalReplication.Mat = int16(SrcRelativeTarg{KernelIndex});
		Source_linear = sub2ind_extended([nChannel nx_ACS ny_ACS], ...
			             AdditionalReplication,0, ...
						 1:nChannel, kernelsize{KernelIndex}(1)+1:nx_ACS-kernelsize{KernelIndex}(2), kernelsize{KernelIndex}(3)+1:ny_ACS-kernelsize{KernelIndex}(4));
		Target_linear = sub2ind_extended([nChannel nx_ACS ny_ACS], ...
			             0,0, ...
						 1:nChannel, kernelsize{KernelIndex}(1)+1:nx_ACS-kernelsize{KernelIndex}(2), kernelsize{KernelIndex}(3)+1:ny_ACS-kernelsize{KernelIndex}(4));

		if(max(Source_linear) > 2^31)
            max(Source_linear)
            display('Change to uint64 in code. Aborting.')
            weights = 0;
            return
		end

		weights{KernelIndex} = zeros([nChannel nChannel*no_SourcePoints nSlice_ACS]);
        for SliceIndex = 1:nSlice_ACS

            % Slice ACS data
            ACS_sliced = ACS(:,:,:,SliceIndex);

            % The Target Points, size: nChannel x kernel repetitions in ACS data
            TargetPoints_ACS = ACS_sliced(Target_linear);
            TargetPoints_ACS = reshape(TargetPoints_ACS, [nChannel numel(TargetPoints_ACS)/nChannel]);

            % The Source Points, size: nChannel*no_SourcePoints x kernel repetitions in ACS data
            SourcePoints_ACS = ACS_sliced(Source_linear);    
            SourcePoints_ACS = reshape(SourcePoints_ACS, [nChannel*no_SourcePoints numel(SourcePoints_ACS)/(nChannel*no_SourcePoints)]);



            % find weights by fitting the source data to target data.
            % The pinv averages over all kernel repetitions in a weighted way ('least square' solution)
            % size: nChannel x nChannel*no_SourcePoints
            weights{KernelIndex}(:,:,SliceIndex) = TargetPoints_ACS * pinv(SourcePoints_ACS); 

        end



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


% Create an extended matrix, extended by the maximum kernelsize. This extension is necessary, so that also the target points at the border of the matrix can be reconstructed, using the kernel.
% Write undersampled data, and zero pad this data

% Find Maximum Kernelsize of all kernels
Maxkernelsize = [kernelsize{:}];
Maxkernelsize = max(reshape(Maxkernelsize, [4 numel(Maxkernelsize)/4]), [], 2);

% Prepare the enlarged/extended Reconstruction Matrix
Reco_dummy = zeros([nChannel nx+sum(Maxkernelsize(1:2)) ny+sum(Maxkernelsize(3:4)) nSlice nTime]);              % This has to be changed for "TODO 2)"
Reco_dummy(:,Maxkernelsize(1)+1:end-Maxkernelsize(2),Maxkernelsize(3)+1:end-Maxkernelsize(4),:,:) = OutData;
Reco_dummy_CheckIfRecoNecess = squeeze(abs(Reco_dummy(1,:,:,1,1)) > 0);
%clear OutData;
%fprintf('\nInitializing took %4.2f s',toc);


for KernelIndex = 1:nKernels
	if(~quiet_flag)
		fprintf('\nDrinking Glass %d\tGUUULLPP . . . UAHHHH . . . ',KernelIndex)
	end
    
	if(KernelCorrespondence(KernelIndex) > 0)
		continue;
	end
	
    % Compute all the linear indices of the target points within Reco_dummy of the processed Kernel
    UndersamplingCell_CurrentTargetPoint = false(size(UndersamplingCell));											% Mark the current target point in the undersampling cell
    UndersamplingCell_CurrentTargetPoint(ElementaryCellLinearIndex(KernelIndex)) = true; 
	UndersamplingCell_CurrentTargetPoint(ElementaryCellLinearIndex(KernelIndex == KernelCorrespondence)) = true;	% And also the points with the same weights
    TargetPoints = false([size(Reco_dummy,2) size(Reco_dummy,3)]);
    
    
    % This has to be changed for "TODO 2)"
    % What is done here?
    % 1) The inner part of the TargetPoints are called, because the target points can only be in that part of Reco_dummy which implements the InData, not the zerofilled part 
    %    (we dont want to reconstruct the zeros, do we?)
    % 2) The UndersamplinCell with the current target point gets replicated to the size of the InData. If the InData is not an integer multiple of the UndersamplingCell then
    %    this will lead to an error here.
    % 3) Then the x- and y-coordinates of the target points are computed.
    UC_dummy = repmat(UndersamplingCell_CurrentTargetPoint, [floor(nx/size(UndersamplingCell,1)) floor(ny/size(UndersamplingCell,2))]);
    UC_dummy = cat(  1,UC_dummy, UC_dummy(1:nx-size(UC_dummy,1),:,:)  );
    UC_dummy = cat(  2,UC_dummy, UC_dummy(:,1:ny-size(UC_dummy,2),:)  );
    %UC_dummy = cat(  3,UC_dummy, UC_dummy(:,:,1:nSlice-size(UC_dummy,3))  );  
    TargetPoints(Maxkernelsize(1)+1:end-Maxkernelsize(2),Maxkernelsize(3)+1:end-Maxkernelsize(4),:) = UC_dummy;
%     TargetPoints(Maxkernelsize(1)+1:end-Maxkernelsize(2),Maxkernelsize(3)+1:end-Maxkernelsize(4),:) = ...
%     repmat(UndersamplingCell_CurrentTargetPoint, [floor(nx/size(UndersamplingCell,1)) floor(ny/size(UndersamplingCell,2)) floor(nSlice/size(UndersamplingCell,3))]);  
    [TargetPoints_x, TargetPoints_y] = find(TargetPoints);
     	
	
	
    % Loop over all target points for the processed kernel
	SourcePoints = zeros([nChannel*size(SrcRelativeTarg{KernelIndex},1) nTime]);	
	for Targetloopy = 1:numel(TargetPoints_x)
        SourceIndices = [TargetPoints_x(Targetloopy) + SrcRelativeTarg{KernelIndex}(:,1), TargetPoints_y(Targetloopy) + SrcRelativeTarg{KernelIndex}(:,2)];
		
		if(sum(Reco_dummy_CheckIfRecoNecess(sub2ind(size(Reco_dummy_CheckIfRecoNecess),SourceIndices(:,1),SourceIndices(:,2)))) < 2)
			continue
		end
		
		
		
		for SliceIndex = 1:nSlice
            % Compute the Source Points for that specific TargetIndex
			for Sourceloopy = 1:size(SourceIndices,1)
                SourcePoints((Sourceloopy-1)*nChannel+1:Sourceloopy*nChannel,:) = Reco_dummy(:,SourceIndices(Sourceloopy,1),SourceIndices(Sourceloopy,2),SliceIndex,:);
			end
            % Reconstruct data by applying weights to Sourcepoints.
            Reco_dummy(:,TargetPoints_x(Targetloopy),TargetPoints_y(Targetloopy),SliceIndex,:)=reshape(weights{KernelIndex}(:,:,SliceIndex)*SourcePoints, [nChannel nTime]); 
		end
	end

end


OutData = Reco_dummy(:,Maxkernelsize(1)+1:end-Maxkernelsize(2),Maxkernelsize(3)+1:end-Maxkernelsize(4),:,:);                                           %Crop out the good data.

if(~quiet_flag)
	fprintf('... %f sec \n',toc)
end






%% ~. Postparations

% Fancy Text Message
if(~quiet_flag)
	fprintf('Mmmmhhh, I love CAIPIRINHA! \n')    
end


