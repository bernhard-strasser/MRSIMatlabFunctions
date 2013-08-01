function [OutData,weights]=opencaipirinha_MRSI_0_2(InData, ACS, UndersamplingCell, MinKernelSrcPts) 
% 
% opengrappa_MRSI_x_y Reconstruct Undersampled MRSI and MRI Data
% 
%  [OutData,weights] = opengrappa(InData,ACS,UndersamplingCell);
%       
%   Input:      
% -     InData             Undersampled Data            (size: [#coils, k_x, k_y, slices, nTime]) (nTime = 1 for MRI)
%                          The undersampled points must be set to zero, NOT MISSING!
% -     ACS                AutoCalibration Signal       (size: [#coils, k_y, k_x, slices])     
% -     UndersamplingCell  Elementary Cell, logical array (or array containing only 1's and 0's) which tells you the measured points (1's) and omitted points (0's). 
%                          Gets replicated to spatial size of InData to define the undersampling pattern. 
%                          Thus spatial size of InData and ACS must be integer multiple of UndersamplingCell.
% -     MinKernelSrcPts    The minimum source points in the kernel that are fitted to the target points. 20 is a good value, as that is standard in GRAPPA.
% 
%   Output:
% -     OutData            Reconstructed Output Data    (size: [#coils, k_x, k_y, slc, nTime]) (nTime = 1 for MRI)
% -     weights            Weights for Reconstruction              
%              
%



% This code is based on the teaching version of the opengrappa function of Felix Breuer that was provided 
% on the parallel imaging workshop in Wuerzburg, 2012, which in turn is based on the opengrappa of Mark Griswold.
%
% 22.06.2008 Felix Breuer (breuer@mr-bavaria.de)
% November 2012 - January 2013 Bernhard Strasser
%
%
%   Please read the license text at the bottom of this program. By using this program, you 
%   implicity agree with the license. 
%
%   The main points of the license:
%
%   1) This code is strictly for non-commercial applications. The code is protected by
%      multiple patents.
%   2) This code is strictly for research purposes, and should not be used in any
%      diagnostic setting.




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
%       1)   Get License.
%       2)   If the InData is not integer muiltiple of UndersamplingCell, it should still be processed:
%            When the data to reconstruct get zerofilled in order to apply the kernel also to the border voxels,
%            copy this extent data instead of the zeros to the zerofilled data.
%       3)   Can something similar be achieved for the ACS data ???
%       4)   If it is not too complicated, compute different kernels for the border voxels. 
%            Because the border voxels are computed with a lot of zeros, but with the "normal" reconstruction weights.
%            Yet, these "normal" weights assume that there is still data in all source points, and not zeros.
%       5)   A ball-shaped kernel instead of a cube might be a little better. Only little advantages expected.
%       6)   Instead of gathering the Source & Target Points for computing the weights in a loop, do it probably
%            by creating a Target-mask and a Source-mask.


%% 0. Preparation

% Assign standard values to variables if nothing is passed to function.

if(~exist('InData','var'))
    display([char(10) 'You should consider inputting data which can be reconstructed. Aborting . . .'])
    return
end

if(~exist('ACS','var'))
    display([char(10) 'I need an Auto Calibration Signal (ACS) for drinking GRAPPA/CAIPIRINHA! Aborting . . .'])
    return
end

if(~exist('ACS','var'))
    display([char(10) 'Please tell me how your data was undersampled by inputting ''UndersamplingCell''. Aborting . . .'])
    return
end

if(~exist('MinKernelSrcPts','var'))
    MinKernelSrcPts = 20;
end





% Further Preparations

% Get the size of both the input data and the autocalibration data
[nChannel,nx,ny, nSlice, nTime] = size(InData);
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





%% 2. Calculate weights

% Fancy Text Message
tic;
fprintf('\nFilling the glasses (Calculating weights)')

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
% For visualization, run the script Caipirinha_VisualizingWeightsComputation_And_Application_1_x.m




% Create the undersampling pattern of the ACS data, to determine the source and target points for computing the weights
UndersamplingPattern_ACS = repmat(UndersamplingCell, [floor(nx_ACS/size(UndersamplingCell,1)) floor(ny_ACS/size(UndersamplingCell,2)) floor(nSlice_ACS/size(UndersamplingCell,3))]);

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
[CentralCellMatIndex_ACS_x CentralCellMatIndex_ACS_y] = ind2sub(size(UndersamplingPattern_ACS), CentralCellLinearIndex_ACS);



% Preallocate kernelsize, kernel and the kernel_TargetPoint (indicating the TargetPoint of the corresponding kernel)
kernelsize = cell([1 nKernels]);               % each cell element is 1x4-matrix, with elements [x_up x_down y_left y_right]
kernel = cell([1 nKernels]);
kernel_TargetPoint = cell([1 nKernels]);
weights = cell([1 nKernels]);

% Iterate over all Kernels
for KernelIndex = 1:nKernels

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
    for kernelcrop = {[1 0 0 0], [0 1 0 0], [0 0 1 0], [0 0 0 1]}
        
        for crop_repeat = 1:max(size(UndersamplingCell))    % In GRAPPA-like kernels, one has to remove several lines/rows, if R_x > 2 | R_y > 2
            kernel_dummy_cropped = kernel_dummy(1 + kernelcrop{1}(1) : end - kernelcrop{1}(2), 1 + kernelcrop{1}(3) : end - kernelcrop{1}(4));

            no_SourcePoints_cropped = sum(sum(kernel_dummy_cropped));

            if(no_SourcePoints_cropped == no_SourcePoints)
                kernel_dummy = kernel_dummy_cropped;
                kernelsize{KernelIndex} = kernelsize{KernelIndex} - kernelcrop{1};
            end
        end
        
    end
    clear kernel_dummy_cropped no_SourcePoints_cropped
    
  
    
    % Test for kernel equality. If the now defined kernel is equal to a formerly processed one, don't do the rest, but copy the weights of this former kernel
    kernel_found = false;
    for test_against_kernel_no = 1:KernelIndex-1
        if(isequal(kernelsize{test_against_kernel_no},kernelsize{KernelIndex}))
            if(isequal(kernel{test_against_kernel_no},kernel_dummy))
                % Copy Weights of kernel{test_against_kernel_no}
                weights{KernelIndex} = weights{test_against_kernel_no};
                kernel_found = true;
                break
            end
        end
    end
    %kernel{KernelIndex} = kernel_dummy;    
    kernel{KernelIndex} = myrepmat_1_0(logical(kernel_dummy), [nChannel size(kernel_dummy)], 2);
    kernel_TargetPoint{KernelIndex} = false(size(kernel_dummy));
    kernel_TargetPoint{KernelIndex}(kernelsize{KernelIndex}(1)+1,kernelsize{KernelIndex}(3)+1) = true;
    kernel_TargetPoint{KernelIndex} = myrepmat_1_0(kernel_TargetPoint{KernelIndex}, [nChannel size(kernel_TargetPoint{KernelIndex})], 2);
    if(kernel_found)
        break
    end
    
    
    
    
    
    % Gather all the Source and Target Points for the fitting process / computing the weights
    SourcePoints_ACS = zeros([sum(sum(sum(sum(kernel{KernelIndex})))) (ny_ACS-kernelsize{KernelIndex}(4)-kernelsize{KernelIndex}(3)) * (nx_ACS-kernelsize{KernelIndex}(1)-kernelsize{KernelIndex}(2))]);
    TargetPoints_ACS = zeros([nChannel (ny_ACS-kernelsize{KernelIndex}(4)-kernelsize{KernelIndex}(3)) * (nx_ACS-kernelsize{KernelIndex}(1)-kernelsize{KernelIndex}(2))]);  
    loopy = 0;                          %This is a very lazy counter. Could be done much faster. 
    for yind=kernelsize{KernelIndex}(3)+1:ny_ACS-kernelsize{KernelIndex}(4)
        for xind=kernelsize{KernelIndex}(1)+1:nx_ACS-kernelsize{KernelIndex}(2)
            loopy=loopy+1;
            

            SourcePattern_ACS = false([nChannel size(UndersamplingPattern_ACS)]);
            SourcePattern_ACS(:,xind-kernelsize{KernelIndex}(1):xind+kernelsize{KernelIndex}(2),yind-kernelsize{KernelIndex}(3):yind+kernelsize{KernelIndex}(4)) = kernel{KernelIndex};
            
            
%             TargetPattern_ACS = false([nChannel size(UndersamplingPattern_ACS)]);
%             TargetPattern_ACS(:,xind-kernelsize{KernelIndex}(1):xind+kernelsize{KernelIndex}(2),yind-kernelsize{KernelIndex}(3):yind+kernelsize{KernelIndex}(4)) = kernel_TargetPoint{KernelIndex};
            
            
            % These are the source points. size: #coils*kernelsize_y*kernelsize_x   x   kernel repetitions in ACS data
            % In the first index just all the ACS points which should be fitted to the target are saved, for one kernel repetition
            % In the second index, all kernel repetitions are saved.
            SourcePoints_ACS(:,loopy) = ACS(SourcePattern_ACS);                   


            % these are the taget points. size: #coils  x  kernel repetitions in ACS data    
            % Again, the first index gives all the target points.
            % The second index are the different target points of the different kernel repetitions within the ACS data.
            TargetPoints_ACS(:,loopy) = ACS(:,xind,yind,1);

        end
    end
   
    % find weights by fitting the source data to target data.
    % The pinv seems to average over all kernel repetitions in a weighted way ('least square' solution)
    weights{KernelIndex} = TargetPoints_ACS * pinv(SourcePoints_ACS); 


    fprintf('... %f sec \n',toc)                                                          
    fprintf('Aaaaahhh! This SMELL! \n')    
    
    
end





%% 3. Apply weights

% Fancy Text Message
tic; fprintf('Drinking CAIPIRINHA (Applyig weights) ')



% Create an extended matrix, extended by the maximum kernelsize. This extension is necessary, so that also the target points at the border of the matrix can be reconstructed, using the kernel.
% Write undersampled data, and zero pad this data

% Find Maximum Kernelsize of all kernels
Maxkernelsize = [kernelsize{:}];
Maxkernelsize = max(transpose(reshape(Maxkernelsize, [4 numel(Maxkernelsize)/4])), [], 1);

% Prepare the enlarged/extended Reconstruction Matrix
Reco_dummy = zeros([nChannel nx+sum(kernelsize{KernelIndex}(1:2)) ny+sum(kernelsize{KernelIndex}(3:4)) nSlice nTime]);              % This has to be changed for "TODO 2)"
Reco_dummy(:,kernelsize{KernelIndex}(1)+1:end-kernelsize{KernelIndex}(2),kernelsize{KernelIndex}(3)+1:end-kernelsize{KernelIndex}(4),:,:) = InData; 



for KernelIndex = 1:nKernels
    fprintf('Gulp ')

    % Computer number of source points for that kernel
    no_SourcePoints = sum(sum(sum(kernel{KernelIndex})));
    
    % Compute all the linear indices of the target points within Reco_dummy of the processed Kernel
    UndersamplingCell_CurrentTargetPoint = false(size(UndersamplingCell));                  % Mark the current target point in the undersampling cell
    UndersamplingCell_CurrentTargetPoint(ElementaryCellLinearIndex(KernelIndex)) = true;       
    TargetPoints = repmat(UndersamplingCell_CurrentTargetPoint, [floor(nx/size(UndersamplingCell,1)) floor(ny/size(UndersamplingCell,2)) floor(nSlice/size(UndersamplingCell,3))]);  % This has to be changed for "TODO 2)"
    TargetPoints = cat(1,zeros(Maxkernelsize(1),size(TargetPoints,2)), TargetPoints, zeros(Maxkernelsize(2),size(TargetPoints,2)));
    TargetPoints = cat(2,zeros(size(TargetPoints,1),Maxkernelsize(3)), TargetPoints, zeros(size(TargetPoints,1),Maxkernelsize(4)));
    [TargetPoints_x TargetPoints_y] = find(TargetPoints);
    %TargetPattern = myrepmat_1_0(TargetPattern, [nSlice size(TargetPattern) nSlice nTime], 2); 
    %TargetPattern = transpose(find(TargetPattern));
     
    
    % Loop over all target points for the processed kernel
    for Targetloopy = 1:numel(TargetPoints_x)
        display([num2str(Targetloopy) ' of ' num2str(numel(TargetPoints_x))])

        % Compute the Source Points for that specific TargetIndex
        SourcePattern = false(size(Reco_dummy));
        kernel_dummy = repmat(kernel{KernelIndex}, [1 1 1 nSlice nTime]);
        SourcePattern(:,TargetPoints_x(Targetloopy)-kernelsize{KernelIndex}(1):TargetPoints_x(Targetloopy)+kernelsize{KernelIndex}(2), ...
                      TargetPoints_y(Targetloopy)-kernelsize{KernelIndex}(3):TargetPoints_y(Targetloopy)+kernelsize{KernelIndex}(4),:,:) = kernel_dummy;
        
        SourcePoints=reshape(Reco_dummy(SourcePattern), [no_SourcePoints nSlice*nTime]);        
        
        % Reconstruct data by applying weights to Sourcepoints.
        Reco_dummy(:,TargetPoints_x(Targetloopy),TargetPoints_y(Targetloopy),:,:)=reshape(weights{KernelIndex}*SourcePoints, [nChannel 1 1 nSlice nTime]); 
    end




end


OutData = Reco_dummy(:,Maxkernelsize(1)+1:nx+Maxkernelsize(2),Maxkernelsize(3)+1:ny+Maxkernelsize(4),:,:);                                           %Crop out the good data.

fprintf('... %f sec \n',toc)




%% ~. Postparations

% Fancy Text Message
fprintf('Mmmmhhh, I love CAIPIRINHA! \n')    



