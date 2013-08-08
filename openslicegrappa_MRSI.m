function [OutData,weights]=openslicegrappa_MRSI(InData, ACS, kernelsize) 
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
%       4)   If it is not too complicated, compute different kernels for the border voxels. 
%            Because the border voxels are computed with a lot of zeros, but with the "normal" reconstruction weights.
%            Yet, these "normal" weights assume that there is still data in all source points, and not zeros.
%       5)   A ball-shaped kernel instead of a cube might be a little better. Only little advantages expected.



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

if(~exist('kernelsize','var'))
    kernelsize = [2 2 2 2];
end
if(numel(kernelsize)==1)
    kernelsize = repmat(kernelsize(1),[1 4]);
elseif(numel(kernelsize) < 4)
    kernelsize = [floor(kernelsize(1)/2) floor(kernelsize(1)/2) floor(kernelsize(2)/2) floor(kernelsize(2)/2)];
end




% Further Preparations


% kernelsize
kernelsize_x = sum(kernelsize(1:2))+1;
kernelsize_y = sum(kernelsize(3:4))+1;

% Get the size of both the input data and the autocalibration data
[nChannel,nx,ny, nSlice, nTime] = size(InData);
[nChannel_ACS,nx_ACS,ny_ACS, nSlice_ACS]=size(ACS);

% Check for dimension size
if(nChannel_ACS~=nChannel)
    disp('Error! The number of coils has to be the same for both inputs! Aborting . . .')
    OutData = InData;
    weights = 0;
    return;
end
if(nSlice == nSlice_ACS)
    disp('Nothing to do, since the InData has the same amount of slices as the ACS data . . .')
    OutData = InData;
    weights = 0;
    return;
end







%% 2. Calculate weights

% Fancy Text Message
tic;
fprintf('\nFilling the glasses (Calculating weights)')

% What does the code do here? Let's assume ... (ohh... did I fall asleep? Hm.)
%
% For visualization, run the script Caipirinha_VisualizingWeightsComputation_And_Application_1_x.m


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


% The Source Points, size: nChannel*no_SrcPts x kernel repetitions in ACS data
ACS_sum = sum(ACS,4);
SourcePoints_ACS = ACS_sum(Source_linear);
clear ACS_sum
SourcePoints_ACS = reshape(SourcePoints_ACS, [nChannel*no_SrcPts numel(SourcePoints_ACS)/(nChannel*no_SrcPts)]);
    
 
% Iterate over all Slices
weights = cell([1 nSlice_ACS]);
for SliceIndex = 1:nSlice_ACS
    
    fprintf('\nCutting slice %d of the lemon.',SliceIndex)
    
    % The Target Points, size: nChannel x kernel repetitions in ACS data    
    ACS_sliced = ACS(:,:,:,SliceIndex);    
    TargetPoints_ACS = ACS_sliced(Target_linear);
    TargetPoints_ACS = reshape(TargetPoints_ACS, [nChannel numel(TargetPoints_ACS)/nChannel]);    
    
    % find weights by fitting the source data to target data.
    % The pinv averages over all kernel repetitions in a weighted way ('least square' solution)
    % size: nChannel x nChannel*no_SrcPts
    weights{SliceIndex} = TargetPoints_ACS * pinv(SourcePoints_ACS); 
    
end

fprintf('\nCutting took\t... %f sec \n',toc)                                                          
fprintf('Aaaaahhh! This SMELL! \n') 





%% 3. Apply weights

% Fancy Text Message
tic; fprintf('\nDrinking GRAPPA with slices of lemon (Applyig weights).')



% Create an extended matrix, extended by the kernelsize. This extension is necessary, so that also the target points at the border of the matrix can be reconstructed, using the kernel.
Reco_dummy = zeros([nChannel nx+sum(kernelsize(1:2)) ny+sum(kernelsize(3:4)) 1 nTime]);
Reco_dummy(:,kernelsize(1)+1:end-kernelsize(2),kernelsize(3)+1:end-kernelsize(4),1,:) = InData;

% Create the OutData which has to have as many slices as the ACS data because we want to reconstruct all those slices.
OutData = zeros([nChannel nx ny nSlice_ACS nTime]);


for SliceIndex = 1:nSlice_ACS
    fprintf('\nDrinking Grappa with lemon slice %d\tGUUULLPP . . . UAHHHH . . . ',SliceIndex)

    
    % Loop over all target points for the processed kernel
    for Target_y = kernelsize(3)+1:ny+kernelsize(3)
        for Target_x = kernelsize(1)+1:nx+kernelsize(1)            

            % Sourcepoints are the points of all channels and all x and y points within the kernel around the target point.
            SourcePoints = reshape(Reco_dummy(:,Target_x-kernelsize(1):Target_x+kernelsize(2),Target_y-kernelsize(3):Target_y+kernelsize(4),1,:), [nChannel*kernelsize_x*kernelsize_y nTime]);

            % Reconstruct data by applying weights to Sourcepoints.
            OutData(:,Target_x-kernelsize(1),Target_y-kernelsize(3),SliceIndex,:)=reshape(weights{SliceIndex}*SourcePoints, [nChannel 1 1 1 nTime]); 

        end
    end


end



fprintf('\nDrinking took\t...\t%f sec.',toc)







%% ~. Postparations

% Fancy Text Message
fprintf('\nMmmmhhh, I love CAIPIRINHA! \n')    



