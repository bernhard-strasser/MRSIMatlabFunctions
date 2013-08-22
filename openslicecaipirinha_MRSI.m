function [OutData,weights]=openslicecaipirinha_MRSI(InData, ACS, FoV_shifts, kernelsize) 
% 
% openslicegrappa_MRSI Reconstruct the Slices of MRSI and MRI Data Ehen Only the Sum of Those Slices Was Measured.
% 
%  [OutData,weights] = openslicegrappa_MRSI(InData,ACS,kernelsize);
%       
%   Input:      
% -     InData             Undersampled Data            (size: [#coils, nx, ny, nSlice = 1, nTime]) (nTime = 1 for MRI)
%                          There must be only 1 slice!
% -     ACS                AutoCalibration Signal       (size: [#coils, nx_ACS, ny_ACS, nSlice_ACS])  
%                          Must have as many slices as the OutData should have. Use similar sequence parameters as for the InData.
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

% Fancy Text Message
fprintf('\nLet the lemon-slice-caipirinha party start!')


% Assign standard values to variables if nothing is passed to function.

if(~exist('InData','var'))
    display([char(10) 'I need Cachaca ( = data which should be processed (InData)) for preparing CAIPIRINHA! Aborting . . .'])
    return
end

if(~exist('ACS','var'))
    display([char(10) 'I need Glasses ( = Auto Calibration Signal (ACS)) for preparing CAIPIRINHA! Aborting . . .'])
    return
end
if(~exist('FoV_shifts','var'))
    FoV_shifts = zeros([size(ACS,4) 2]);
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
    disp([char(10) 'Error! The number of coils has to be the same for both inputs! Aborting . . .'])
    OutData = InData;
    weights = 0;
    return;
end
if(nSlice == nSlice_ACS)
    disp([char(10) 'Nothing to do, since the InData has the same amount of slices as the ACS data . . .'])
    OutData = InData;
    weights = 0;
    return;
end








%% 1. Shift the ACS FoV


ACS = kSpace_FoVShift(ACS,FoV_shifts);




%% 2. Calculate weights

% Fancy Text Message
tic;
fprintf('\nFilling the glasses and cutting lemon slices (Calculating weights)')

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
tic; fprintf('\nDrinking Caipirinha with slices of lemon (Applyig weights).')



% Create an extended matrix, extended by the kernelsize. This extension is necessary, so that also the target points at the border of the matrix can be reconstructed, using the kernel.
Reco_dummy = zeros([nChannel nx+sum(kernelsize(1:2)) ny+sum(kernelsize(3:4)) 1 nTime]);
Reco_dummy(:,kernelsize(1)+1:end-kernelsize(2),kernelsize(3)+1:end-kernelsize(4),1,:) = InData;
InPlane_UndersamplingPattern = abs(squeeze(Reco_dummy(1,:,:,1,1))) > 0;
[TargetPoints_x TargetPoints_y] = find(InPlane_UndersamplingPattern);   % Find the target points that should be reconstructed.


% Create the OutData which has to have as many slices as the ACS data because we want to reconstruct all those slices.
OutData = zeros([nChannel nx ny nSlice_ACS nTime]);


for SliceIndex = 1:nSlice_ACS
    fprintf('\nDrinking Caipirinha with lemon slice %d\tGUUULLPP . . . UAHHHH . . . ',SliceIndex)

    
    % Loop over all target points for the processed kernel
    for Targetloopy = 1:numel(TargetPoints_x)

        % Sourcepoints are the points of all channels and all x and y points within the kernel around the target point.
        SourcePoints = reshape(Reco_dummy(:,TargetPoints_x(Targetloopy)-kernelsize(1):TargetPoints_x(Targetloopy)+kernelsize(2), ...
        TargetPoints_y(Targetloopy)-kernelsize(3):TargetPoints_y(Targetloopy)+kernelsize(4),1,:), [nChannel*kernelsize_x*kernelsize_y nTime]);

        % Reconstruct data by applying weights to Sourcepoints.
        OutData(:,TargetPoints_x(Targetloopy)-kernelsize(1),TargetPoints_y(Targetloopy)-kernelsize(3),SliceIndex,:)=reshape(weights{SliceIndex}*SourcePoints, [nChannel 1 1 1 nTime]); 
                
    end


end



fprintf('\nDrinking took\t...\t%f sec.',toc)




%% 4. Undo FoV-Shifts

OutData = kSpace_FoVShift(OutData,-FoV_shifts);




%% ~. Postparations

% Fancy Text Message
fprintf('\nMmmmhhh, I love Caipirinha! \n')    



