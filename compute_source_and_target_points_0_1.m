function [Source_linear,Target_linear]=compute_source_and_target_points_0_1(size_RefData, kernel,kernelsize UndersamplingCell, MinKernelSrcPts) 
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





%% 0. Preparation

% Assign standard values to variables if nothing is passed to function.

if(~exist('InData','var'))
    display([char(10) 'You should consider inputting data which can be reconstructed. Aborting . . .'])
    return
end


[nChannel,nx,ny,nz,nTime] = size_RefData;

    
%% 1. Compute Source Points

% Compute the Source and Target Points as linear indices
nx_wo_border = nx - sum(kernelsize(1:2));
ny_wo_border = ny - sum(kernelsize(3:4));

% Source Points
% Create channel info
Source_Channels = int16(transpose(1:nChannel));
Source_Channels = repmat(Source_Channels, [1 no_SourcePoints nx_wo_border ny_wo_border]);

% Create spatial info
Source_x = int16(kernelsize(1)+1:nx-kernelsize(2));
Source_y = int16(kernelsize(3)+1:ny-kernelsize(4));

% Copy Source to Target Points
Target_x = Source_x;
Target_y = Source_y;

% Apply Relative Info
Source_x = repmat(transpose(Source_x), [1 no_SourcePoints]);
Source_y = repmat(transpose(Source_y), [1 no_SourcePoints]);
Source_x = Source_x + repmat(reshape(int16(SrcRelativeTarg{1}(:,1)),[1 no_SourcePoints]),[size(Source_x,1) 1]);
Source_y = Source_y + repmat(reshape(int16(SrcRelativeTarg{1}(:,2)),[1 no_SourcePoints]),[size(Source_y,1) 1]);

% Replicate spatial info
Source_x = repmat(Source_x, [1 1 nChannel ny_wo_border]);
Source_y = repmat(reshape(Source_y, [ny_wo_border no_SourcePoints]), [1 1 nChannel nx_wo_border]);

% Reorder Source Points
Source_x = permute(Source_x, [3 2 1 4]);
Source_y = permute(Source_y, [3 2 4 1]);


% Target Points
Target_Channels = reshape(Source_Channels(:,1,:,:), [nChannel nx_wo_border ny_wo_border]);
Target_x = repmat(Target_x,[nChannel 1 ny_wo_border]);
Target_y = repmat(reshape(Target_y,[1 1 numel(Target_y)]),[nChannel nx_wo_border 1]);


% Linear Indices
Target_linear = sub2ind([nChannel nx ny], uint32(reshape(Target_Channels, [1 numel(Target_Channels)])), uint32(reshape(Target_x, [1 numel(Target_x)])), uint32(reshape(Target_y, [1 numel(Target_y)])));    
Source_linear = sub2ind([nChannel nx ny], uint32(reshape(Source_Channels, [1 numel(Source_Channels)])), uint32(reshape(Source_x, [1 numel(Source_x)])), uint32(reshape(Source_y, [1 numel(Source_y)])));
% For Reconstructing MRSI data, uint32 has to be changed to uint64
    
    
    

    
    
    


%% ~. Postparations

    



