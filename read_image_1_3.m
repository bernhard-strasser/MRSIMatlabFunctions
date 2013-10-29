function [image, image_kspace] = read_image_1_3(image_path,DesiredSize,interpol_method,flip,fredir_shift,Hamming_flag,EllipticalFilterSize, phase_encod_dir, sum_averages_flag)
%
% read_image_x_x Read in csi-data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function only decides if input file is .dat file or DICOM according to its ending. It then calls the read_image_dat_x_x or 
% read_image_dicom_x_x functions. Refer to these for more info.
%
%
% [csi,csi_kspace] = read_image_1_2(image_path,DesiredSize,interpol_method,flip,fredir_shift,Hamming_flag,EllipticalFilterSize, phase_encod_dir)
%
% Input: 
% -         image_path                  ...     Path of MRI file.
% -         DesiredSize                 ...     [desired size in freq encod dir, desired size in phase encod dir] 
%                                               if scalar is inputted, assume both are the same, if DesiredSize(1) = 0, use sizes written in file header.
% -         interpol_method             ...     'ZeroFilling' for zerofilling/Truncation in kspace, 
%                                               any other string for bicubic interpolation in image space
% -         flip                        ...     If the FoV is rotated by 180° (e.g. for correcting for Gradient Delays), flip=1 to unflip image
% -         fredir_shift                ...     If you shifted the FoV in frequency encoding direction, the phase of the image changes dependent
%                                               on this shift in cm divided by the voxel size in frequency encoding direction. This is the case
%                                               at least for GRE images. Input this round(shift/voxel size) - ratio to correct for that phase change
%                                               Set fredir_shift = 0 for no correction
% -         Hamming_flag                ...     If 1, apply Hamming filter in k-space
% -         EllipticalFilterSize        ...     If >0, cut out an circle in k-space with radius EllipticalFilterSize
% -         phase_encod_dir             ...     If the phase encoding direction is in right-left direction, the image is rotated by 90°.
% -         sum_averages_flag           ...     If = 1, the averages will be summed                                      If you want to undo this rotation, set phase_encod_dir = 'RL', otherwise set it to anything else.
%
% Output:
% -         image                         ...     Output data in image domain. size: channel x ROW x COL x SLC
% -         image_kspace                  ...     Output data in k-space.      size: channel x ROW * OverSampling x COL x SLC
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: read_image_dat_2_1   --> memused_linux_1_0, read_ascconv_1_2, Analyze_image_mdh_1_2, EllipticalFilter_1_0, HammingFilter_1_3
%                  read_image_dicom_1_0 -->         

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.





%% 0. Assign standard values


% Split Input file in magnitude and phase file
Image_MagPha_file = regexp(image_path,' ','split');

% Read ascconv header
ParList_asc = read_ascconv_1_3(Image_MagPha_file{1});

% Assign standard values to variables if nothing is passed to function.
if(~exist('DesiredSize','var') || DesiredSize(1) == 0)
    DesiredSize(1) = ParList_asc.ROW_raw;           % There is some problem with the read_ascconv function. For image files, the FinalMatrixSize
    DesiredSize(2) = ParList_asc.COL_raw;           % Is written where for CSI files the raw-sizes are written...
end
if(numel(DesiredSize) == 1)
    DesiredSize(2) = DesiredSize;                   % If only one size is inputted, assume that both sizes should be the same.
end
clear ParList_asc

% Assign standard values to variables if nothing is passed to function.
if(~exist('interpol_method','var'))
    interpol_method = 'ZeroFilling';
end
if(~exist('flip','var'))
    flip = 0;
end
if(~exist('fredir_shift','var'))
    fredir_shift = 0;
end
if(~exist('Hamming_flag','var'))
    Hamming_flag = 0;
end
if(~exist('EllipticalFilterSize','var'))
    EllipticalFilterSize = 0;
end
if(~exist('phase_encod_dir','var'))
    phase_encod_dir = 'AP';
end
if(~exist('sum_averages_flag','var'))
    sum_averages_flag = 1;
end    
    
    
    
%% 1. Decide whether input file is DICOM or DAT; Call then appropriate function
    
if(numel(strfind(image_path, '.dat')) > 0)
    
    if(nargout > 1)
        [image,image_kspace] = read_image_dat_2_4(image_path,DesiredSize,interpol_method,flip,fredir_shift,Hamming_flag,EllipticalFilterSize, phase_encod_dir,sum_averages_flag);
    else
        image = read_image_dat_2_4(image_path,DesiredSize,interpol_method,flip,fredir_shift,Hamming_flag,EllipticalFilterSize, phase_encod_dir,sum_averages_flag);        
    end
    
else
    
    if(nargout > 1)
        [image, image_kspace] = read_image_dicom_1_2(image_path,DesiredSize,interpol_method,flip,fredir_shift,Hamming_flag,EllipticalFilterSize, phase_encod_dir);
    else
        image = read_image_dicom_1_2(image_path,DesiredSize,interpol_method,flip,fredir_shift,Hamming_flag,EllipticalFilterSize, phase_encod_dir);        
    end
        
        
end



