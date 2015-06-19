function mask_grown = MaskGrow(mask,NoOfVoxels,SliceBySlice_flag)
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
% -         NoiseCorrMat                ...     If size(NoiseCorrMat) = [cha cha]: the k-space Data gets decorrelated with this matrix. 
%                                               If NoiseCorrMat = 0, or not existant: No Noise Decorrelation performed
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

if(~exist('mask','var'))
    mask_grown = 0;
	fprintf('Input needed.')
	return
end

   
if(~exist('NoOfVoxels','var') || NoOfVoxels < 1)
    NoOfVoxels = 1;
end
if(~exist('SliceBySlice_flag','var'))
    SliceBySlice_flag = true;
end



%% Grow

% Very slow implementation but works
mask_grown = mask;

if(SliceBySlice_flag)
	mask2 = EllipticalFilter(ones(size(mask(:,:,1))), [1 2], [1 1 1 NoOfVoxels],1);
else
	mask2 = EllipticalFilter(ones(size(mask)), [1 2 3], [1 1 1 NoOfVoxels],1);
end

for z = 1:size(mask,3)
	if(SliceBySlice_flag)
		mask4 = mask(:,:,z);
	end
	for x = 1:size(mask,1)
		for y = 1:size(mask,2)
			if(mask(x,y,z) == 1)
				continue;
			end			

			% create a second mask with radius NoOfVoxels
			ShiftBy = [x-(ceil(size(mask2,1)/2))-1, y-(ceil(size(mask2,1)/2))-1, z-(ceil(size(mask2,1)/2))-1];
			ShiftBy = ShiftBy(1:numel(size(mask2)));
			mask3 = NonCircShift(mask2,ShiftBy);
			if(~isempty(find(mask4(logical(mask3)),1)))
				mask_grown(x,y,z) = 1;
			end

		end
	end
end

	



   