function [image_dicom, image_dicom_kspace] = read_image_dicom_1_1(image_path,DesiredSize,interpol_method,flip,fredir_shift,Hamming_flag,EllipticalFilterSize, phase_encod_dir)
%
% read_image_dat_x_x Read in image-data in raw format
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MR images (or at least some of them...) in Siemens raw data format.
%
%
% [image,image_kspace] = read_image_dat_2_1(image_path, DesiredSize,interpol_method,flip,fredir_shift,phase_encod_dir)
%
% Input: 
% -         image_path                  ...     Path of MRI file.
% -         DesiredSize                 ...     [desired size in freq encod dir, desired size in phase encod dir] 
%                                               if scalar is inputted, assume both are the same
% -         interpol_method             ...     'ZeroFilling' for zerofilling/Truncation in kspace, 
%                                               any other string for bicubic interpolation in image space
% -         flip                        ...     If the FoV is rotated by 180° (e.g. for correcting for Gradient Delays), flip=1 to unflip image
% -         fredir_shift                ...     If you shifted the FoV in frequency encoding direction, the phase of the image changes dependent
%                                               on this shift in cm divided by the voxel size in frequency encoding direction. This is the case
%                                               at least for GRE images. Input this shift/voxel size - ratio to correct for that phase change
%                                               Set fredir_shift = 0 for no correction
% -         Hamming_flag                ...     If 1, apply Hamming filter in k-space
% -         EllipticalFilterSize        ...     If >0, cut out an circle in k-space with radius EllipticalFilterSize
% -         phase_encod_dir             ...     If the phase encoding direction is in right-left direction, the image is rotated by 90°.
%                                               If you want to undo this rotation, set phase_encod_dir = 'RL', otherwise set it to anything else.
%
% Output:
% -         image                         ...     Output data in image domain. size: channel x ROW x COL x SLC
% -         image_kspace                  ...     Output data in k-space.      size: channel x ROW * OverSampling x COL x SLC
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: 













%% 0. Preparations

% Read ascconv header
ParList_asc = read_ascconv_1_3(image_path);
%RemoveOversampling = ParList_asc.RemoveOversampling;
BaseResolution = ParList_asc.ROW_raw;
SLC = ParList_asc.SLC;

% Assign standard values to variables if nothing is passed to function.
if(~exist('DesiredSize','var'))
    DesiredSize(1) = ParList_asc.ROW_raw;           % There is some problem with the read_ascconv function. For image files, the FinalMatrixSize
    DesiredSize(2) = ParList_asc.COL_raw;           % Is written where for CSI files the raw-sizes are written...
end
if(numel(DesiredSize) == 1)
    DesiredSize(2) = DesiredSize;                   % If only one size is inputted, assume that both sizes should be the same.
end
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
clear ParList_asc



% Split Input file in magnitude and phase file
Image_MagPha_file = regexp(image_path,'\t','split');
Image_MagPha_file = regexp(Image_MagPha_file,' ','split');
Image_MagPha_file = regexp(Image_MagPha_file,',','split');



%% 1. READ DATA


% INITIALIZE
file_index = 0;
Image_Matrix_AbsPha = zeros(2,total_channel_no,ROW,COL,SLC); % create 5D matrix (mag/phase ,channel,ROW,COL,SLC


for Image_file = Image_MagPha_file

    file_index = file_index + 1;
    if(mod(file_index,2) == 1)
        MagPha_Index = 1;
    else
        MagPha_Index = 2;
    end
    
    
    % OPEN FILE, read headersize, go to the point where the real data starts
    fid = fopen(Image_file{1},'r');
    %headersize = fread(fid,1, 'int16')
    %fseek(fid, headersize,'bof');             % go to beginning of dataheader/data for reading_data = fread(fid,'float32');
    fseek(fid, -(ROW*COL*SLC*2), 'eof');
    
    % READ IN
    k = 0;
    for z=1:SLC
        k = k + 1;
        Image_Matrix_AbsPha(MagPha_Index,1,:,:,z) = fread(fid, [ROW,COL], 'uint16');
    end
    fclose(fid);
    
    
%     % Read data ALTERNATIVE:
%     Image_data = fread(fid,'uint16');
%     fclose(fid);
%
%     k=0;
%     for z=1:SLC
%         for y=1:ROW
%             for x=1:COL
%                 k=k+1;
%                 Image_Matrix_AbsPha(MagPha_Index,1,x,y,z) = Image_data(k);                   % Number 1: Channel Index, not yet implemented for multichannel
%             end
%         end
%     end

end


image_dicom = zeros([total_channel_no,ROW,COL,SLC]);
image_dicom(:,:,:,:) = Image_Matrix_AbsPha(1,:,:,:,:) .* exp(1i*pi*(Image_Matrix_AbsPha(2,:,:,:,:)*2 - 4096)/4096);
image_dicom = reshape(circshift(squeeze(image_dicom), [-1 0]), [1, ROW, COL, SLC]);




%% 2. Resample image if necessary


image_dicom_resized = zeros(total_channel_no,fredir_desired,phadir_desired);
for channel = 1:total_channel_no
    image_dicom_resized(channel,:,:) = imresize(squeeze(image_dicom(channel,:,:)),[fredir_desired, phadir_desired],'bicubic');
end

image_dicom = image_dicom_resized;
image_dicom_kspace = 0;




%% 3. Add phase to correct for FoV-shifts in frequency encoding direction


if(exist('fredir_shift','var') && ne(fredir_shift,0))                                                    
    image_dicom = image_dicom * exp(1i*deg2rad(fredir_shift*(360/2*ROW)));   % 2*ROW: Actually fredir_oversampling * ROW should be here, but this value is not written in DICOM file ?!      
end                                                                                                    


