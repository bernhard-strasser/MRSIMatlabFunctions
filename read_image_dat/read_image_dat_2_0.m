function [image,image_kspace] = read_image_dat_2_0(image_file, fredir_desired, phadir_desired,interpol_method,fredir_shift,phase_encod_dir,flip)
%
% read_csi_dat_x_x Read in csi-data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function only decides if input file is .dat file or DICOM according to its ending. It then calls the read_csi_dat_x_x or read_csi_dicom_x_x
% functions. Refer to these for more info
%
%
% [csi,csi_kspace] = read_csi_dat_1_10(csi_path, zerofill_to_nextpow2_flag, zerofilling_fact, x_shift,y_shift)
%
% Input: 
% -         csi_path                    ...     Path of MRS(I) file.
% -         zerofill_to_nextpow2_flag   ...     Flag, if the MRSI data should be zerofilled to the next power of 2 in k-space (e.g. 42x42 sampled --> zf to 64x64?)
% -         zerofilling_fact            ...     Factor with which the MRSI data should be zerofilled in k-space for interpolation (e.g. zerofill from 64x64 to 128x128)
% -         Hadamard_flag               ...     If data is multislice hadamard encoded, perform hadamard-decoding function
% -         x_shift                     ...     Shift the MRSI data in the left-right direction ( = row direction of matrix) by x_shift voxels
% -         y_shift                     ...     Shift the MRSI data in anterior-posterior direction ( = column direction of matrix) by y_shift voxels
%
% Output:
% -         csi                         ...     Output data in image domain. In case of Single Voxel Spectroscopy, this is the only output
% -         csi_kspace                  ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x vecSize
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: read_csi_dat_1_10, Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m, read_csi_dicom_1_1






%% 0. Definitions, Conventions, Preparations



% fredir = frequency encoding direction; phadir = phase encoding direction


% Read ascconv header
ParList_asc = read_ascconv_1_2(image_file);
RemoveOversampling = ParList_asc.RemoveOversampling;
BaseResolution = ParList_asc.ROW_raw;

% Assign standard values to variables if nothing is passed to function.
if(~exist('fredir_desired','var'))
    fredir_desired = ParList_asc.ROW_raw;
end 
if(~exist('phadir_desired','var'))
    phadir_desired = ParList_asc.COL_raw;
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
if(~exist('phase_encod_dir','var'))
    phase_encod_dir = 'AP';
end
clear ParList_asc



% READ MDH INFORMATION
ParList = Analyze_image_mdh_1_2(image_file);
phadir_measured = ParList.phadir_measured;
fredir_measured = ParList.fredir_measured;
Averages = ParList.Averages;
total_k_points = ParList.total_k_points;
total_channel_no = ParList.total_channel_no;
%k_center_fredir = ParList.k_center_fredir;
clear ParList

OversamplingFactor = fredir_measured / BaseResolution;



% %% 1. Define k-Space centers
% 
% if(strcmpi(interpol_method,'ZeroFilling'))
%     % Define k-space center of measured k-space
%     kSpaceCenter_fredir_mdh = floor(fredir_measured/2) + 1;
%     kSpaceCenter_phadir_mdh = floor(phadir_measured/2) + 1;
% 
%     % Define k-space center of desired k-space
%     kSpaceCenter_fredir_desired = floor(fredir_measured/BaseResolution * fredir_desired / 2) + 1;       % There might be oversampling; We have to keep this oversampling, and remove it in image domain
%     kSpaceCenter_phadir_desired = floor(phadir_desired/2) + 1;
% 
%     % Define shift between measured and desired k-space center
%     kSpaceCenter_fredir_shift = kSpaceCenter_fredir_mdh - kSpaceCenter_fredir_desired;
%     kSpaceCenter_phadir_shift = kSpaceCenter_phadir_mdh - kSpaceCenter_phadir_desired;
% else
%     % Define shift between measured and desired k-space center
%     kSpaceCenter_fredir_shift = 0;
%     kSpaceCenter_phadir_shift = 0;
% end






%% 1. READ ALL DATA


raw_image_fid = fopen(sprintf('%s', image_file),'r');
headersize = fread(raw_image_fid,1, 'int32');
fseek(raw_image_fid, headersize,'bof');

image_kspace = zeros(total_channel_no,fredir_measured,phadir_measured,SLC,Averages);


for ADC_MeasNo = 1:Averages*total_k_points*total_channel_no
    
    chak_header = fread(raw_image_fid, 64, 'int16');
    k_slice = chak_header(19) + 1;                                   % What about 3d-data sets?                      
    k_line = chak_header(17) + 1;
    channel_no = chak_header(63) + 1;
    Avg = chak_header(24) + 1;                                       % Averages

    
    chak_data = fread(raw_image_fid, fredir_measured*2, 'float32');
    image_real = chak_data(1:2:end);
    image_imag = chak_data(2:2:end);
    image_complex = complex(image_real,image_imag); 
    image_kspace(channel_no,:,k_line,k_slice,Avg) = image_complex; 
    

end


fclose(raw_image_fid);




%% 2. Unflip flipped image in k-space

if(flip == 1)
   image_kspace = flipdim(flipdim(image_kspace,2),3);
   image_kspace = circshift(image_kspace, [0 1 1 0]);
end





%% 3. Correct phase due to shift in freq encoding direction
% Strange phenomenon for images: if you shift them in frequency encoding direction, the phase of the image changes. This undoes the phase change.


if(ne(fredir_shift,0))                                                    
    image_kspace = image_kspace * exp(1i*deg2rad(fredir_shift*(360/2^nextpow2(fredir_measured))));     % Is it really 2^nextpow2(fredir_measured) and not just fredir_measured?   
end                                                                                                      





%% 3. Zerofill data to size [fredir_desired*oversampling, phadir_desired]

if(strcmpi(interpol_method,'ZeroFilling') && (fredir_desired*OversamplingFactor ~= fredir_measured || phadir_desired ~= phadir_measured))
    
    % 3.1 Define kSpace Centers
    % Define k-space center of measured k-space
    kSpaceCenter_fredir_mdh = floor(fredir_measured/2) + 1;
    kSpaceCenter_phadir_mdh = floor(phadir_measured/2) + 1;

    % Define k-space center of desired k-space
    kSpaceCenter_fredir_desired = floor(OversamplingFactor * fredir_desired / 2) + 1;       % There might be oversampling; We have to keep this oversampling, and remove it in image domain later
    kSpaceCenter_phadir_desired = floor(phadir_desired/2) + 1;

    
    
    
    % 3.2 Perform Zerofilling / kSpace-Truncation in Both Directions
    
    % fredir
    if(fredir_desired * OversamplingFactor > fredir_measured)
        image_kspace_zf = zeros(total_channel_no,fredir_desired * OversamplingFactor,phadir_measured,SLC,Averages);
        % save whole measured k-space in part of desired k-space
        image_kspace_zf(:,kSpaceCenter_fredir_desired - fredir_measured/2 : kSpaceCenter_fredir_desired + fredir_measured/2 - 1,:,:,:) = image_kspace;
        image_kspace = image_kspace_zf;
        clear image_kspace_zf
    elseif(fredir_desired * OversamplingFactor < fredir_measured)
        % save part of measured k-space in whole desired k-space
        image_kspace = image_kspace(:,kSpaceCenter_fredir_mdh - fredir_desired*OversamplingFactor/2 : kSpaceCenter_fredir_mdh + fredir_desired*OversamplingFactor/2 - 1,:,:,:);
    else
        image_kspace = image_kspace;
    end
    
    % phadir
    if(phadir_desired > phadir_measured)
        image_kspace_zf = zeros(total_channel_no,fredir_desired * OversamplingFactor,phadir_desired,SLC,Averages);        
        image_kspace(:,:,kSpaceCenter_phadir_desired - phadir_measured/2 : kSpaceCenter_phadir_desired + phadir_measured/2 - 1,:,:) = image_kspace_zf;        
    elseif(phadir_desired < phadir_measured)
        image_kspace = image_kspace_zf(:,:,kSpaceCenter_phadir_mdh - phadir_desired/2 : kSpaceCenter_phadir_mdh + phadir_desired/2 - 1,:,:);
    else
        image_kspace = image_kspace_zf;
    end   
    
    
end





%% 5. Apply Elliptical Filter






%% 6. Apply Hamming Filter






%% 7. FFT FROM K-SPACE TO DIRECT SPACE

image = ifftshift(ifftshift(image_kspace,2),3);
image = fft(fft(image,[],2),[],3);
image = fftshift(fftshift(image,2),3);





%% 8. REMOVE OVERSAMPLING IN IMAGE DOMAIN

% TRUNCATE IN XSPACE FROM [fredir_desired*oversampling_fredir, phadir_desired] TO [fredir_desired, phadir_desired]
image = reshape(image(:,fredir_left_border_xspace:fredir_right_border_xspace,:,:), [total_channel_no phadir_measured phadir_measured, SLC]);    %truncate left and right borders




%% 9. INTERPOLATE DATA





% TRUNCATE IN XSPACE FROM [fredir_desired*oversampling_fredir, phadir_desired] TO [fredir_desired, phadir_desired]
image = reshape(image(:,fredir_left_border_xspace:fredir_right_border_xspace,:,:), [total_channel_no phadir_measured phadir_measured, SLC]);    %truncate left and right borders


% interpolate data

image_resized = zeros(total_channel_no,fredir_desired,phadir_desired,SLC);
for slice_index = 1:SLC
    for channel = 1:total_channel_no
        image_resized(channel,:,:,slice_index) = imresize(squeeze(image(channel,:,:,slice_index)),[fredir_desired, phadir_desired],'bicubic');
    end
end
image = image_resized;





%% 10. PHASE ENCODING RL CORRECTION, FLIP LEFT & RIGHT (Physicians...)
% If the phase encoding is set to RL, the image is rotated --> rotate back

if(strcmp(phase_encod_dir,'RL'))
   
    for channel_no = 1:total_channel_no
        image(channel_no,:,:) = flipdim(rot90(squeeze(image(channel_no,:,:)),-1),2);
    end
    
else
    image = flipdim(image,2);
    
end







