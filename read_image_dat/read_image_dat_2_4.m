function [image,image_kspace] = read_image_dat_2_4(image_path, DesiredSize,interpol_method,flip,fredir_shift,Hamming_flag,EllipticalFilterSize, phase_encod_dir,sum_averages_flag)
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
% -         sum_averages_flag           ...     If = 1, the averages will be summed
%
% Output:
% -         image                         ...     Output data in image domain. size: channel x ROW x COL x SLC x Averages
% -         image_kspace                  ...     Output data in k-space.      size: channel x ROW * OverSampling x COL x SLC x Averages
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux_1_0, read_ascconv_1_3, Analyze_image_mdh_1_2, EllipticalFilter_1_0, HammingFilter_1_3

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.




%% 0. Definitions, Conventions, Preparations


% fredir = frequency encoding direction; phadir = phase encoding direction
tic

% Find out memory used by MATLAB
memused_before = memused_linux_1_0(1); 




% Read ascconv header
ParList_asc = read_ascconv_1_3(image_path);
%RemoveOversampling = ParList_asc.RemoveOversampling;
AsymmetricEcho = ParList_asc.AsymmetricEcho;                                     
BaseResolution = ParList_asc.ROW_raw;
SLC = ParList_asc.SLC;

% Assign standard values to variables if nothing is passed to function.
if(~exist('DesiredSize','var') || DesiredSize(1) == 0)
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
if(~exist('sum_averages_flag','var'))
    sum_averages_flag = 1;
end
%clear ParList_asc



% READ MDH INFORMATION
ParList = Analyze_image_mdh_1_3(image_path);
phadir_measured = ParList.phadir_measured;
if(AsymmetricEcho)
    fredir_measured = 2^nextpow2(ParList.fredir_measured);  % When Gradient echo image with asymmetric echo is acquired, the k-space has to be either 
else                                                        % zerofilled or cutted, because otherwise the k-space center is not in the center
    fredir_measured = ParList.fredir_measured;
end
fredir_ReallyMeasured = ParList.fredir_measured;
Averages = ParList.Averages;
total_k_points = ParList.total_k_points;
total_channel_no = ParList.total_channel_no;
%k_center_fredir = ParList.k_center_fredir;
clear ParList


% Definitions
OversamplingFactor = fredir_measured / BaseResolution;
%OversamplingFactor_ReallyMeasured = fredir_ReallyMeasured / BaseResolution;


%% 1. READ DATA


raw_image_fid = fopen(sprintf('%s', image_path),'r');
headersize = fread(raw_image_fid,1, 'int32');
fseek(raw_image_fid, headersize,'bof');

image_kspace = zeros(total_channel_no,fredir_ReallyMeasured,phadir_measured,SLC,Averages);

for ADC_MeasNo = 1:Averages*total_k_points*total_channel_no
    
    chak_header = fread(raw_image_fid, 64, 'int16');
    k_slice = chak_header(19) + 1;                                   % What about 3d-data sets?                      
    k_line = chak_header(17) + 1;
    channel_no = chak_header(63) + 1;
    Avg = chak_header(18) + 1;                                       % Averages

    chak_data = fread(raw_image_fid, fredir_ReallyMeasured*2, 'float32');
    image_real = chak_data(1:2:end);
    image_imag = chak_data(2:2:end);
    image_complex = complex(image_real,image_imag); 
    image_kspace(channel_no,:,k_line,k_slice,Avg) = image_complex;   % AsymmetricEcho...: If due to asymmetric echo
                                                                                                % the array had to be zerofilled in fredir
end

if(sum_averages_flag)
    image_kspace = sum(image_kspace,5)/Averages;
    Averages = 1;
end

clear chak_header k_slice k_line channel_no Avg chak_data image_real image_imag image_complex
fclose(raw_image_fid);




%% 2. Perform zerofilling in case of Asymmetric Echo
% If asymmetric echo is allowed, the fredir size is probably 212 and the real k-space center is not in the center of the matrix --> zerofill array
% To correct for that

AsymmetricEchoZeros = zeros([total_channel_no fredir_measured-fredir_ReallyMeasured size(image_kspace,3) SLC Averages]);
image_kspace = cat(2,AsymmetricEchoZeros,image_kspace);
clear AsymmetricEchoZeros





%% 3. Unflip flipped image in k-space

if(flip == 1)
   image_kspace = flipdim(flipdim(image_kspace,2),3);
   image_kspace = circshift(image_kspace, [0 1 1 0 0]);
end





%% 4. Correct phase due to shift in freq encoding direction
% Strange phenomenon for images: if you shift them in frequency encoding direction, the phase of the image changes. This undoes the phase change.


if(ne(fredir_shift,0))                                                    
    image_kspace = image_kspace * exp(1i*deg2rad(fredir_shift*(360/2^nextpow2(fredir_measured))));     % Is it really 2^nextpow2(fredir_measured) and not just fredir_measured?   
end                                                                                                      





%% 5. Zerofill/Truncate kSpace to size [DesiredSize(1)*oversampling, DesiredSize(2)]

if(strcmpi(interpol_method,'ZeroFilling'))

    
    % Perform Zerofilling / kSpace-Truncation in Phase Encoding Direction
    
    % Zerofilling
    if(DesiredSize(2) > phadir_measured)
        
        % Define zero-array that should be concatenated to left and right. Overall, DesiredSize(2) - phadir_measured zeros are added
        zerofill_zeros = zeros([total_channel_no size(image_kspace,2) ceil((DesiredSize(2) - phadir_measured)/2) SLC Averages]);        
        image_kspace = cat(3,zerofill_zeros,image_kspace,zerofill_zeros);
        clear zerofill_zeros
        
    % Truncating 
    elseif(DesiredSize(2) < phadir_measured)
        
        % Define kSpace Center
        kSpaceCenter_phadir_mdh = floor(phadir_measured/2) + 1;
        % save part of measured k-space in whole desired k-space
        image_kspace = image_kspace(:,:,kSpaceCenter_phadir_mdh - floor(DesiredSize(2)/2) : kSpaceCenter_phadir_mdh + ceil(DesiredSize(2)/2) - 1,:,:);
        clear kSpaceCenter_phadir_mdh
        
    end 
    
    
    
    % Perform Zerofilling / kSpace-Truncation in Frequency Encoding Direction
    
    % Zerofilling
    if(DesiredSize(1) * OversamplingFactor > fredir_measured)
        
        % Define zero-array that should be concatenated to up and down of array. Overall, DesiredSize(1)*OversamplingFactor - fredir_measured zeros are added
        zerofill_zeros = zeros([total_channel_no ceil((DesiredSize(1)*OversamplingFactor-fredir_measured)/2) size(image_kspace,3) SLC Averages]);
        image_kspace = cat(2,zerofill_zeros,image_kspace,zerofill_zeros);
        clear zerofill_zeros
        
    % Truncating
    elseif(DesiredSize(1) * OversamplingFactor < fredir_measured)
        
        % Define kSpace Center
        kSpaceCenter_fredir_mdh = floor(fredir_measured/2) + 1;
        % save part of measured k-space in whole desired k-space
        image_kspace = image_kspace(:,kSpaceCenter_fredir_mdh - floor(DesiredSize(1)*OversamplingFactor/2) : kSpaceCenter_fredir_mdh + ceil(DesiredSize(1)*OversamplingFactor/2) - 1,:,:,:);
        clear kSpaceCenter_fredir_mdh
        
    end  
    
    
end




%% 6. Reorder Multislice Data if Necessary

if(SLC > 1 && ParList_asc.InterleavedAcquisition)
    image_kspace = reorder_multislice_image_1_0(image_kspace,4);
end



%% 7. Apply Elliptical Filter

if(EllipticalFilterSize > 0)
    image_kspace = EllipticalFilter_1_1(image_kspace,[2 3],[OversamplingFactor 1 1 EllipticalFilterSize],1);               % Apply Elliptical Filter in directions 2,3 
end




%% 8. Apply Hamming Filter

if(Hamming_flag)
    image_kspace = HammingFilter_1_3(image_kspace,[2 3],1);
end
    




%% 9. FFT FROM K-SPACE TO DIRECT SPACE

image = ifftshift(ifftshift(image_kspace,2),3);
image = fft(fft(image,[],2),[],3);
image = fftshift(fftshift(image,2),3);

if(nargout < 2)
    clear image_kspace
end








%% 10. Remove oversampling in image domain

if(OversamplingFactor ~= 1)
    
    image_center = floor(size(image,2) / 2) + 1;
    left_border = image_center - size(image,2) / (OversamplingFactor*2);                      % = image center of oversampled data - half the size after Oversampling-removal
    right_border = image_center + size(image,2) / (OversamplingFactor*2) - 1;                 % = same but + half of the size after Oversampling-removal - 1 (-1 because of image center)
    image = reshape(image(:,left_border:right_border,:,:,:), [total_channel_no, right_border-left_border+1, size(image,3), SLC, Averages]);  % truncate left and right borders
    clear image_center left_border right_border
    
end





%% 11. Interpolate data in image domain

if(~strcmpi(interpol_method,'ZeroFilling'))
    
    % Perform interpolation slice by slicee and channel by channel. This is ok for multi-slice images, but not ok for real 3d-data sets.
    image_resized = zeros(total_channel_no,DesiredSize(1),DesiredSize(2),SLC);
    for slice_index = 1:SLC
        for channel = 1:total_channel_no
            for Avg = 1:Averages
                image_resized(channel,:,:,slice_index,Avg) = imresize(squeeze(image(channel,:,:,slice_index,Avg)),[DesiredSize(1), DesiredSize(2)],'bicubic');
            end
        end
    end
    image = image_resized;
    clear image_resized
    
end




%% 12. PHASE ENCODING RL CORRECTION, FLIP LEFT & RIGHT (Physicians...)

if(strcmp(phase_encod_dir,'RL')) % If the phase encoding is set to RL, the image is rotated --> rotate back
   
    for slice_index = 1:SLC
        for channel_no = 1:total_channel_no
            for Avg = 1:Averages
                image(channel_no,:,:,slice_index,Avg) = flipdim(rot90(squeeze(image(channel_no,:,:,slice_index,Avg)),-1),2);
            end
        end
    end
    
else
    image = flipdim(image,2);
    
end







%% 13. Postparations

memused_after = memused_linux_1_0(1); 
display([char(10) 'The function used ' num2str(memused_after-memused_before) '% of the total memory.'])

fprintf('\nExecution of the function took %10.6f seconds.\n',toc) 


