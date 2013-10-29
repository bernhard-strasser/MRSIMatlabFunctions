function [image_dat,image_dat_kspace] = read_image_dat_1_1(image_dat_file,total_channel_no, ROW, ROW_oversampling,COL,x_shift,y_shift,ROW_zf_factor)
%% 0. Preparations
SLC = 1;

ROW_os_nextpow2 = 2^nextpow2(ROW_oversampling);     % the data is sometimes oversampled from 128 to 212 and sometimes to 256 etc. So in general the data gets zerofilled if that is the case
%ROW_zf_factor = 32;                                 % defines the zerofilling factor, with which the hybridspace (in x_direction x_space in y_direction kspace),
                                                    % gets oversampled in order to find out the point in kspace with the highest magnitude
pause on                                            % enables the usage of pause

%% 1. READ DATA

% %read in ROW, COL, ROW_oversampling DOESNT WORK
raw_image_fid = fopen(sprintf('%s', image_dat_file),'r');
% fseek(raw_image_fid, 3600,'bof');
% chak_header = fread(raw_image_fid, 100, 'int8=>char');
% chak_header = chak_header';
% ROW = str2double(chak_header(41:43))
% COL = chak_header(77:79)
% fseek(raw_image_fid, 476832,'bof');
% chak_header = fread(raw_image_fid, 64, 'int16');
% ROW_oversampling = chak_header(15)
% k_x_center = chak_header(33) + 1

%read data
fseek(raw_image_fid, -256*total_channel_no-total_channel_no*COL*(128+ROW_oversampling*2*4),'eof');
chak_header = fread(raw_image_fid, 64, 'int16');
k_x_center = chak_header(33)+1;
fseek(raw_image_fid, -256*total_channel_no-total_channel_no*COL*(128+ROW_oversampling*2*4),'eof');

image_dat_kspace = zeros(total_channel_no,ROW_os_nextpow2,COL,1);

for k_line = 1:COL
    for channel_no = 1:total_channel_no

        chak_header = fread(raw_image_fid, 64, 'int16');
        chak_data = fread(raw_image_fid, ROW_oversampling*2, 'float32');
        image_real = chak_data(1:2:end);
        
        image_imag = chak_data(2:2:end);
        image_complex = complex(image_real,image_imag); 
        image_dat_kspace(channel_no,ROW_os_nextpow2/2-(k_x_center-1)+1:end,k_line,1) = image_complex; 
    end
end

fclose(raw_image_fid);


% % DEBUG MODE
% % figure
% % imagesc(abs(squeeze(image_dat_kspace)))
% % pause
% display 'after read in'
% size(image_dat_kspace)
% % DEBUG MODE END


%% CORRECTION: 2.1. DO X_SHIFT IN K-SPACE

if(exist('x_shift','var') && ne(x_shift,0))
    image_dat_kspace = repmat(exp(1i*x_shift*2*pi/ROW_os_nextpow2*(1:ROW_os_nextpow2)),[total_channel_no,1,COL,SLC]) .* image_dat_kspace;   % shift 1 256 voxels instead of 1 128 voxel; The resolut.
    image_dat_kspace = image_dat_kspace * exp(-1i*deg2rad(x_shift*180));                                                                    % does not increase, there is just crap on the borders
end                                                                                                                                         % so you still have to shift 1 voxel.


if(exist('y_shift','var') && ne(y_shift,0))   
    image_dat_kspace = repmat(reshape(exp(1i*y_shift*2*pi/COL*(1:COL)), [1,1,COL,SLC]),[total_channel_no,ROW_os_nextpow2,1,SLC]) .* image_dat_kspace;
    image_dat_kspace = image_dat_kspace * exp(-1i*deg2rad(y_shift*(180 + 360/COL)));                                                        %shift results in offset in phase, this corrects that
end



%% CORRECTION 2.2: ZERO PADDING IN HYBRIDSPACE (1. spatial dim in x_space, 2. in k) AT BEGINNING AND END OF HYBRID SPACE AND SHIFT KSPACE

if(ROW_zf_factor > 1)
    
    % TRANSFORM TO HYBRIDSPACE: FOURIER TRANSFORM 1. SPATIAL DIM
    image_dat_hybridspace = fftshift(fft(ifftshift(image_dat_kspace,2),ROW_os_nextpow2,2),2);

    % ZERO PADDING SO THAT RESULT IS ROW_zf_factor LARGER THAN ORIGINAL 
    image_dat_hybridspace_zf = cat(2,zeros(total_channel_no,ROW_os_nextpow2*(ROW_zf_factor-1)/2,COL), image_dat_hybridspace, zeros(total_channel_no,ROW_os_nextpow2*(ROW_zf_factor-1)/2,COL));

    % GO FROM HYBRIDSPACE BACK TO K_SPACE BY IFFT
    image_dat_kspace_zf = fftshift(ifft(ifftshift(image_dat_hybridspace_zf,2),ROW_os_nextpow2*ROW_zf_factor,2),2);


    % GO FROM HYBRIDSPACE FURTHER TO X_SPACE BY FFT OF 2. SPATIAL DIMENSION
    % image_dat_xspace_zf = fftshift(fft(ifftshift(image_dat_hybridspace_zf,3),COL,3),3);


%     % DEBUG MODE
%     display 'before zero padding'
%     size(image_dat_hybridspace)
%     display 'after zero padding'
%     size(image_dat_hybridspace_zf)
%     display 'in kspace'
%     size(image_dat_kspace_zf)
% 
%     figure
%     imagesc(abs(squeeze(image_dat_kspace_zf(1,:,:))))
%     figure
%     imagesc(rot90(abs(squeeze(image_dat_hybridspace_zf(1,:,:))),-1))
%     pause
%     % DEBUG MODE END
    
%% CORRECTION 2.3: SHIFT KSPACE
    %shift kspace because center of kspace is not point with highest abs value otherwise

    for channel_no = 1:total_channel_no

        abs_kspace = abs(squeeze(image_dat_kspace_zf(channel_no,:,:)));

        % FIND POSITION OF HIGHEST KSPACE MAGNITUDE
        max_kspace_value = max(max(abs_kspace));                                    % compute the maximum
        [max_ROW_index, max_COL_index] = find(abs_kspace == max_kspace_value);      % find the position

        kspace_center = [ROW_os_nextpow2*ROW_zf_factor/2+1, COL/2+1];               % compute the position of the k_space_center
        kspace_shift = [0 (kspace_center - [max_ROW_index, max_COL_index])];        % difference of the k_space_centre to the position of the highest k_space_value: this is the value with which the kspace should be shifted
        % CONSIDER THAT THE KSPACE IS 256x128, K-SPACE CENTER = [129,65]; THE GREATEST MAGNITUDE VALUE IN KSPACE IS AT [127,64]; THE DISTANCE IS [129-127,65-64]=[2,1]. SO WE HAVE TO CIRCSHIFT THE MATRIX 2 ROWS DOWN AND 1 COL 
        % TO THE LEFT; THAT MEANS: CIRCSHIFT(MATRIX,[2,1])
        image_dat_kspace(channel_no,:,:) = circshift(image_dat_kspace_zf, kspace_shift);


    % DEBUG MODE
    %     display 'after kspace shift'
    %     size(image_dat_kspace_zf)
    %     figure
    %     imagesc(abs(squeeze(image_dat_kspace_zf(1,:,:))))
    %     pause
    % DEBUG MODE END

    end

    
    
    
end

%% 3. FFT AND TRUNCATE IMAGE AGAIN, FLIP LEFT AND RIGHT

image_dat = ifftshift(ifftshift(image_dat_kspace,2),3);
image_dat = fft(fft(image_dat,ROW_os_nextpow2*ROW_zf_factor,2),COL,3);
image_dat = fftshift(fftshift(image_dat,2),3);

truncate_ROW_voxels = ROW_os_nextpow2*ROW_zf_factor - ROW;
left_border = truncate_ROW_voxels/2 + 1;
right_border = ROW_os_nextpow2*ROW_zf_factor - truncate_ROW_voxels/2;

image_dat = reshape(image_dat(:,left_border:right_border,:), [total_channel_no ROW COL, 1]);    %truncate left and right borders

% display 'after shift in kspace'
% size(image_dat_kspace)


image_dat = flipdim(image_dat,2);


%image_dat = circshift(image_dat, [0 1 0 0]);