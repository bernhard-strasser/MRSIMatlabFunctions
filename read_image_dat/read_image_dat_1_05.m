function [image_mat_complex_dat] = read_image_dat_1_1(image_dat_file,total_channel_no, ROW,COL, freq_encod_os, phase_encoding_dir, x_shift,y_shift)
%% 0. Preparations
SLC = 1;

oversampling = (2^nextpow2(ROW_oversampling)-ROW)/ROW + 1;

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

image_mat_complex_dat_kspace = zeros(total_channel_no,2^nextpow2(ROW_oversampling),COL,1);

for k_line = 1:COL
    for channel_no = 1:total_channel_no

        chak_header = fread(raw_image_fid, 64, 'int16');
        chak_data = fread(raw_image_fid, ROW_oversampling*2, 'float32');
        image_real = chak_data(1:2:end);
        
        image_imag = chak_data(2:2:end);
        image_complex = complex(image_real,image_imag); 
        image_mat_complex_dat_kspace(channel_no,(2^nextpow2(ROW_oversampling))/2-(k_x_center-1)+1:end,k_line) = image_complex; 
    end
end

fclose(raw_image_fid);


%% 2. SHIFT IN K-SPACE, FFT, FLIP LEFT AND RIGHT

%shift to the left and then down

% total_channel_no = 4;
% xxx = real(repmat(reshape(exp(1i*y_shift*2*pi/COL*(1:COL)), [1,1,COL,SLC]),[total_channel_no,oversampling*ROW,1,SLC]));
% 
% size(reshape(exp(1i*y_shift*2*pi/COL*(1:COL)), [1,1,COL,SLC]))
% size(repmat(reshape(exp(1i*y_shift*2*pi/COL*(1:COL)), [1,1,COL,SLC]),[total_channel_no,oversampling*ROW,1,SLC]))
% size(reshape(repmat(reshape(exp(1i*y_shift*2*pi/COL*(1:COL)), [1,1,COL,SLC]),[total_channel_no,oversampling*ROW,1,SLC]), [total_channel_no,oversampling*ROW,COL,SLC]))
% squeeze(xxx(1:4,1,50:54,1))
% squeeze(xxx(1,4:8,30:34,1))

if(exist('x_shift','var'))
    image_mat_complex_dat_kspace = repmat(exp(1i*x_shift*2*pi/ROW*(1:oversampling*ROW)),[total_channel_no,1,COL,SLC]) .* image_mat_complex_dat_kspace;         %shift 2 256 voxels instead of 1 128 voxel; oversampling 
    %image_mat_complex_dat_kspace = repmat(exp(-1i*x_shift*2/2*pi/ROW*(1:oversampling*ROW)),[total_channel_no,1,COL,SLC]) .* image_mat_complex_dat_kspace;         %shift 2 256 voxels instead of 1 128 voxel; oversampling 
    %image_mat_complex_dat_kspace = image_mat_complex_dat_kspace * exp(-1i*deg2rad(x_shift*oversampling*(180 + 360/(ROW*oversampling))));                       %cancels in x_shift*os/(ROW*os)
    image_mat_complex_dat_kspace = image_mat_complex_dat_kspace * exp(-1i*deg2rad(x_shift*1.5*oversampling*(180 + 360/(ROW*oversampling))));                       %cancels in x_shift*os/(ROW*os)
end

if(exist('y_shift','var'))   
    image_mat_complex_dat_kspace = repmat(reshape(exp(1i*y_shift*2*pi/COL*(1:COL)), [1,1,COL,SLC]),[total_channel_no,oversampling*ROW,1,SLC]) .* image_mat_complex_dat_kspace;
    image_mat_complex_dat_kspace = image_mat_complex_dat_kspace * exp(-1i*deg2rad(y_shift*(180 + 360/COL))); %shift results in offset in phase, this corrects that
end

image_mat_complex_dat = ifftshift(ifftshift(image_mat_complex_dat_kspace,2),3);
image_mat_complex_dat = fft(fft(image_mat_complex_dat,2^nextpow2(ROW_oversampling),2),COL,3);
image_mat_complex_dat = fftshift(fftshift(image_mat_complex_dat,2),3);
image_mat_complex_dat = flipdim(image_mat_complex_dat,2);
truncate_ROW_voxels = (2^nextpow2(ROW_oversampling) - ROW)/2;
left_border = truncate_ROW_voxels + 1;
right_border = 2^nextpow2(ROW_oversampling) - truncate_ROW_voxels;
image_mat_complex_dat = reshape(image_mat_complex_dat(1,left_border:right_border,:,1), [total_channel_no ROW COL 1]);    %truncate left and right borders
%image_mat_complex_dat = circshift(image_mat_complex_dat, [0 1 0 0]);