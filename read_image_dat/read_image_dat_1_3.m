%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%    FUNCTION TO READ IN THE .dat IMAGING FILES    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [image_dat,image_dat_kspace] = read_image_dat_1_3(image_dat_file,total_channel_no, ROW, ROW_oversampling,COL,x_shift,y_shift)
%% 0. Preparations
SLC = 1;

ROW_os_nextpow2 = 2^nextpow2(ROW_oversampling);     % the data is sometimes oversampled from 128 to 212 and sometimes to 256 etc. So in general the data gets zerofilled if that is the case
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

%image_dat_kspace = circshift(image_dat_kspace, [0 -15 0 0]);


% % DEBUG MODE
% figure
% imagesc(abs(squeeze(image_dat_kspace)))
% pause
% display 'after read in'
% size(image_dat_kspace)
% % DEBUG MODE END


%% CORRECTION 2.1: SET K_SPACE POINTS THAT ARE MUCH TOO HIGH (NOT IN CENTER!) TO VALUE OF NEIGHBOURING K-POINTS

size_k_middle = 0.03;
replace_sensitivity = 0.8;   % higher values --> less points get changed
averaging_dist_to_kcenter = 1;

change_k_points = zeros(total_channel_no, ROW_os_nextpow2, COL);

mean_of_kspace_middlepoints = mean(mean(abs(image_dat_kspace(:,ROW_os_nextpow2/2-averaging_dist_to_kcenter:ROW_os_nextpow2/2+averaging_dist_to_kcenter,COL/2-averaging_dist_to_kcenter:COL/2-averaging_dist_to_kcenter)),2),3);
change_k_points(abs(image_dat_kspace) > repmat(replace_sensitivity*mean_of_kspace_middlepoints, [1 ROW_os_nextpow2 COL])) = 1;

%change_k_points(abs(image_dat_kspace) > repmat(replace_sensitivity*abs(image_dat_kspace(:,ROW_os_nextpow2/2+1,COL/2+1)), [1 ROW_os_nextpow2 COL])) = 1;

change_k_points(:,floor(ROW_os_nextpow2/2+1-size_k_middle*ROW_os_nextpow2):ceil(ROW_os_nextpow2/2+1+size_k_middle*ROW_os_nextpow2),floor(COL/2+1-size_k_middle*COL):ceil(COL/2+1+size_k_middle*COL)) = 0;


% %DEBUG MODE: WHICH K_SPACE POINTS GET CHANGED, WHICH NOT?
% 
% dont_change = zeros(total_channel_no, ROW_os_nextpow2, COL);
% dont_change(:,floor(ROW_os_nextpow2/2+1-size_k_middle*ROW_os_nextpow2):ceil(ROW_os_nextpow2/2+1+size_k_middle*ROW_os_nextpow2),floor(COL/2+1-size_k_middle*COL):ceil(COL/2+1+size_k_middle*COL)) = 1;
% 
% image_dat_kspace_changed_points = image_dat_kspace;
% image_dat_kspace_changed_points(logical(change_k_points)) = 10;
% 
% image_dat_kspace_dont_change_points = image_dat_kspace;
% image_dat_kspace_dont_change_points(logical(dont_change)) = 10;
% 
% figure
% imagesc(abs(squeeze(image_dat_kspace_changed_points(1,:,:))), [0 4E-3])
% title('kspace and changed points')
% 
% figure
% imagesc(abs(squeeze(image_dat_kspace_dont_change_points(1,:,:))), [0 4E-3])
% title('kspace middle dont change points')
% 
% figure
% imagesc(squeeze(change_k_points(1,:,:)))
% title('kspace only changed points')
% pause
% 
% % DEBUG MODE END


% FIND OUT THE INDICES OF THE K_POINTS TO CHANGE
% linear_index = find(change_k_points == 1);
% [change_k_channel, change_k_row, change_k_col] = ind2sub(size(change_k_points),linear_index); 



image_dat_kspace(logical(change_k_points)) = 0;

% figure
% imagesc(abs(squeeze(image_dat_kspace(2,:,:))))



%% CORRECTION: 2.2. DO X_SHIFT IN K-SPACE

if(exist('x_shift','var') && ne(x_shift,0))
    image_dat_kspace = repmat(exp(1i*x_shift*2*pi/ROW_os_nextpow2*(1:ROW_os_nextpow2)),[total_channel_no,1,COL,SLC]) .* image_dat_kspace;  % Shift in x-direction according to Fourier Shift Theorem.
    image_dat_kspace = image_dat_kspace * exp(-1i*deg2rad(x_shift*(360/ROW_os_nextpow2)));                                                 % Shift results in offset in phase, this corrects that
end                                                                                                                                        % by multiplying with constant phase.


if(exist('y_shift','var') && ne(y_shift,0))   
    image_dat_kspace = repmat(reshape(exp(1i*y_shift*2*pi/COL*(1:COL)), [1,1,COL,SLC]),[total_channel_no,ROW_os_nextpow2,1,SLC]) .* image_dat_kspace;
    image_dat_kspace = image_dat_kspace * exp(-1i*deg2rad(y_shift*(360/COL)));                                                        
end


%% 3. FFT AND TRUNCATE IMAGE AGAIN, FLIP LEFT AND RIGHT


image_dat_kspace2 = image_dat_kspace(:,ROW_os_nextpow2/4+1:3*ROW_os_nextpow2/4,COL/4+1:3*COL/4,:);
image_dat = ifftshift(ifftshift(image_dat_kspace2,2),3);
image_dat = fft(fft(image_dat,ROW,2),COL/2,3);
image_dat = fftshift(fftshift(image_dat,2),3);

truncate_ROW_voxels = ROW_os_nextpow2/2 - ROW/2; % 128 - 64 = 64
left_border = truncate_ROW_voxels/2 + 1; 
right_border = ROW_os_nextpow2/2 - truncate_ROW_voxels/2; % 128 - 64/2 = 96

image_dat = reshape(image_dat(:,left_border:right_border,:), [total_channel_no ROW/2 COL/2, 1]);    %truncate left and right borders

% display 'after shift in kspace'
% size(image_dat_kspace)


image_dat = flipdim(image_dat,2);




% image_dat = ifftshift(ifftshift(image_dat_kspace,2),3);
% image_dat = fft(fft(image_dat,ROW_os_nextpow2,2),COL,3);
% image_dat = fftshift(fftshift(image_dat,2),3);
% 
% truncate_ROW_voxels = ROW_os_nextpow2 - ROW;
% left_border = truncate_ROW_voxels/2 + 1;
% right_border = ROW_os_nextpow2 - truncate_ROW_voxels/2;
% 
% image_dat = reshape(image_dat(:,left_border:right_border,:), [total_channel_no ROW COL, 1]);    %truncate left and right borders
% 
% % display 'after shift in kspace'
% % size(image_dat_kspace)
% 
% 
% image_dat = flipdim(image_dat,2);
% 
% 
% %image_dat = circshift(image_dat, [0 1 0 0]);