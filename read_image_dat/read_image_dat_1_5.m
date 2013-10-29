%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%    FUNCTION TO READ IN THE .dat IMAGING FILES    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% fredir = frequency encoding direction; phadir = phase encoding direction



function [image_dat,image_dat_kspace] = read_image_dat_1_5(image_dat_file, fredir_desired, phadir_desired,kspace_corr_flag,fredir_shift,phase_encod_dir)






%% 0. Preparations


pause on


% READ SIZE OF HEADER

[headersize,k_fredir_center,fredir_measured,phadir_measured,total_channel_no] = read_image_dat_meas_header_1_0(image_dat_file);

raw_image_fid = fopen(sprintf('%s', image_dat_file),'r');


% SLC, FREQUENCY ENCOD DIR OVERSAMPLING, OVERSAMPLING IN FREQ ENCOD DIR, PAUSE
%SLC = 1;
fredir_os_nextpow2 = 2^nextpow2(fredir_measured);   % the data is sometimes oversampled from 128 to 212 and sometimes to 256 etc. So in general the data gets zerofilled if that is the case
oversampling_fredir = (fredir_os_nextpow2-phadir_measured)/phadir_measured+1;


% if user didn't parse fredir_desired and phadir_desired

if(~exist('fredir_desired','var'))
    fredir_desired = phadir_measured;
end
    
if(~exist('phadir_desired','var'))
    phadir_desired = phadir_measured;
end  




% DEFINE CENTER OF ENLARGED KSPACE (NEXTPOW2)
k_center = [fredir_os_nextpow2/2+1 phadir_measured/2+1];


% LEFT/RIGHT BORDERS OF TRUNCATION IN KSPACE
fredir_left_border_kspace = (fredir_os_nextpow2-fredir_desired*oversampling_fredir)/2+1;
fredir_right_border_kspace = (fredir_os_nextpow2+fredir_desired*oversampling_fredir)/2;

 
% NUMBER OF TRUNCATION VOXELS IN FREQ ENCOD DIR, LEFT/RIGHT BORDERS
truncate_fredir_voxels = fredir_desired*(oversampling_fredir-1);        
fredir_left_border_xspace = truncate_fredir_voxels/2 + 1;
fredir_right_border_xspace =  fredir_desired*oversampling_fredir - truncate_fredir_voxels/2;


%% 1. READ DATA

fseek(raw_image_fid, headersize,'bof');             %go back to beginning of dataheader/data for readin
image_dat_kspace = zeros(total_channel_no,fredir_os_nextpow2,phadir_measured,1);

for k_line = 1:phadir_measured
    for channel_no = 1:total_channel_no
        fseek(raw_image_fid,128,'cof');
        %chak_header = fread(raw_image_fid, 64, 'uint16');
        chak_data = fread(raw_image_fid, fredir_measured*2, 'float32');
        image_real = chak_data(1:2:end);
        
        image_imag = chak_data(2:2:end);
        image_complex = complex(image_real,image_imag); 
        image_dat_kspace(channel_no,fredir_os_nextpow2/2-(k_fredir_center-1)+1:end,k_line,1) = image_complex; 
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




%% CORRECTION 2.1: SET K_SPACE POINTS THAT ARE MUCH TOO HIGH (NOT IN CENTER!) TO VALUE OF NEIGHBOURING K-POINTS



if(exist('kspace_corr_flag','var') && logical(kspace_corr_flag))
    size_k_middle = 0.03;
    replace_sensitivity = 0.8;   % higher values --> less points get changed
    averaging_dist_to_kcenter = 1;

    change_k_points = zeros(total_channel_no, fredir_os_nextpow2, phadir_measured);

    mean_of_kspace_middlepoints = mean(mean(abs(image_dat_kspace(:,fredir_os_nextpow2/2-averaging_dist_to_kcenter:fredir_os_nextpow2/2+averaging_dist_to_kcenter, ...
    phadir_measured/2-averaging_dist_to_kcenter:phadir_measured/2-averaging_dist_to_kcenter)),2),3);
    change_k_points(abs(image_dat_kspace) > repmat(replace_sensitivity*mean_of_kspace_middlepoints, [1 fredir_os_nextpow2 phadir_measured])) = 1;

    %change_k_points(abs(image_dat_kspace) > repmat(replace_sensitivity*abs(image_dat_kspace(:,fredir_os_nextpow2/2+1,phadir_measured/2+1)), [1 fredir_os_nextpow2 phadir_measured])) = 1;

    change_k_points(:,floor(fredir_os_nextpow2/2+1-size_k_middle*fredir_os_nextpow2):ceil(fredir_os_nextpow2/2+1+size_k_middle*fredir_os_nextpow2), ...
    floor(phadir_measured/2+1-size_k_middle*phadir_measured):ceil(phadir_measured/2+1+size_k_middle*phadir_measured)) = 0;


    % %DEBUG MODE: WHICH K_SPACE POINTS GET CHANGED, WHICH NOT?
    % 
    % dont_change = zeros(total_channel_no, fredir_os_nextpow2, phadir_measured);
    % dont_change(:,floor(fredir_os_nextpow2/2+1-size_k_middle*fredir_os_nextpow2):ceil(fredir_os_nextpow2/2+1+size_k_middle*fredir_os_nextpow2),floor(phadir_measured/2+1-size_k_middle*phadir_measured): ...
    % ceil(phadir_measured/2+1+size_k_middle*phadir_measured)) = 1;
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
    % [change_k_channel, change_k_fredir, change_k_phadir] = ind2sub(size(change_k_points),linear_index); 



    image_dat_kspace(logical(change_k_points)) = 0;

    % figure
    % imagesc(abs(squeeze(image_dat_kspace(2,:,:))))

end



%% CORRECTION: 2.2. DO fredir_SHIFT IN K-SPACE


if(exist('fredir_shift','var') && ne(fredir_shift,0))                                                    % shift 1 256 voxels instead of 1 128 voxel; The resolut.
    image_dat_kspace = image_dat_kspace * exp(1i*deg2rad(fredir_shift*(360/fredir_os_nextpow2)));        % does not increase, there is just crap on the borders
end                                                                                                      % so you still have to shift 1 voxel.




%% 3. FFT AND TRUNCATE IMAGE AGAIN, FLIP LEFT AND RIGHT


% TRUNCATE IN KSPACE TO [fredir_desired*oversampling_fredir, phadir_desired], FFT TO DIRECT SPACE
image_dat_kspace2 = image_dat_kspace(:,fredir_left_border_kspace:fredir_right_border_kspace,k_center(2)-phadir_desired/2:k_center(2)-1+phadir_desired/2,:);
image_dat = ifftshift(ifftshift(image_dat_kspace2,2),3);
image_dat = fft(fft(image_dat,[],2),[],3);
image_dat = fftshift(fftshift(image_dat,2),3);


% TRUNCATE IN XSPACE FROM [fredir_desired*oversampling_fredir, phadir_desired] TO [fredir_desired, phadir_desired]
image_dat = reshape(image_dat(:,fredir_left_border_xspace:fredir_right_border_xspace,:), [total_channel_no fredir_desired phadir_desired, 1]);    %truncate left and right borders



%% 4. ROTATE IF IMAGE IS FLIPPED BECAUSE OF REVERSED ENCODING DIRECTIONS; FLIP BECAUSE LEFT RIGHT PHYSICIAN MIX UP
if(exist('phase_encod_dir','var') && strcmp(phase_encod_dir,'RL'))
   
    for channel_no = 1:total_channel_no
        image_dat(channel_no,:,:) = flipdim(rot90(squeeze(image_dat(channel_no,:,:)),-1),2);
    end
    
else
    image_dat = flipdim(image_dat,2);
end


% image_dat = circshift(image_dat, [0 1 0 0]);
pause off