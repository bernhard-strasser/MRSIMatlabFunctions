function [image_dicom, image_dicom_kspace] = read_image_dicom_1_0(image_dicom_file,fredir_desired, phadir_desired,fredir_shift,total_channel_no,ROW,COL,SLC)
%% 0. Preparations

% if(~exist('zerofilling_fact','var'))
%     zerofilling_fact = 1;
% end
% 
% if(~exist('x_shift','var'))
%     x_shift = 0;
% end
% 
% if(~exist('y_shift','var'))
%     y_shift = 0;
% end


%[ROW,COL,SLC,vecSize,total_channel_no,total_k_points,headersize] = read_csi_dat_meas_header_1_0(csi_dat_file);




% Split Input file

Image_MagPha_file = regexp(image_dicom_file,'\t','split');
%Image_Mag_file = Image_MagPha_file{1};
%Image_Pha_file = Image_MagPha_file{2};



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


