function [ROW,COL,SLC,vecSize,total_channel_no,total_k_points,headersize] = read_csi_dat_meas_header_1_0(csi_dat_file)
%% 0. Preparations


raw_csi_fid = fopen(sprintf('%s', csi_dat_file),'r');

% READ HEADERSIZE
headersize = fread(raw_csi_fid,1, 'uint32');




%% 1. READ VECSIZE, TOTAL CHANNELS FROM FIRST DATA-HEADER

fseek(raw_csi_fid,headersize,'bof'); 
chak_header = fread(raw_csi_fid, 64, 'uint16');

vecSize = chak_header(15);

total_channel_no = chak_header(16);



%% 2. READ TOTAL k-points AND NUMBER OF SLICES FROM END OF FILE
fseek(raw_csi_fid, -256,'eof');
chak_header = fread(raw_csi_fid, 23, 'uint16');

total_k_points = chak_header(5)-1;

if(chak_header(23) > 0)
    SLC = chak_header(23);
else
    SLC = 1;
end


%% 3. READ MATRIX SIZE (ROW AND COLUMN NUMBERS) FROM MAXIMUM K-POINT
fseek(raw_csi_fid,headersize,'bof');
ROW = 1;
COL = 1;
for i = 1:total_k_points

        chak_header = fread(raw_csi_fid, 64, 'uint16');
        if(chak_header(17) > ROW)
            ROW = chak_header(17);
        end
        if(chak_header(22) > COL)
            COL = chak_header(22);
        end
        fseek(raw_csi_fid,vecSize*2*4*total_channel_no + (total_channel_no-1)*128,'cof');

end

ROW = ROW + 1;
COL = COL + 1;

fclose(raw_csi_fid);
