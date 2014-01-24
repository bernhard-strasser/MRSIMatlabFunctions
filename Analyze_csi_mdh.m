function [headersize,total_channel_no,ROW,COL,SLC,vecSize,total_k_points,kline_min,kphase_min] = read_csi_dat_meas_header_1_2(csi_dat_file)
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

if(chak_header(19) > 0)
    SLC = chak_header(19)+1;
else
    SLC = 1;
end


%% 3. READ MATRIX SIZE (ROW AND COLUMN NUMBERS) FROM MAXIMUM K-POINT
fseek(raw_csi_fid,headersize,'bof');
kline_max = 0; kphase_max = 0; kline_min = 99999; kphase_min = 99999; 

for i = 1:total_k_points
        chak_header = fread(raw_csi_fid, 64, 'uint16');
        if(chak_header(17) > kline_max)
            kline_max = chak_header(17);
        end
        if(chak_header(17) < kline_min)
            kline_min = chak_header(17);
        end        
        
        if(chak_header(22) > kphase_max)
            kphase_max = chak_header(22);
        end
        if(chak_header(22) < kphase_min)
            kphase_min = chak_header(22);
        end        
        fseek(raw_csi_fid,vecSize*2*4*total_channel_no + (total_channel_no-1)*128,'cof');
end

ROW = kline_max - kline_min + 1;
COL = kphase_max - kphase_min + 1;

if(ROW > 1 && ~isa(log2(ROW),'integer'))
    ROW = ROW + 1;
end
if(COL > 1 && ~isa(log2(COL),'integer'))
    COL = COL + 1;
end



%% 4. POSTPARATIONS

fclose(raw_csi_fid);
