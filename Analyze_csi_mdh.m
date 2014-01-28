function ParList = Analyze_csi_mdh(csi_path,AnalyzeWholekSpace_flag)
%
% Analyze_csi_mdh_x_x Analyze measurement data header of Siemens csi raw data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function analyzes the mdh of Siemens raw csi data to find out the Parameters listed under Output
%
%
% ParList = read_Analyze_csi_mdh_1_3(csi_path, AnalyzeWholekSpace_flag)
%
% Input: 
% -         csi_path                        ...     Path of MRS(I) file.
% -         AnalyzeWholekSpace_flag         ...     Determines if over whole k-space is looped to find out the matrix sizes ROW and COL.
%
% Output:
% -         ParList                         ...     Structure giving all the Parameters. It contains:
%           -- ParList.vecSize                      - VectorSize in spectroscopic dimension
%           -- ParList.total_channel_no             - number of receive-channels
%           -- ParList.total_k_points               - obvious
%           -- ParList.SLC                          - Number of Slices
%           -- ParList.Averages                     - Number of acquired averages.
% if AnalyzeWholekSpace_flag = true
%           ++ ParList.ROW_measured                 - Number of measured rows (lines), these gives the number of the REALLY measured lines
%           ++ ParList.COL_measured                 - Number of measured columns (phase_encoding)
%           ++ ParList.kline_min                    - Minimum value of k-point in ROW-direction. All values start from 1, not from 0 like in C!
%           ++ ParList.kline_max                    - Maximum value of k-point in ROW-direction.
%           ++ ParList.kphase_min                   - obvious
%           ++ ParList.kphase_max                   - obvious
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None





%% 0. Preparations

% Assign standard values to variables if nothing is passed to function.
if(~exist('AnalyzeWholekSpace_flag','var'))
    AnalyzeWholekSpace_flag = true;
end

% Open file
raw_csi_fid = fopen(sprintf('%s', csi_path),'r');

% READ HEADERSIZE
headersize = fread(raw_csi_fid,1, 'uint32');




%% 1. READ FROM FIRST DATA-HEADER

fseek(raw_csi_fid,headersize,'bof'); 
chak_header = fread(raw_csi_fid, 64, 'uint16');

% VectorSize & Total Channel Number
ParList.vecSize = chak_header(15);




%% 2. READ FROM END OF FILE

fseek(raw_csi_fid, -256,'eof');
chak_header = fread(raw_csi_fid, 64, 'uint16');

% total measured k-space points & Number of Slices

if(chak_header(19) > 0)
    ParList.SLC = chak_header(19)+1;
else
    ParList.SLC = 1;
end
ParList.Averages = chak_header(24) + 1;
ParList.total_k_points = (chak_header(5)-1) / ParList.Averages;			% This is wrong! If there are more ADCs in the file (e.g. from Prescans) the kPoints cannot be computed like that!
ParList.total_ADC_meas = chak_header(5)-1;
ParList.total_channel_no = chak_header(16);




%% 3. READ MATRIX SIZE (ROW AND COLUMN NUMBERS) FROM MAXIMUM K-POINT

if(AnalyzeWholekSpace_flag)
    
    fseek(raw_csi_fid,headersize,'bof');
    kline_max = 0; kphase_max = 0; kline_min = 99999; kphase_min = 99999; 

    for i = 1:ParList.total_k_points*ParList.Averages
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
            fseek(raw_csi_fid,ParList.vecSize*2*4*ParList.total_channel_no + (ParList.total_channel_no-1)*128,'cof');
    end
    ParList.ROW_measured = kline_max - kline_min + 1;
    ParList.COL_measured = kphase_max - kphase_min + 1;
    ParList.kline_min = kline_min + 1;
    ParList.kline_max = kline_max + 1;
    ParList.kphase_min = kphase_min + 1;
    ParList.kphase_max = kphase_max + 1;
    
end





%% 4. POSTPARATIONS

fclose(raw_csi_fid);
