function ParList = Analyze_image_mdh(image_path)
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

% Furhter remarks: The number of the beast is always greater than the number you think of.
% Die Vereinigung abzählbar unendlich vieler abzählbar unendlicher Mengen ist abzählbar unendlich (Karl Auinger).



%% 0. Preparations

% fredir = frequency encoding direction; phadir = phase encoding direction


% Open file
raw_image_fid = fopen(sprintf('%s', image_path),'r');

% READ SIZE OF HEADER
headersize = fread(raw_image_fid,1, 'int32');



%% 1. Read from first mdh

fseek(raw_image_fid, headersize,'bof');
chak_header = fread(raw_image_fid, 39, 'uint16');


ParList.fredir_measured = chak_header(15);      % Number of data points in freq encod direction
ParList.total_channel_no = chak_header(16);
ParList.k_center_fredir = chak_header(33)+1;    % Center of k-space in frequency encoding direction
ParList.phadir_measured = chak_header(39)*2;    % Number of measured lines in phase encoding direction



%% 2. Read from last mdh

fseek(raw_image_fid, -256,'eof');
chak_header = fread(raw_image_fid, 64, 'uint16');


ParList.SLC = chak_header(19)+1;
ParList.Averages = chak_header(18) + 1;
ParList.total_k_points = (chak_header(5)-1) / ParList.Averages;     % Number of total k-points
ParList.total_ADC_meas = chak_header(5)-1;                          % Number of total ADC-measurements





%% 3. Template for Analyzing every mdh.

% phadir_measured = 1;
% fseek(raw_image_fid, headersize,'bof');
% while(~isempty(chak_header))
%     chak_header = fread(raw_image_fid, 64, 'uint16');
%     if(chak_header(39) > phadir_measured)
%         phadir_measured = chak_header(39);
%     end
%     fseek(raw_image_fid,fredir_measured*2*4*total_channel_no + 0*(total_channel_no-1)*128,'cof');
% end






%% 4. Postparations

fclose(raw_image_fid);

