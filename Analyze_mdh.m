function ParList = Analyze_csi_mdh(csi_path,AnalyzeWholeMDH)
%
% Analyze_csi_mdh_x_x Analyze measurement data header of Siemens csi raw data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function analyzes the mdh of Siemens raw csi data to find out the Parameters listed under Output
%
%
% ParList = read_Analyze_csi_mdh_1_3(csi_path, AnalyzeWholeMDH)
%
% Input: 
% -         csi_path                        ...     Path of MRS(I) file.
% -         AnalyzeWholeMDH         ...     Determines if over whole k-space is looped to find out the matrix sizes ROW and COL.
%
% Output:
% -         ParList                         ...     Structure giving all the Parameters. It contains:
%           -- ParList.vecSize                      - VectorSize in spectroscopic dimension
%           -- ParList.total_channel_no             - number of receive-channels
%           -- ParList.total_k_points               - obvious
%           -- ParList.SLC                          - Number of Slices
%           -- ParList.Averages                     - Number of acquired averages.
% if AnalyzeWholeMDH = true
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
if(~exist('AnalyzeWholeMDH','var'))
    AnalyzeWholeMDH = 0;
end

% Open file
fid = fopen(sprintf('%s', csi_path),'r');

% READ HEADERSIZE
headersize = fread(fid,1, 'uint32');




%% 1. READ FROM FIRST DATA-HEADER

fseek(fid,headersize,'bof');
chak_header = zeros([1 64]);



% Read first mdh
chak_header(1:5) = fread(fid, 5, 'uint32');
EvalInfoMask = fread(fid, 64, 'ubit1');
chak_header(8:57) = fread(fid, 50, 'int16');
EvalInfoMask_First = Associate_EvalInfoMask(EvalInfoMask);


% VectorSize & Total Channel Number
ParList.(EvalInfoMask_First).vecSize = chak_header(8);




%% 2. READ FROM END OF FILE

% Read last mdh
fseek(fid, -256,'eof');
chak_header = fread(fid, 64, 'uint16');

ParList.General.total_ADC_meas = chak_header(5)-1;



%% Find out which Datasets are in the file

if(mod(AnalyzeWholeMDH,2) == 1)
	
	
	
	fseek(fid, headersize,'bof');
	EvalInfoMask_Cur = EvalInfoMask_First;
	step = 12;
	while(~strcmpi(EvalInfoMask_Cur,'ACQEND') || ~feof(fid))



		fseek_ok = fseek(fid, step*(64*2 + ParList.(EvalInfoMask_Cur).vecSize*2*4),'cof');
		if(fseek_ok == -1)
			step = 1;
			continue;
		end
		
		
		
		chak_header(1:7) = fread(fid, 7, 'uint32');
		chak_header(8) = fread(fid, 1, 'uint16');		
		fseek(fid, -7*4-1*2,'cof');



		if(strcmp(Associate_EvalInfoMask(Interpret_EvalInfoMask(chak_header(6:7))),'0'))
			fseek(fid, -step*(64*2 + ParList.(EvalInfoMask_Cur).vecSize*2*4),'cof');
			step = 1;
		elseif(strcmp(Associate_EvalInfoMask(Interpret_EvalInfoMask(chak_header(6:7))),'ACQEND'))
			break;
		elseif(~strcmp(Associate_EvalInfoMask(Interpret_EvalInfoMask(chak_header(6:7))),EvalInfoMask_Cur))
			EvalInfoMask_Cur = Associate_EvalInfoMask(Interpret_EvalInfoMask(chak_header(6:7)));
			step = 12;
			ParList.(EvalInfoMask_Cur).vecSize = chak_header(8);
		end


	end
	
	
	
end



%% 3. Find other Datasets inside file

% if(strcmp(EvalInfoMask_First,EvalInfoMask_Last))
% 	ParList.DataSetNames = EvalInfoMask_First;
% else
% 	ParList.DataSetNames{1} = EvalInfoMask_First;
% 	ParList.DataSetNames{2} = EvalInfoMask_Last;
% 	
% 
% 	
% 	
% 	
% end



%% 3. READ MATRIX SIZE (ROW AND COLUMN NUMBERS) FROM MAXIMUM K-POINT

if(AnalyzeWholeMDH > 1)
    
    fseek(fid,headersize,'bof');
    kline_max = 0; kline_min = 99999; 
	kphase_max = 0; kphase_min = 99999;
	kpart_max = 0; kpart_min = 99999;
	slice_max = 0; slice_min = 99999;
	echo_max = 0; echo_min = 99999;
	avg_max = 0; avg_min = 99999;
	rep_max = 0; rep_min = 99999;
	
	SampleList = zeros([1 10^5]);
	InfoList = zeros([8 10^5]);
	
    for i = 1:ParList.total_k_points*ParList.Averages
            chak_header = fread(fid, 64, 'uint16');
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
            fseek(fid,ParList.vecSize*2*4*ParList.total_channel_no + (ParList.total_channel_no-1)*128,'cof');
    end
    ParList.ROW_measured = kline_max - kline_min + 1;
    ParList.COL_measured = kphase_max - kphase_min + 1;
    ParList.kline_min = kline_min + 1;
    ParList.kline_max = kline_max + 1;
    ParList.kphase_min = kphase_min + 1;
    ParList.kphase_max = kphase_max + 1;
    
end





%% 4. POSTPARATIONS

fclose(fid);
