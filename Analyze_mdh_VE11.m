function ParList = Analyze_mdh_VE11(csi_path,AnalyzeWholeMDH)
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
%           -- ParList.total_ADC_meas               - obvious
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

% Read SuperRaidFileHeader
% Read MrParcRaidFileHeader
ID = fread(fid,1, 'uint32');
NoOfMeas = fread(fid,1, 'uint32');

% Read MrParcRaidFileEntries
for CurMeas = 1:NoOfMeas
    MeasID{CurMeas} = fread(fid,1, 'uint32');
    FileID{CurMeas} = fread(fid,1, 'uint32');
    MeasOffset{CurMeas} = fread(fid,1, 'uint64');
    MeasLength{CurMeas} = fread(fid,1, 'uint64');
    PatName{CurMeas} = transpose(fread(fid,64, 'uint8=>char'));
    ProtName{CurMeas} = transpose(fread(fid,64, 'uint8=>char'));
end
ParList.General.LastMeasOffset = MeasOffset{end};


% Goto last measurement, skip all previous ones (they are automatic "pre-scans" like AdjCoilSens, and not of interest right now)
fseek(fid,MeasOffset{end},'bof');

% READ HEADERSIZE
headersize = fread(fid,1, 'uint32');
ParList.General.headersize = headersize;



%% 1. READ FROM FIRST MDH

fseek(fid,headersize-4,'cof');      % -4: Because we read in the headersize of 4 bytes
chak_header = zeros([1 64]);



% Read first mdh
chak_header(1:5) = fread(fid, 5, 'uint32');
fseek(fid,5*4,'cof');

EvalInfoMask = fread(fid, 64, 'ubit1');
chak_header(8:57) = fread(fid, 50, 'int16');
EvalInfoMask_First = Associate_EvalInfoMask(EvalInfoMask);


% VectorSize & Total Channel Number
ParList.(EvalInfoMask_First).Samples = chak_header(8);
ParList.(EvalInfoMask_First).total_channel_no = chak_header(9);



%% 2. READ FROM LAST MDH

% Read last mdh
fseek(fid, -512,'eof');                             % 512 works, but why 512? It should be actually fseek(fid, -224-(512-mod(MeasLength{end},512)),'eof');
chak_header = fread(fid, 3, 'uint32');              % (512-mod(MeasLength{end},512)) to go back the "Alignment to 512 bytes filled with 0", and the -224 going back  
ParList.General.total_ADC_meas = chak_header(3)-1;  % the 32 bytes channel header, and the 192 bytes Scan Header...



%% Find out which Datasets are in the file

if(mod(AnalyzeWholeMDH,2) == 1)
	
	
	
	fseek(fid, headersize+MeasOffset{end},'bof');
	EvalInfoMask_Cur = EvalInfoMask_First;
	step_big = 1;
	step = step_big;
	while(~strcmpi(EvalInfoMask_Cur,'ACQEND') || ~feof(fid))



		fseek_ok = fseek(fid, step*(192 + ParList.(EvalInfoMask_First).total_channel_no*(32 + ParList.(EvalInfoMask_Cur).Samples*2*4)),'cof');
		if(fseek_ok == -1)
			step = 1;
			continue;
		end
		
		
		chak_header(1:5) = fread(fid, 5, 'uint32');
        fseek(fid,5*4,'cof');
		EvalInfoMask = fread(fid, 64, 'ubit1');
		chak_header(8) = fread(fid, 1, 'uint16');		
		fseek(fid, -5*4-5*4-64*1/8-1*2,'cof');


		Assoc = Associate_EvalInfoMask(EvalInfoMask);
		if(strcmp(Assoc,'None') || strcmp(Assoc,'Other'))
			fseek(fid, -step*(64*2 + ParList.(EvalInfoMask_Cur).Samples*2*4),'cof');
			step = 1;
		elseif(strcmp(Assoc,'ACQEND'))
			break;
		elseif(~strcmp(Assoc,EvalInfoMask_Cur))
			EvalInfoMask_Cur = Assoc;
			step = step_big;
			ParList.(EvalInfoMask_Cur).Samples = chak_header(8);
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






%% 3. 

if(AnalyzeWholeMDH > 0)
    
    chak_header = zeros([1 64]);
    fseek(fid,headersize+MeasOffset{end},'bof');
    for i = 1:ParList.General.total_ADC_meas
        
        
        	fseek(fid,+10*4,'cof');	EvalInfoMask_loop = fread(fid,64,'ubit1');            
            Assoc = Associate_EvalInfoMask(EvalInfoMask_loop);
            
            
            if(~isfieldRecursive(ParList,Assoc,'NoOfADCs'))
                ParList.(Assoc).NoOfADCs = 0;
                
                chak_header(8:57) = fread(fid, 50, 'int16');                               
                % Assume these dont change for different ADCs of this type
                ParList.(Assoc).total_channel_no = chak_header(9);
                ParList.(Assoc).Samples = chak_header(8);
                ParList.(Assoc).SamplesAfterEcho = chak_header(8) - chak_header(23)*2;
                fseek(fid,-50*2,'cof');
            end
            fseek(fid,-10*4-8,'cof');
            
            ParList.(Assoc).NoOfADCs = ParList.(Assoc).NoOfADCs+1;
            
            % Jump forward to next mdh
            fseek(fid,192 + ParList.(Assoc).total_channel_no*(32 + ParList.(Assoc).Samples*2*4),'cof');
    end

    for fieldys = transpose(fieldnames(ParList))
        
        CurField = fieldys{:};
        if(strcmpi(CurField,'General'))
            continue
        end
        
        ParList.(CurField).TotalRawPointsPerChannel = ParList.(CurField).NoOfADCs * ParList.(CurField).Samples;
        ParList.(CurField).TotalRawPoints = ParList.(CurField).TotalRawPointsPerChannel * ParList.(CurField).total_channel_no;
    end
    
    
    
    
end






%% 4. READ MATRIX SIZE (ROW AND COLUMN NUMBERS) FROM MAXIMUM K-POINT
% THIS IS LEGACY CODE, NO LONGER MAINTAINED!!!

if(AnalyzeWholeMDH > 3)
    
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
	
    for i = 1:ParList.total_ADC_meas*ParList.Averages
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
            fseek(fid,ParList.Samples*2*4*ParList.total_channel_no + (ParList.total_channel_no-1)*128,'cof');
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
