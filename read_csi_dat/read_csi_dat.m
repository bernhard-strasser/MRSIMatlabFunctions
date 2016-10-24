function [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% read_csi_dat Read in raw data from Siemens
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         file                    ...     Path of MRS(I) file.
% -         DesiredSize             ...     If you want to perform zerofilling or use only a part of the kspace, set
%                                           DesiredSize to your wanted [kx,ky,kz], e.g. [64 64 1].
% -         ReadInDataSets          ...     
%
% Output:
% -         kSpace                      ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x Samples x Averages
% -         Info                        ...     
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.





%% 0. Preparations


% Find out memory used by MATLAB
memused_before = memused_linux(1); 

if(exist('DesiredSize','var') && numel(DesiredSize) == 1 && DesiredSize == 0)
	clear DesiredSize;
end
if(exist('DesiredSize','var') && ~isstruct(DesiredSize))
	DesiredSize2 = DesiredSize; clear DesiredSize; DesiredSize.ONLINE = DesiredSize2; clear DesiredSize2;
end
if(~exist('ReadInDataSets','var'))
	ReadInDataSets = 'All';
end

Info.mdhEntryNames = {'channel','kx','ky','kz','slice','echo','avg','rep','samples','samples before echo','ida','idb','idc','idd','FreeIcePara1','FreeIcePara2','FreeIcePara3','FreeIcePara4'};



%% 1. Gather information from header



ParList = read_ascconv(file);
Info.General.Ascconv = ParList;




% Info about sizes
% First try: Get info from wipMemBlock_tfree
if(ParList.wipMemBlock_tFree ~= 0)
	% Only evaluate the tFree string if it is proper MATLAB CODE (doesnt throw any errors)
	try
		eval(ParList.wipMemBlock_tFree);
	catch
		fprintf(['\nThe wipMemBlock.tFree entry is not MATLAB code. Consider writing the info about\n' ...
		'the sizes of those scans in the wipMemBlock.tFree!\n\n'])
		ParList.wipMemBlock_tFree = 0;
	end
end

% Second try: Try to get Prescans Info from wipMemBlock and ONLINE Info from normal ascconv header and mdh.
% Prescan Info
if(isfieldRecursive(ParList,'WipMemBlockInterpretation','Prescan'))
	if(isfieldRecursive(ParList.WipMemBlockInterpretation,'NOISEADJSCAN','Dwelltime') && ParList.WipMemBlockInterpretation.Prescan.NOISEADJSCAN.Dwelltime > 0)
		Info.NOISEADJSCAN = ParList.WipMemBlockInterpretation.Prescan.NOISEADJSCAN;
	end
	if(isfield(ParList.WipMemBlockInterpretation,'ONLINE'))
		Info.ONLINE = ParList.WipMemBlockInterpretation.Prescan.ONLINE;
	end
	if(isfieldRecursive(ParList.WipMemBlockInterpretation.Prescan,'PATREFANDIMASCAN','Dwelltime') && ParList.WipMemBlockInterpretation.Prescan.PATREFANDIMASCAN.Dwelltime > 0)
		Info.PATREFANDIMASCAN = ParList.WipMemBlockInterpretation.Prescan.PATREFANDIMASCAN;
	end
end
	

% ONLINE Info
if(~isfield(ParList,'wipMemBlock_tFree') || (isfield(ParList,'wipMemBlock_tFree') && ParList.wipMemBlock_tFree == '0'))

	Info.ONLINE.nReadEnc = ParList.nFreqEnc;
	Info.ONLINE.nPhasEnc = ParList.nPhasEnc;
	Info.ONLINE.nPartEnc = ParList.nPartEnc;
	Info.ONLINE.nReadEncFinalMatrix = ParList.nFreqEnc_FinalMatrix;
	Info.ONLINE.nPhasEncFinalMatrix = ParList.nPhasEnc_FinalMatrix;
	%Info.ONLINE.nPartEncFinalMatrix = ParList.nSLC_FinalMatrix;
	if(Info.ONLINE.nReadEncFinalMatrix == 0)
		Info.ONLINE.nReadEncFinalMatrix = Info.ONLINE.nReadEnc;
	end
	if(Info.ONLINE.nPhasEncFinalMatrix == 0)
		Info.ONLINE.nPhasEncFinalMatrix = Info.ONLINE.nPhasEnc;
	end	
% 	if(Info.ONLINE.nPartEncFinalMatrix == 0)
% 		Info.ONLINE.nPartEncFinalMatri = Info.ONLINE.nPartEnc;
% 	end	
	
	Info.ONLINE.kSpaceShiftDueToICEZeroFill = zeros([1 5]);
	
% This is legacy code... Probably the Zerofilling should be done by PreProcessMRIData_Wrapper only, not here and there... ?!
% 	if(???)
% 		Info.ONLINE.kSpaceShiftDueToICEZeroFill(1) = floor(Info.ONLINE.nReadEncFinalMatrix/2) - floor(Info.ONLINE.nReadEnc/2);
% 		Info.ONLINE.kSpaceShiftDueToICEZeroFill(2) = floor(Info.ONLINE.nPhasEncFinalMatrix/2) - floor(Info.ONLINE.nPhasEnc/2);
% 		%Info.ONLINE.kSpaceShiftDueToICEZeroFill(3) = floor(Info.ONLINE.nPartEncFinalMatrix/2) - floor(Info.ONLINE.nPartEnc/2);
% 	end
	Info.ONLINE.nAverages = ParList.nAverages;
	Info.ONLINE.Dwelltime = ParList.Dwelltimes(1);	% For now only take the first one. Assume that all slices have the same dwelltime
	
	if(isfieldRecursive(ParList,'WipMemBlockInterpretation','OneDCaipi','NoOfMeasSlices'))
		Info.ONLINE.nSLC = ParList.WipMemBlockInterpretation.OneDCaipi.NoOfMeasSlices;
	else
		Info.ONLINE.nSLC = ParList.nSLC;
	end
end


% Analyze mdh.
ParList = Analyze_mdh(file,0);
%Info.ONLINE.total_channel_no = ParList.total_channel_no;
Info.General.total_ADCs = ParList.General.total_ADC_meas;
clear ParList;


% Otherwise read data without memory preallocation. This can be slow.




%% 2. READ DATA

tic
fprintf('\nReading data\t\t\t...')
        
file_fid = fopen(sprintf('%s', file),'r');
headersize = fread(file_fid,1, 'uint32');
fseek(file_fid, headersize,'bof'); 

chak_header = zeros([1 64]);

kSpace.ONLINE = NaN;
ACQEND_flag = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%	Loop over different measurement sets	%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SkipMeassageAlreadyPrinted = false;
while(~ACQEND_flag)

    % Read first mdh
    chak_header(1:5) = fread(file_fid, 5, 'uint32');
    EvalInfoMask = fread(file_fid, 64, 'ubit1');
	CurChak = fread(file_fid, 64-14, 'int16');
    fseek(file_fid, -128,'cof');
    
    % Stop if ACQEND was found
	if(EvalInfoMask(1) || isempty(CurChak))
        ACQEND_flag = true;
        continue;
	end
    
	chak_header(8:64-7) = CurChak;

	
    % Set CurrentMeasSet
	CurrentMeasSet = Associate_EvalInfoMask(EvalInfoMask);
	
	if((sum(strcmpi(ReadInDataSets,'All')) || sum(strcmpi(CurrentMeasSet,ReadInDataSets))) )
		fprintf('\nRead\t%s data.', CurrentMeasSet)
	else
		if(~SkipMeassageAlreadyPrinted)
			fprintf('\nSkip\t%s data.', CurrentMeasSet)
			SkipMeassageAlreadyPrinted = true;
		end
		fseek(file_fid,(128+chak_header(8)*2*4)*chak_header(9),'cof');
		continue;
	end	
	SkipMeassageAlreadyPrinted = false;
	
	% Initialize Current Info
	if(~isfield(Info,CurrentMeasSet))
		Info.(CurrentMeasSet) = struct('nReadEnc',1,'nPhasEnc',1,'nPartEnc',1,'nSLC',1, 'nAverages',1);
	end
	
	dim = 0;
	for fieldly = {'nReadEnc','nPhasEnc','nPartEnc','nSLC','nAverages'}
		dim = dim + 1;
		% Initialize missing fields with ones
		if(~isfield(Info.(CurrentMeasSet), fieldly{:}))
			Info.(CurrentMeasSet).(fieldly{:}) = 1;
			DoShift = false;
		else
			DoShift = true;
		end		
		if(DoShift && exist('DesiredSize','var') && isfield(DesiredSize,CurrentMeasSet) && numel(DesiredSize.(CurrentMeasSet)) >= dim)
			kSpaceShift.(CurrentMeasSet)(dim) = floor(Info.(CurrentMeasSet).(fieldly{:})/2) - floor(DesiredSize.(CurrentMeasSet)(dim)/2);
			Info.(CurrentMeasSet).(fieldly{:}) = DesiredSize.(CurrentMeasSet)(dim);			
		else
			kSpaceShift.(CurrentMeasSet)(dim) = 0;
			DesiredSize.(CurrentMeasSet)(dim) = 99999;		% So that it has no effect			
		end
		if(isfield(Info.(CurrentMeasSet),'kSpaceShiftDueToICEZeroFill') && Info.(CurrentMeasSet).kSpaceShiftDueToICEZeroFill(dim) > 0)
			kSpaceShift.(CurrentMeasSet)(dim) = kSpaceShift.(CurrentMeasSet)(dim) - Info.(CurrentMeasSet).kSpaceShiftDueToICEZeroFill(dim);
		end
	end
	

	Info.(CurrentMeasSet).total_channel_no = chak_header(9);
	Info.(CurrentMeasSet).Samples = chak_header(8) - chak_header(38) *2;
	
	if( (strcmpi(CurrentMeasSet,'ONLINE') || strcmpi(CurrentMeasSet,'NOISEADJSCAN')) && Info.General.Ascconv.vecSize > 1 )
		vecSize = Info.(CurrentMeasSet).Samples;
	else
		vecSize = 1;
		Info.(CurrentMeasSet).nReadEnc = Info.(CurrentMeasSet).Samples;
	end
	
	
	% Allocate memory
	if( ~isfield(kSpace,CurrentMeasSet) || isnan(kSpace.(CurrentMeasSet)) )					% By this statement, interleaved Prescan measurements are not allowed (data will be overwritten)
		kSpace.(CurrentMeasSet) = zeros([Info.(CurrentMeasSet).total_channel_no,Info.(CurrentMeasSet).nReadEnc * Info.(CurrentMeasSet).nPhasEnc * ...
		Info.(CurrentMeasSet).nPartEnc * Info.(CurrentMeasSet).nSLC * vecSize * Info.(CurrentMeasSet).nAverages]);
	end
	Info.(CurrentMeasSet).mdhInfo = zeros([18 Info.General.total_ADCs*Info.(CurrentMeasSet).total_channel_no]);
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%	Loop over all measurements	%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BreakOutOfPrison = false;
	CurPoint = 0;
	for ADC_MeasNo = 1 : Info.General.total_ADCs;
		
		if(BreakOutOfPrison)
			break;
		end
		
		for channel_no = 1:Info.(CurrentMeasSet).total_channel_no
			
			% Read again EvalInfoMask.
			fseek(file_fid,+5*4,'cof');	EvalInfoMask_loop = fread(file_fid,64,'ubit1');
			CurInfoPt = (ADC_MeasNo-1)*Info.(CurrentMeasSet).total_channel_no + channel_no;
			
			% If the EvalInfoMask has changed within the loop, break out of loops.
			CurrentMeasSet_loop = Associate_EvalInfoMask(EvalInfoMask_loop);
			if(~strcmpi(CurrentMeasSet,CurrentMeasSet_loop))
				fseek(file_fid,-7*4,'cof');
				Info.(CurrentMeasSet).mdhInfo(:,CurInfoPt:end) = [];
				BreakOutOfPrison = true;
				break;
			end
			
			% Read rest of the mdh
			chak_header(8:57) = fread(file_fid, 50, 'int16');
			
			Info.(CurrentMeasSet).mdhInfo(1,CurInfoPt) = chak_header(56) + 1;										% channel
			Info.(CurrentMeasSet).mdhInfo(2,CurInfoPt) = chak_header(10) + 1 - kSpaceShift.(CurrentMeasSet)(1);		% kx
			Info.(CurrentMeasSet).mdhInfo(3,CurInfoPt) = chak_header(15) + 1 - kSpaceShift.(CurrentMeasSet)(2);		% ky
			Info.(CurrentMeasSet).mdhInfo(4,CurInfoPt) = chak_header(13) + 1 - kSpaceShift.(CurrentMeasSet)(3);		% kz
			Info.(CurrentMeasSet).mdhInfo(5,CurInfoPt) = chak_header(12) + 1;										% slice (also hada step)
			Info.(CurrentMeasSet).mdhInfo(6,CurInfoPt) = chak_header(14) + 1;										% echo
			Info.(CurrentMeasSet).mdhInfo(7,CurInfoPt) = chak_header(17) + 1;										% avg
			Info.(CurrentMeasSet).mdhInfo(8,CurInfoPt) = chak_header(16) + 1;										% rep
			Info.(CurrentMeasSet).mdhInfo(9,CurInfoPt) = chak_header(8);											% samples
			Info.(CurrentMeasSet).mdhInfo(10,CurInfoPt) = chak_header(38) *2;										% samples before echo
			Info.(CurrentMeasSet).mdhInfo(11,CurInfoPt) = chak_header(19) + 1;										% ida
			Info.(CurrentMeasSet).mdhInfo(12,CurInfoPt) = chak_header(20) + 1;										% idb
			Info.(CurrentMeasSet).mdhInfo(13,CurInfoPt) = chak_header(21) + 1;										% idc
			Info.(CurrentMeasSet).mdhInfo(14,CurInfoPt) = chak_header(22) + 1;										% idd
			Info.(CurrentMeasSet).mdhInfo(15,CurInfoPt) = chak_header(34);											% FreeIcePara1
			Info.(CurrentMeasSet).mdhInfo(16,CurInfoPt) = chak_header(35);											% FreeIcePara2	
			Info.(CurrentMeasSet).mdhInfo(17,CurInfoPt) = chak_header(36);											% FreeIcePara3
			Info.(CurrentMeasSet).mdhInfo(18,CurInfoPt) = chak_header(37);											% FreeIcePara4	

			% Read & Assign Data	% Read real & imaginary (--> Info.(CurrentMeasSet).Samples*2) measured points
            chak_data = fread(file_fid, Info.(CurrentMeasSet).mdhInfo(9,CurInfoPt)*2, 'float32'); 
			chak_data = chak_data(Info.(CurrentMeasSet).mdhInfo(10,CurInfoPt)*2+1:end); 
			chak_data = complex(chak_data(1:2:end),chak_data(2:2:end));
			kSpace.(CurrentMeasSet)(channel_no,CurPoint+1 : CurPoint+numel(chak_data)) = chak_data;

			
			
			% Check if this was the last measurement of scan
			% Temporarily disabled because the flag is set for all HadamardSteps. Has to be changed.
% 			if(EvalInfoMask_loop(9))		% LASTSCANINMEAS
% 				BreakOutOfPrison = true;
% 				break;
% 			end		
			

		end 
		CurPoint = CurPoint + numel(chak_data);
	end

	
	
end

if(numel(kSpace.ONLINE) == 1 && isnan(kSpace.ONLINE))
	kSpace = rmfield(kSpace,'ONLINE');
end
    
fclose(file_fid);

fprintf('\n\t\t\t\t...took\t%10.6f seconds',toc)       

if(numel(fields(kSpace)) > 1 && numel(fields(Info)) + 1 < numel(fields(kSpace)) )
	fprintf(['\nNo wipMemBlock.tFree and wipMemBlock.alFree[50-55] entries found. If you have several datasets\n(like Prescans) in your raw data, consider writing the info about\n' ...
	'the sizes of those scans in the wipMemBlock.tFree or wipMemBlock.alFree[50-55]!\n\n'])
end
%fprintf('\n\n\nChange ''Interpret_EvalInfoMask.m'' function to rename data sets.')



%% 3. Reshape Data

tic
fprintf('\n\nReshaping data\t\t\t...')

for CurrentMeasSet2 = transpose(fields(kSpace))
	CurrentMeasSet = CurrentMeasSet2{1};
	fprintf('\nReshape\t%s data.', CurrentMeasSet)

	% Test if each ADC has unique mdh. If not, write warning.
	if(size(unique(transpose(Info.(CurrentMeasSet).mdhInfo),'rows'),1) < size(Info.(CurrentMeasSet).mdhInfo,2))
		fprintf('\n\nWARNING:\t%s has non-unique mdh entries for different ADCs. Data will be overwritten.\n\n', CurrentMeasSet)	
	end
	
	% What's this for?
	if(sum(  Info.(CurrentMeasSet).mdhInfo(6,:) ~= Info.(CurrentMeasSet).mdhInfo(5,:)  ) == 0)
		Info.(CurrentMeasSet).mdhInfo(6,:) = 1;
	end
	
	maxi = max(transpose(Info.(CurrentMeasSet).mdhInfo),[],1);
	Temp = cell([1 maxi(6)]);
	mdhInfo_reshaped = cell([1 maxi(6)]);
	
	% Find out Samples of each echo
	minecho = min(Info.(CurrentMeasSet).mdhInfo(6,:));
	samples = zeros([1 maxi(6)-minecho+1]);
	for echo = minecho : maxi(6)
		echoplace = find(Info.(CurrentMeasSet).mdhInfo(6,:) == echo);
		echoplace = echoplace(1);
		samples(echo) = Info.(CurrentMeasSet).mdhInfo(9,echoplace) - Info.(CurrentMeasSet).mdhInfo(10,echoplace);
	end
	
	% Correct for the problem in 2D-GRAPPA and certain 2D-CAIPI Patterns, where the border voxels are not measured and therefore the matrix has a wrong size.
	if(isfieldRecursive(Info,'General','Ascconv','WipMemBlockInterpretation','TwoDCaipi','Skip_Matrix') && strcmp(CurrentMeasSet,'ONLINE'))
		maxi(2:4) = [Info.ONLINE.nReadEnc Info.ONLINE.nPhasEnc Info.ONLINE.nPartEnc];
	end
	
	% Correct for the problem that in circular encoding e.g. 63 points are measured, but the matrix should be 64!
	for dim = 2:4
		if(2^nextpow2(maxi(dim)) - maxi(dim) == 1)
			maxi(dim) = maxi(dim) + 1;
		end
	end
	
	% Initialize
	for echo = minecho:maxi(6)
		Temp{echo} = zeros([Info.(CurrentMeasSet).total_channel_no maxi(2) maxi(3) maxi(4) maxi(5) samples(echo) maxi(7) maxi(8) maxi(11) maxi(12)]);
		mdhInfo_reshaped{echo} = zeros([1 maxi(2) maxi(3) maxi(4) maxi(5) 1 maxi(7) maxi(8) maxi(11) maxi(12) size(Info.(CurrentMeasSet).mdhInfo,1)]);	% Reshape mdhInfo: Before: NoOfADCs x 18
	end																										% Now: {echo}[cha x kx x ky x kz x slc x samples x avg x rep x ADCNo x TempIntNo x 18]

	CurPoint = 0; 
	for i = 1 : Info.(CurrentMeasSet).total_channel_no : size(Info.(CurrentMeasSet).mdhInfo,2)
		% {echo}[cha x kx x ky x kz x slc x samples x avg x rep x ADCNo x TempIntNo]
		CurEco = Info.(CurrentMeasSet).mdhInfo(6,i); Curkx = Info.(CurrentMeasSet).mdhInfo(2,i); Curky = Info.(CurrentMeasSet).mdhInfo(3,i); Curkz = Info.(CurrentMeasSet).mdhInfo(4,i);
		CurSlc = Info.(CurrentMeasSet).mdhInfo(5,i); CurAvg = Info.(CurrentMeasSet).mdhInfo(7,i); CurRep = Info.(CurrentMeasSet).mdhInfo(8,i);
		CurADC = Info.(CurrentMeasSet).mdhInfo(11,i); CurTempIntNo = Info.(CurrentMeasSet).mdhInfo(12,i);
		
		Temp{CurEco}(:,Curkx, Curky, Curkz, CurSlc, :, CurAvg, CurRep, CurADC,CurTempIntNo) = ...
		kSpace.(CurrentMeasSet)(:,CurPoint+1:CurPoint+Info.(CurrentMeasSet).mdhInfo(9,i)-Info.(CurrentMeasSet).mdhInfo(10,i));
	
		mdhInfo_reshaped{CurEco}(1,Curkx, Curky, Curkz, CurSlc, 1, CurAvg, CurRep, CurADC,CurTempIntNo,:) = Info.(CurrentMeasSet).mdhInfo(:,i);
	
		CurPoint = CurPoint + (Info.(CurrentMeasSet).mdhInfo(9,i)-Info.(CurrentMeasSet).mdhInfo(10,i));
	end
	Info.(CurrentMeasSet).mdhInfo = mdhInfo_reshaped; clear mdhInfo_reshaped; 
	
	for echo = minecho:maxi(6)
		if(size(Temp{echo},3) == 1 && (size(Temp{echo},2) > 1 || size(Temp{echo},4) > 1) )			% This should basically mean: If data is imaging data. However it is a little cheated
			Temp{echo} = permute(Temp{echo},[1 6 2 4 5 3 7 8 9 10]);
		end
	end
	
	kSpace.(CurrentMeasSet) = Temp;		
	
end
fprintf('\n\t\t\t\t...took\t%10.6f seconds',toc)       




%% 7. Postparations

memused_after = memused_linux(1); 
display([char(10) 'The function used ' num2str(memused_after-memused_before) '% of the total memory.'])


