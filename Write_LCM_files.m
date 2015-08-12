function Write_LCM_files(InArray,Paths,MetaInfo,ControlInfo,mask,CPU_cores)
%
% Write_LCM_files Write the files necessary for LCModel fittting.
%
% This function was written by Bernhard Strasser, 2011 and 2012, revised August 2012.
%
%
% This function creates the files that are necessary to start the LCModel fitting from the terminal. It creates the .control files with LCModel
% instructions how to fit the data, and the .RAW files containing the FIDs. The function also creates bash-script files which contains
% the bash-commands to start the LCModel fitting sequentially for each FID of the InArray. The variable 'CPU_cores' determines
% to how many bash-script files these commands are splitted. If you execute all these files at the same time (e.g. by using 
% "bash bash_file_1 & bash bash_file_2 ...") the fitting process is split to several cores (LCModel normally uses only 1 core).
% 
%
% Write_LCM_files(InArray,Paths,MetaInfo,ControlInfo,mask,CPU_cores)
%
% InArray:              Complex array with time domain data, e.g. [64x64x1x2048] CSI Matrix
%						OR
%                       Struct with fields .csi and .watref
% Paths:                Struct with Path info:
%                       -   Paths.out_dir 
%                       -   Paths.basis_file 
%                       -   Paths.LCM_ProgramPath 
%                       -   Paths.batchdir: Folder to which the bash scripts to start LCModel processing should be written
% MetaInfo:             Struct with metainfo about InArray: 
%                       -   MetaInfo.DimNames:       The dimension names of InArray. E.g. {'x','y','z'}. LCModel is instructed to create the output files as 
%                                                    "pat_name_DimNames1a_DimNames2b_DimNames3c" where a,b,c are numbers (e.g. the voxel [32 32 1]).
%                       -   MetaInfo.Dimt1:          Tells the function in which dimension of InArray the FID's are saved, i.e. which dimension is the
%                                                    t1-dimension / vecSize-dimension
%                       -   MetaInfo.pat_name:       The patient name. Used for the output file names
%                       -   MetaInfo.LarmorFreq:     The Larmor frequency which should be read out of the DICOM or raw-data header 
%                       -   MetaInfo.dwelltime:      The dwelltime, i.e. the duration between two consecutive time points.
% ControlInfo:          Struct or a string referring .m file with control parameters for LCModel fitting. This file / struct is executed at the end of the Control-Parameter settings,
%                       and can thus overwrite ALL control parameters. Be careful!
%						You can also pass over a path with control parameters and additional parameters. The path is passed over by field .Path. The other parameters normal, e.g. .Others1 = ...
%						In this case, specify which one should have priority by providing
%						ControlInfo.Priority = 'Fields' or [...] = 'File'.
% mask:                 Only write files for voxels with 1 in mask; If all should be processed: mask = ones(...) (same size as InArray but no vecSize dim)
%
% CPU_cores:            Determines to how many batch scripts the individual voxels are split, so that those batch scripts can be run in parallel.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux_1_0,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks:





%% 0. Definitions, Preparations, Initialization


% VarNames:             Names of the dimensions of InArray, e.g. {'ROW','COL','SLC','vecSize'}
% mask:                 Only write files for voxels with 1 in mask; If all should be processed: mask = ones(...) (same size as InArray but no vecSize dim)
% out_dir:              Path to where all files should be written
% basis_file:           Path of basis file with which LCModel should fit the spectra
% pat_name:             Name of Patient
% MetaInfo.LarmorFreq:  Frequency where the scanner detected the water-resonance
% MetaInfo.dwelltime:   Time interval between two consecutive FID-points
% zero_order_phase:     Prior Knowledge about the zero order phase of the spectra
% sddegz:               (StandardDeviationDegreeZero: StandardDeviation controlling how much LCModel is allowed to deviate from zero_order_phase
% first_order_phase:    Prior Knowledge about the first order phase of the spectra
% sddegp:               Same as sddegz, but for first_order_phase
% subbas_flag:          flag controlling whether the baseline should be subtracted in the .PS files of LCM
% use_phantom_flag:     If data is from a phantom, some parameters can be different (e.g. phantoms have different metabolites)
% LCM_ProgramPath:      Path of the lcmodel file, e.g. SOME_PATH/.lcmodel/bin/lcmodel






%% 0.1 Stop if not enough info is passed to function
if(nargin < 4)
    display('You must at least input InArray, Paths and MetaInfo.')
    return
end

% if InArray contains WaterReference
if(isstruct(InArray))
	WaterReference = InArray.watref;
	InArray = InArray.csi;
end

% The function can only handle Input Array with 6 dimensions. Quit function if input is larger.
DimNumber_Input = numel(size(InArray));
if(DimNumber_Input > 6)
    display('Dimension of Input Array too large. Rewrite function ''Write_LCM_files'' or shrink Input Array. char(10) Program quitted without output.')
    return
end




%% 0.2 Assign standard values to MetaInfo, mask and CPUcores if nothing is passed to the function.


if(~isfield(MetaInfo,'Dimt1'))
    MetaInfo.Dimt1 = find(max(size(InArray)) == size(InArray));
end
if(~isfield(MetaInfo,'pat_name'))
    MetaInfo.pat_name = 'NoName';
end
if(~isfield(MetaInfo,'DimNames'))
    MetaInfo.DimNames = {'x','y','z'};
end
if(~exist('mask','var'))
    mask = ones(size(squeeze(squeeze_single_dim(InArray,MetaInfo.Dimt1))));
end
if(~exist('CPU_cores','var'))
    CPU_cores = 8;
end
if(exist('ControlInfo','var') && isnumeric(ControlInfo))
	clear ControlInfo;				% Easier to handle this way
end






%% 0.4 DEFINITIONS of Standard Control Info



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                            THESE ARE STANDARD VALUES.                                                               %%%%%
%%%%%   DON'T CHANGE THESE VALUES IF YOU NEED DIFFERENT ADJUSTMENTS, BUT INSTEAD DEFINE THESE VARIABLES IN "ControlParameters" ACCORDING TO YOUR NEEDS.   %%%%%
%%%%%               ONLY CHANGE THE VALUES, IF YOU THINK DIFFERENT STANDARD VALUES ARE APPROPRIATE (BECAUSE THESE ARE TOO SPECIAL E.G.).                  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Referencing, see LCModel Manual page 105
ControlWrite.DOREFS = {'DOREFS(1) = T','DOREFS(2) = F'};           % T: Use standard water referencing, F: No Other metabolites used for referencing
% ControlWrite.NREFPK = {'NREFPK(2) = 1'};                         % NREFPK(JCCF): number of delta functions used for Reference JCCF, e.g. 1 delta function for NAA at 2.01; Each delta function is numbered by JREF
% ControlWrite.PPMREF = {'PPMREF(2,1) = 2.01'};                    % PPMREF(JREF,JCCF): What chemical shift should the delta function with number JREF for the Reference with number JCCF have?
% ControlWrite.SHIFMNMX = {'SHIFMN(2) = -0.1','SHIFMX(2) = 0.1'};  % SHIFMN(JCCF): Minimum chemshift to search for peak. LCModel searches from PPMREF(2,1) + SHIFMN(2) to PPMREF(2,1) + SHIFMX(2)
                                                                  % SHIFMN should be negative. 

                                                                  
% Water Scaling, Absolute Quantification
ControlWrite.WSMET = 'WSMET = ''DSS''';                          % This tells LCModel what to use for scaling the absolute fitting concentrations. 
ControlWrite.WSPPM = 'WSPPM = 0.0';                              % The chemical shift of WSMET   
ControlWrite.N1HMET = 'N1HMET = 9';                              % The number of protons contributing to the signal



% Plotting Parameters
ControlWrite.SUBBAS =  'SUBBAS = T';                             % Subtracts the baseline from the spectra
ControlWrite.NEACH =  'NEACH = 99';                              % "the number of metabolites for which individual plots are to be made." (LCM Manual p. 118)
ControlWrite.WDLINE =  {'WDLINE(6) = 0'};                        % Set the fine grid lines to thickness = 0. 



% Analysis Window
ControlWrite.PPMST = 'PPMST = 4.2';                              % Fit data in chemical shift region [PPMEND, PPMST], PPMST > PPMEND
ControlWrite.PPMEND = 'PPMEND = 1.8';



% Zero Order Phase
ControlWrite.DEGZER =  'DEGZER = 0';                             % zero order phase prior knowledge, set to zero for no prior knowledge
ControlWrite.SDDEGZ =  'SDDEGZ = 999';                           % standard deviation of DEGZER, set to 999 for no prior knowledge.

% First order phase
ControlWrite.DEGPPM =  'DEGPPM = 0';                             % 1st order phase prior knowledge, set to zero for no prior knowledge
ControlWrite.SDDEGP =  'SDDEGP = 20';                            % standard deviation of DEGPPM; note that LCM varies the phase a lot: E.g. for sddegp=1 a total 1.order_phase > 15 is not rare!
%%%%% FIRST ORDER PHASE:
%%%%% if 2 metabolites have circle freq w1, w2 and a phase 0. order of ph0 and they are measured with acquisition delay t
%%%%% then their phases in rad are:
%%%%% phw1(t) = ph0 + w1*t
%%%%% phw2(t) = ph0 + w2*t
%%%%% so their phase difference is delta_ph[rad] = (w2-w1)*t
%%%%% the acq delay is then t = delta_ph[rad]/(w2-w1) in seconds, converted to degree/ppm it is 
%%%%% t [s] = delta_ph[rad]/(w2-w1) = delta_ph[rad]/((f2-f1)*2pi) = delta_ph[rad]/(ppm_difference*297.223*10^6*2pi);  (297.223 @ 7T !)
%%%%% so assuming ppm_difference = 2ppm and a delta_ph[rad]=0.05236rad  equal to 3deg --> 1stOrderPhase = 0.0262 rad/ppm equal to 1.5 deg/ppm
%%%%% this leads to an acq delay of: t ~ 0.014 ms; 
%%%%% so such an deviation (because of wrong basis file or technical inaccuracy) can be compensated





% Basis Set Parameters
ControlWrite.NSIMUL = 'NSIMUL = 0';                              % Don't Simulate additional Basis-spectra that are not in the Basis Set

ControlWrite.NOMIT =  'NOMIT = 5';                               % Number of Metabolites within the Basis Set that should be omitted from the analysis
ControlWrite.CHOMIT =  {'CHOMIT(1) = ''Cho''','CHOMIT(2) = ''Act''','CHOMIT(3) = ''mm3''','CHOMIT(4) = ''mm4''','CHOMIT(5) = ''Glc_B'''}; % Names of omitted metabolites



% ControlWrite.NUSE =  'NUSE = 2'; 
% ControlWrite.CHUSE1 = {'CHUSE(1)=''Act''', 'CHUSE(2)=''Lac'''};% Only Use the following metabolites in the Preliminary Analysis.




% Controls for creating different files
ControlWrite.LTABLE =  'LTABLE = 7';         % Create a .table file
ControlWrite.LCSV =  'LCSV = 0';             % Don't create a .CSV file 
ControlWrite.LCOORD =  'LCOORD = 9';         % Create a Coord file






% OVERWRITE VARIABLES GIVEN BY ControlInfo
if(exist('ControlInfo','var'))
	if(isstruct(ControlInfo))
		struct_fieldnames = fieldnames(ControlInfo);
		struct_fieldnames = struct_fieldnames(cellfun('isempty',regexp(struct_fieldnames,'Path|Priority')));  % Exclude Path and Priority, they should never be given 
		if(isfield(ControlInfo,'Priority') && strcmpi(ControlInfo.Priority,'File'))		% Run first the fields, because then the file can overwrite those values 	
			% Copy existing fields of ControlInfo to ControlWrite struct.
			for field_no = 1:numel(struct_fieldnames); 
				eval(['ControlWrite.' struct_fieldnames{field_no} ' = ControlInfo.' struct_fieldnames{field_no} ';']);
			end
			% Run the file with path ControlParameters if it is a file.
			if(isfield(ControlInfo,'Path') && exist(ControlInfo.Path,'file'))
				eval(['run ' ControlInfo.Path]);	
			end

		else																			% even if Priority is not defined, do it this way as default behavior.
			% Run the file with path ControlParameters if it is a file.
			if(isfield(ControlInfo,'Path') && exist(ControlInfo.Path,'file'))
				eval(['run ' ControlInfo.Path]);	
			end
			% Copy existing fields of ControlInfo to ControlWrite struct.
			for field_no = 1:numel(struct_fieldnames); 
				eval(['ControlWrite.' struct_fieldnames{field_no} ' = ControlInfo.' struct_fieldnames{field_no} ';']);
			end				
		end
		
	else
		if(exist(ControlInfo,'file'))
			% Run the file with path ControlParameters if it is a file.
			eval(['run ' ControlInfo]);	
		end
		
	end
	
end


clear ControlInfo












%% 0.3 Create files, directories, Compute variables


% Create directory
if(~exist(Paths.out_dir,'dir'))
    mkdir(Paths.out_dir);
end
if(~exist([Paths.out_dir '/CoordFiles'],'dir'))
    mkdir([Paths.out_dir '/CoordFiles']);
end

vecSize = size(InArray,MetaInfo.Dimt1);
TotalVoxelNo = sum(reshape(mask,1,[]));


% creating one batch file PER CPU core for LCmodel processing
batch_fids = zeros([1 CPU_cores]);
for core = 1:CPU_cores
    if(exist(sprintf('%s/lcm_process_core_%02d.sh',Paths.batchdir,core),'file'))
        %delete(sprintf('%s/lcm_process_core_%02d.sh',Paths.batchdir,core));
    end
    batch_fids(core) = fopen(sprintf('%s/lcm_process_core_%02d.sh',Paths.batchdir,core),'a');
	fprintf(batch_fids(core), 'echo -e "Starting batch %d, pid = $$, ppid = $PPID"\nsleep 1\n',core);  
end
fprintf(batch_fids(1), 'echo -e ''\nLCModel Processing\t...\t0 %%''\n');  

% Find out the dimensions of InArray over which should be looped (e.g. for [size(InArray) = [32 32 1 1024] these dimensions would be [1 2 3 5 6]
% The dimensions 5 and 6 will be ignored because size(InArray,5) = 1 = size(InArray,6) with the above InArray)
ArrayDimIndices = 1:6;
ArrayDimIndices(MetaInfo.Dimt1) = [];

% Define the string to access the InArray, e.g. 'VarInd1,VarInd2,VarInd3,:'
InArray_AccessStr = cell([1 6]);
InArray_AccessStr(ArrayDimIndices) = {'VarInd1','VarInd2','VarInd3','VarInd4','VarInd5'};
InArray_AccessStr{MetaInfo.Dimt1} = ':';

for CommaLoop = 1:5
    InArray_AccessStr{CommaLoop} = [InArray_AccessStr{CommaLoop} ','];
end
InArray_AccessStr = cell2mat(InArray_AccessStr);
















%% 1. Loop over all entries of the InArray except vecSize entry
no_vox=0;
Percentage_written = false([1 10]); Percentage_written(1) = true;
for VarInd1 = 1:size(InArray,ArrayDimIndices(1))
    for VarInd2 = 1:size(InArray,ArrayDimIndices(2))
        for VarInd3 = 1:size(InArray,ArrayDimIndices(3))
            for VarInd4 = 1:size(InArray,ArrayDimIndices(4))
                for VarInd5 = 1:size(InArray,ArrayDimIndices(5))
            
                    
                    
                    
                    
                    %% 1.1 Loop Preparations
                    %  create files only for voxels inside mask
                    if ( mask(VarInd1,VarInd2,VarInd3,VarInd4,VarInd5) == 0 )
                        continue
                    end
                    single_voxel = eval([ 'squeeze(InArray(' InArray_AccessStr '));' ]);
					


                    % Create the names for control files, e.g. ROW1_COL3 to get the right names for the output_files
                    Filename = MetaInfo.pat_name;

                    VarInd1leadingzero = sprintf('%02d', VarInd1);
                    VarInd2leadingzero = sprintf('%02d', VarInd2);            
                    VarInd3leadingzero = sprintf('%02d', VarInd3);           
                    VarInd4leadingzero = sprintf('%02d', VarInd4);           
                    VarInd5leadingzero = sprintf('%02d', VarInd5);            

                    for DimNamesIndex = 1:numel(MetaInfo.DimNames)

                        Filename = [ Filename '_' MetaInfo.DimNames{DimNamesIndex} num2str(eval([ 'VarInd' num2str(DimNamesIndex) 'leadingzero' ])) ];
                    end


                    
                    
                    
                    
                    %% 1.2 write one RAW file per voxel %%%%%

					
					% Water Reference
					if( exist('WaterReference','var') && (~exist('single_voxel_watref','var') || sum(size(WaterReference)>1) > 1) )	 % If WaterReference does not exist --> immediately jump over
																																	 % If it exists and this is first time it would be called.
						single_voxel_watref = eval([ 'squeeze(WaterReference(' InArray_AccessStr '));' ]);							 % (single_voxel_watref does not exist yet), do it.
																																	 % If it exists, and this is not first time, only do it if
						% Metabo																									 % if WaterReference has more than 1 dimension with size > 1
						voxel_raw_out_watref = sprintf('%s/%s_watref.RAW', Paths.out_dir, Filename);								 % (i.e. not just a vector with 1 dim and singleton dims).
						fid = fopen(voxel_raw_out_watref,'w+');

						%write some header info
						fprintf(fid, ' $SEQPAR\n');
						fprintf(fid, ' hzpppm=%d\n',MetaInfo.LarmorFreq/1000000);
						fprintf(fid, ' $END\n');
						fprintf(fid, ' $NMID\n');
						fprintf(fid, ' fmtdat=''(2e14.5)''\n');
						fprintf(fid, ' $END\n');

						%write data into file
						single_voxel_watref = reshape(single_voxel_watref,[vecSize 1]);
						single_voxel_sep_watref=[real(single_voxel_watref)'; imag(single_voxel_watref)'];
						fprintf(fid, '%14.5e%14.5e\n', single_voxel_sep_watref);

						fclose(fid);
						
					end
					
					
					
					
					% Metabo
                    voxel_raw_out = sprintf('%s/%s.RAW', Paths.out_dir, Filename);
                    fid = fopen(voxel_raw_out,'w+');

                    %write some header info
                    fprintf(fid, ' $SEQPAR\n');
                    fprintf(fid, ' hzpppm=%d\n',MetaInfo.LarmorFreq/1000000);
                    fprintf(fid, ' $END\n');
                    fprintf(fid, ' $NMID\n');
                    fprintf(fid, ' fmtdat=''(2e14.5)''\n');
                    fprintf(fid, ' $END\n');

                    %write data into file
                    single_voxel = reshape(single_voxel,[vecSize 1]);
                    single_voxel_sep=[real(single_voxel)'; imag(single_voxel)'];
                    fprintf(fid, '%14.5e%14.5e\n', single_voxel_sep);

                    fclose(fid);


                    %             % DEBUG MODE
                    %             size(single_voxel)
                    %             size(single_voxel_sep)
                    %             figure; plot(1:2048,real(single_voxel), 'b', 1:2048, imag(single_voxel), 'g'); waitforbuttonpress
                    %             % DEBUG MODE END

                    
                    
                    

                    %% 1.3 write one CONTROL file per voxel

                    % Open Control File
                    voxel_control_out = sprintf('%s/%s.control', Paths.out_dir, Filename);
                    control_fid = fopen(voxel_control_out,'w+');

                    
                    
                    % General Info, Water Frequency, Dwelltime, ...
                    fprintf(control_fid, ' $LCMODL\n');
                    fprintf(control_fid, ' OWNER=''MR Exzellenzzentrum, Radiodiagnostik, MUW''\n');
                    fprintf(control_fid, ' Title=''%s''\n', Filename);
                    fprintf(control_fid, ' HZPPPM=%d, DELTAT=%d, NUNFIL=%i\n',MetaInfo.LarmorFreq/1000000, MetaInfo.dwelltime/1000000000, vecSize);
                    % If for each slice a different basis-file was provided, use all these. Otherwise use only the one and only basis_file{1}
                    if(numel(Paths.basis_file) == size(InArray,ArrayDimIndices(3)))
                        fprintf(control_fid, ' FILBAS=''%s''\n', Paths.basis_file{VarInd3});
                    else
                        fprintf(control_fid, ' FILBAS=''%s''\n', Paths.basis_file{1});                
                    end
                    
                    
                    % Referencing
                    fprintf(control_fid, ' %s\n',ControlWrite.DOREFS{1});
                    fprintf(control_fid, ' %s\n',ControlWrite.DOREFS{2});                   
                    if(isfield(ControlWrite,'NREFPK'))
                        for dumli = 1:numel(ControlWrite.NREFPK)
                            fprintf(control_fid, ' %s\n',ControlWrite.NREFPK{dumli}); 
                        end
                    end
                    if(isfield(ControlWrite,'PPMREF'))
                        for dumli = 1:numel(ControlWrite.PPMREF)                        
                            fprintf(control_fid, ' %s\n',ControlWrite.PPMREF{dumli});   
                        end
                    end
                    if(isfield(ControlWrite,'SHIFMNMX'))
                        for dumli = 1:numel(ControlWrite.SHIFMNMX)                                               
                            fprintf(control_fid, ' %s\n',ControlWrite.SHIFMNMX{dumli});  
                        end    
                    end
                    
                    
                    % Water Scaling, Absolute Quantification
                    
					if(exist('WaterReference','var'))
						fprintf(control_fid, ' DOWS = T\n');								% DO Water Scaling         
						fprintf(control_fid, ' FILH2O=''%s''\n',voxel_raw_out_watref);		% Path of water scaling raw file
					else
	                    fprintf(control_fid, ' %s\n',ControlWrite.WSMET);					% If water referencing is done, this is not so useful. The defaults are ok (use Cr at 3.027 with 3 1H's)    
						fprintf(control_fid, ' %s\n',ControlWrite.WSPPM);                                          
						fprintf(control_fid, ' %s\n',ControlWrite.N1HMET);					
					end
                    
                    
                    % Plotting Parameters
                    fprintf(control_fid, ' %s\n',ControlWrite.SUBBAS);                   
                    fprintf(control_fid, ' %s\n',ControlWrite.NEACH);  
                    for dumli = 1:numel(ControlWrite.WDLINE)
                        fprintf(control_fid, ' %s\n', ControlWrite.WDLINE{dumli});
                    end
                    
                    
                    % Analysis Window
                    fprintf(control_fid, ' %s\n',ControlWrite.PPMST); 
                    fprintf(control_fid, ' %s\n',ControlWrite.PPMEND); 
                    if(isfield(ControlWrite,'PPMGAP'))
                        for dumli = 1:numel(ControlWrite.PPMGAP)
                            fprintf(control_fid, ' %s\n',ControlWrite.PPMGAP{dumli});
                        end
                    end

                    
                    % Zero order phase
                    fprintf(control_fid, ' %s\n',ControlWrite.DEGZER);
                    fprintf(control_fid, ' %s\n',ControlWrite.SDDEGZ);      

                    % First order phase
                    fprintf(control_fid, ' %s\n',ControlWrite.DEGPPM); 
                    fprintf(control_fid, ' %s\n',ControlWrite.SDDEGP);             

                    
                    
                    % Basis Set Parameters
                    fprintf(control_fid, ' %s\n', ControlWrite.NSIMUL);                   
                    fprintf(control_fid, ' %s\n', ControlWrite.NOMIT); 
                    for dumli = 1:numel(ControlWrite.CHOMIT)
                        fprintf(control_fid, ' %s\n', ControlWrite.CHOMIT{dumli});
                    end
                    if(isfield(ControlWrite,'NUSE'))
                        fprintf(control_fid, ' %s\n', ControlWrite.NUSE); 
                        for dumli = 1:numel(ControlWrite.CHUSE1)
                            fprintf(control_fid, ' %s\n', ControlWrite.CHUSE1{dumli});
                        end
                    end
                    if(isfield(ControlWrite,'NKEEP'))
                        fprintf(control_fid, ' %s\n', ControlWrite.NKEEP);
                        for dumli = 1:numel(ControlWrite.CHKEEP)
                            fprintf(control_fid, ' %s\n', ControlWrite.CHKEEP{dumli});
                        end
                    end

                    
                    
                    % Other Parameters. These have to be set in the format ControlWrite.Others[n] (e.g. ControlWrite.Others1), where n is a natural number.
                    fields = fieldnames(ControlWrite);
                    fields(cellfun(@isempty,(regexp(fields,'Others')))) = [];       % Kill all fields that dont have 'Others' in their names.
                    for field_loopy = transpose(fields)
                        fprintf(control_fid, ' %s\n', eval(['ControlWrite.' field_loopy{:}]));
                    end
                    
                    
                    
                    % Controls for creating different files

                    fprintf(control_fid, ' %s\n',ControlWrite.LTABLE);
                    fprintf(control_fid, ' %s\n',ControlWrite.LCSV);
                    fprintf(control_fid, ' %s\n',ControlWrite.LCOORD);   
                    
                    if(strcmp(ControlWrite.LTABLE, 'LTABLE = 7'))
                        fprintf(control_fid, ' FILTAB=''%s/%s.table''\n',Paths.out_dir, Filename);
                    end
                    if(strcmp(ControlWrite.LCSV, 'LCSV = 11'))
                        fprintf(control_fid, ' FILCSV=''%s/%s.CSV''\n',Paths.out_dir, Filename);
                    end
                    if(strcmp(ControlWrite.LCOORD, 'LCOORD = 9'))
                        fprintf(control_fid, ' FILCOO=''%s/CoordFiles/%s.coord''\n', Paths.out_dir, Filename);            
                    end                    
                    fprintf(control_fid, ' FILRAW=''%s/%s.RAW''\n',Paths.out_dir, Filename);
                    fprintf(control_fid, ' FILPS=''%s/%s.ps''\n',Paths.out_dir, Filename);

                    
                    
                    % THE END
                    fprintf(control_fid, ' $END\n');

                    
                    % Close Control File
                    fclose(control_fid);



                    %% 1.4 write commands to batch files - split work to CPU_cores
                    
                    
                    
                    core_tmp = mod(no_vox,CPU_cores) + 1;
                    Percentage_done = round(100*no_vox/TotalVoxelNo);
                    % The if conditions: 1) Print in steps of 10%. 2) If for this percentage (e.g. 20%) the command was already written, dont write it again.
                    % 3) Dont write it for 0%, because a different text should be displayed there (see 0.3). 4) Dont write it for 100%, because this case is treated below.
                    if(mod(Percentage_done,10) == 0 && int8(Percentage_done) ~= 0 && int8(Percentage_done) ~= 100 && ~Percentage_written(int8(Percentage_done/10+1)))
                        Percentage_written(int8(Percentage_done/10+1)) = true;
                        fprintf(batch_fids(core_tmp), 'echo ''\t\t\t...\t%d %%''\n', round(100*no_vox/TotalVoxelNo));
                    end
                    fprintf(batch_fids(core_tmp), ' %s < %s 2>/dev/null\n', Paths.LCM_ProgramPath, voxel_control_out);  % Suppress the LCModel output
                    

                    no_vox = no_vox+1;
            
                    
                    
                    
                end  % VarInd5
            end      % VarInd4
        end          % VarInd3
    end              % VarInd2
end                  % VarInd1


%% Close all batch-fids

fprintf(batch_fids(CPU_cores), 'echo ''\t\t\t...\t100 %%''\n');  
fclose('all');