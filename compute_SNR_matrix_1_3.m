function [SNR_mat,Shift_mat] = compute_SNR_matrix_1_3(csi_mat_time,zerofilling,water_frequency,dwelltime,mask,out_dir,NoiseRegion_CS, ...
                               SeekNAA_PeakRegionWidth_IfWaterFound,SeekNAA_PeakRegionWidth_IfWaterNotFound,SeekNAA_BasisRegionPeakdistance_CS,Phantom_flag,SeekNAA_InRealOrAbs,quiet)
%
% Compute the SNR for each voxel of a CSI-matrix. The matrix must have the size (ROW,COL,SLC,vecSize), so it must already be summed up over the channels and it must be phased.
% Function will write a error-log-file and plots concerning the NAA_seeking process and the Noise_seeking process if a out_dir is given

%% 0. COVENTIONS, PREPARATIONS, DEFINITIONS

% CONVENTIONS

% CS = Chemical Shift, so the value val_CS is a Chemical Shift
% SP = Spectral Point(s), so the value val_SP is not a Chemical Shift, but the number of points inside the spectrum, counting: left to right (e.g. if you have a vecSize = 1024, then the last point (LoFi) would have 1024 
% SS = Spectral Signal, so the value val_SS is a Signal of the spectrum 
% HiFi = HighField or Upfield, i.e. the Chemical Shift at high field
% LoFi = Low Field or Downfield
% If vectors of CS or SP are given, the first entry is always the HiFi (low SP), the second the LoFi (high SP)
%
% Compute Spectral Point out of Chemical Shift
% Example:  NoiseRegion_SP(2) = find(min(abs(CS_vec - NoiseRegion_CS(2))) == abs(CS_vec - NoiseRegion_CS(2)));
% What is done here? 
% 1) Compute the differences of the CS_vector (= all possible CS-values) and the given CS-value (CS_vec - NoiseRegion_CS(2))
% 2) Find the value inside the CS_vector that is closest to the given CS_value (min(abs(CS_vec - NoiseRegion_CS(2))))
% 3) Create a array where the position of this closest point is given in the CS_vector (min(abs(...)) == abs(...)) by a logical 1 (all other entries of the array are zero)
% 4) Find the position of the 1



% PREPARATIONS


if(~exist('quiet','var'))
    quiet = 'loud';
end
if(~exist('SeekNAA_BasisRegionPeakdistance_CS','var'))
    SeekNAA_BasisRegionPeakdistance_CS = [0.16 0.16];
end
if(~exist('SeekNAA_PeakRegionWidth_IfWaterFound','var'))
    SeekNAA_PeakRegionWidth_IfWaterFound = 0.7;
end
if(~exist('SeekNAA_PeakRegionWidth_IfWaterNotFound','var'))
    SeekNAA_PeakRegionWidth_IfWaterNotFound = 0.7;
end
if(~exist('NoiseRegion_CS','var'))
    NoiseRegion_CS = [6 7.5];
end

pause on




% DEFINITIONS

% 0.1 ROW, COL, SLC, vecSize

ROW = size(csi_mat_time,1);
COL = size(csi_mat_time,2);
SLC = size(csi_mat_time,3);
vecSize = size(csi_mat_time,4);


% 0.1 chemshift vector
CS_vec = compute_chemshift_vector_1_0(water_frequency,dwelltime/10^9,vecSize);                   % Compute the Chemical Shift: For each point of spectrum this function gives the corresponding chemical shift
CS_vec_zf = compute_chemshift_vector_1_0(water_frequency,dwelltime/10^9,vecSize*zerofilling);    % dwelltime gets not increased with zero_filling: zeroes get just added at the END of vector

% 0.2 Variables for seeking water
SeekWater_PeakRegion_LoFi_CS = 4.65-0.2;            % initial LoFi (low chemical shift, low frequency) chemshift region to search for water, CS = ChemicalShift, so this is not measured in spectral points but ppm
SeekWater_PeakRegion_HiFi_CS = 4.65+0.2;            % initial HiFi region
SeekWater_PeakRegion_StepSize_CS = 0.04;                % Stepsize with which the PeakRegion gets increased in each loop
SeekWater_TotalSteps = 1;                               % Total Steps for increasing the PeakRegion    
SeekWater_SNRthreshold = 15.0;                          % Threshold for SNR in order to consider peak as water peak
SeekWater_BasisRegionWidth_CS = 0.08;                   % The Basis values will be searched in a certain region. This value defines how large this region should be in ppm (chemshift)
SeekWater_BasisRegionPeakdistance_HiFi_CS = 0.1;     % Compute the "Basis"-value in the HiFi region compared to peak at a Chemical Shift distance of this value
SeekWater_BasisRegionPeakdistance_LoFi_CS = 0.1;   % Same for LoFi

% Convert CS --> SP
SeekWater_BasisRegionWidth_SP = find(min(abs(CS_vec_zf(1) - SeekWater_BasisRegionWidth_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekWater_BasisRegionWidth_CS - CS_vec_zf));
SeekWater_BasisRegionPeakdistance_HiFi_SP = find(min(abs(CS_vec_zf(1) - SeekWater_BasisRegionPeakdistance_HiFi_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekWater_BasisRegionPeakdistance_HiFi_CS - CS_vec_zf));
SeekWater_BasisRegionPeakdistance_LoFi_SP = find(min(abs(CS_vec_zf(1) - SeekWater_BasisRegionPeakdistance_LoFi_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekWater_BasisRegionPeakdistance_LoFi_CS - CS_vec_zf));



% 0.3 Variables for seeking NAA
IdealWater_CS = 4.65;                                   % Chemical Shift of water in ideal case                
IdealNAA_CS = 2.01;                                     % Chemical Shift of NAA in ideal case
% SeekNAA_PeakRegionWidth_IfWaterFound = 0.7;            % Seek for NAA in the region (IdealNAA_CS - Chemshift of WaterPeak +- this value) If WaterPeak was found
% SeekNAA_PeakRegionWidth_IfWaterNotFound = 0.7;          % Same but if WaterPeak was not found
SeekNAA_SNRthreshold = 2;                               % Threshold for SNR in order to consider peak as NAA peak
SeekNAA_BasisRegionWidth_CS = 0.05;                     % The Basis values will be searched in a certain region. This value defines how large this region should be in ppm (chemshift)
SeekNAA_BasisRegionPeakdistance_CS = sort(SeekNAA_BasisRegionPeakdistance_CS,'descend');
%SeekNAA_BasisRegionPeakdistance_CS(1) = 0.16;       % Compute the "Basis"-value in the HiFi region compared to peak at a Chemical Shift distance of this value (((0.08)))
%SeekNAA_BasisRegionPeakdistance_CS(2) = 0.16;     % Same for LoFi
%SeekNAA_InRealOrAbs = 'abs';

% Convert CS --> SP
SeekNAA_BasisRegionWidth_SP = find(min(abs(CS_vec_zf(1) - SeekNAA_BasisRegionWidth_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekNAA_BasisRegionWidth_CS - CS_vec_zf));
SeekNAA_BasisRegionPeakdistance_HiFi_SP = find(min(abs(CS_vec_zf(1) - SeekNAA_BasisRegionPeakdistance_CS(1) - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekNAA_BasisRegionPeakdistance_CS(1) - CS_vec_zf));
SeekNAA_BasisRegionPeakdistance_LoFi_SP = find(min(abs(CS_vec_zf(1) - SeekNAA_BasisRegionPeakdistance_CS(2) - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekNAA_BasisRegionPeakdistance_CS(2) - CS_vec_zf));



% 0.4 Create error log file
if(exist('out_dir','var'))
    mkdir(out_dir);
    error_log_file = sprintf('%s/compute_SNR_error_log.txt',out_dir);
    error_log_fid = fopen(error_log_file, 'w');
    fprintf(error_log_fid,'Error Log File\nThese errors, sorted by x,y,z of matrix were found:\n\n');
end


% 0.5 initialization

SNR_mat = zeros(ROW,COL,SLC);
Shift_mat = zeros(ROW,COL,SLC);


% 0.6 Variables for Noise computation

% NoiseRegion_CS(2) = 6.0;   % LoFi region where to compute the standard deviation of the noise    6  6
% NoiseRegion_CS(1) = 7.5;   % Same for HiFi                                                     7.5  10.6
% Noise_SmoothValue = 29;
% Noise_ApodValue = 20;

% Compute Specpoints out of Chemical Shift
NoOfNoiseRegions = size(NoiseRegion_CS,2);
NoiseRegion_SP = cell(1,NoOfNoiseRegions);
NoiseRegion_SP_vec = cell(1,NoOfNoiseRegions);
for NoiseIndex = 1:NoOfNoiseRegions
    NoiseRegion_CS{NoiseIndex} = sort(NoiseRegion_CS{NoiseIndex},'descend');
    NoiseRegion_SP{NoiseIndex}(1) = find(min(abs(CS_vec - NoiseRegion_CS{NoiseIndex}(1))) == abs(CS_vec - NoiseRegion_CS{NoiseIndex}(1)));
    NoiseRegion_SP{NoiseIndex}(2) = find(min(abs(CS_vec - NoiseRegion_CS{NoiseIndex}(2))) == abs(CS_vec - NoiseRegion_CS{NoiseIndex}(2)));
    NoiseRegion_SP_vec{NoiseIndex} = NoiseRegion_SP{NoiseIndex}(1) : NoiseRegion_SP{NoiseIndex}(2);
end
NoiseSubplot_Row = round(sqrt(NoOfNoiseRegions));
NoiseSubplot_Col = ceil(sqrt(NoOfNoiseRegions));


% 0.7 Variables for plotting the Metabolite Region
MetaboRegion_HiFi_CS = 4.0;
MetaboRegion_LoFi_CS = 1.0;
% Compute Specpoints out of Chemical Shift
MetaboRegion_SP = [find(min(abs(CS_vec - MetaboRegion_HiFi_CS)) == abs(CS_vec - MetaboRegion_HiFi_CS)), find(min(abs(CS_vec - MetaboRegion_LoFi_CS)) == abs(CS_vec - MetaboRegion_LoFi_CS))];

% Total Voxels to compute & Initialize CurrentVoxel
TotalVoxels = sum(sum(mask));
CurrentVoxel = 1;


% 0.8 out_dir
if(exist('out_dir','var'))
    out_dir_failed = [out_dir, '/failed'];
    out_dir_succeeded = [out_dir, '/succeeded'];
    mkdir(out_dir_failed)
    mkdir(out_dir_succeeded)
end


%% 0.D1 DEBUG MODE: How close are the CS_vec values to the SeekWater_PeakRegion_LoFi_CS and HiFi?; is the zerofilling enough?

% diff_SeekWaterLoFi_to_chemshiftvec = CS_vec_zf(min(abs(CS_vec_zf - SeekWater_PeakRegion_LoFi_CS)) == abs(CS_vec_zf - SeekWater_PeakRegion_LoFi_CS)) - SeekWater_PeakRegion_LoFi_CS
% diff_SeekWaterHiFi_to_chemshiftvec = CS_vec_zf(min(abs(CS_vec_zf - SeekWater_PeakRegion_HiFi_CS)) == abs(CS_vec_zf - SeekWater_PeakRegion_HiFi_CS)) - SeekWater_PeakRegion_HiFi_CS
% clear diff_SeekWaterLoFi_to_chemshiftvec diff_SeekWaterHiFi_to_chemshiftvec




%% 0.D2 DEBUG MODE2: Apodize the FID to see what is in the Noise & Metabolite region


% csi_mat_time_apod = transpose(squeeze(csi_mat_time(22,16,SLC,:))) .* exp(-(1:vecSize)*0.00);
% csi_mat_freq_apod = fftshift(fft(csi_mat_time_apod));
% % csi_mat_freq = fftshift(fft(squeeze(csi_mat_time(14,10,SLC,:))));
% % figure
% % plot(CS_vec(MetaboRegion_SP(1):MetaboRegion_SP(2)),real(csi_mat_freq_apod(MetaboRegion_SP(1):MetaboRegion_SP(2))) ,'r','Linewidth',1.9)
% % hold on
% % plot(CS_vec(MetaboRegion_SP(1):MetaboRegion_SP(2)),real(csi_mat_freq(MetaboRegion_SP(1):MetaboRegion_SP(2))))
% % hold off
% set(gca,'XDir','reverse');
% figure
% plot(CS_vec(NoiseRegion_SP_vec),real(csi_mat_freq_apod(NoiseRegion_SP_vec)))
% set(gca,'XDir','reverse');
% pause


%% 1. zerofilling & Time Fourier Transform

% squeeze nur in bestimmte richtung moeglich?

csi_mat_time_zf = zeros(ROW,COL,SLC,zerofilling*vecSize);
csi_mat_time_zf(:,:,:,1:vecSize) = csi_mat_time;

csi_mat_freq_zf = fftshift(fft(csi_mat_time_zf, [],4),4);
csi_mat_freq = fftshift(fft(csi_mat_time, [],4),4);


%% Debug Mask
% mask = zeros(64);
% mask(28,10) = 1;
% mask(24,14) = 1;
%mask(23,12) = 1;

%% 2. Process each voxel individually: LOOPS

for x = 1:ROW
	for y = 1:COL
        for z = 1:SLC

           
            
            
            %% 4. Voxel-specific Preparations
            
            
            if(mask(x,y) == 0)
				continue
            end
            
            if(~strcmp(quiet,'quiet'))
                display([ char(10) 'Processing Voxel x = ' num2str(x) ', y = ' num2str(y) ', z = ' num2str(z) ', Voxel ' num2str(CurrentVoxel) ' of ' num2str(TotalVoxels) ' . . .'])  
            end
            
            
			spectrum_zf = transpose(squeeze(csi_mat_freq_zf(x,y,z,:)));
            spectrum_zf_real = real(spectrum_zf);
            spectrum_zf_abs = abs(spectrum_zf);
            spectrum = transpose(squeeze(csi_mat_freq(x,y,z,:)));
            spectrum_real = real(spectrum);
            %spectrum_imag = imag(spectrum);           
            %spectrum_abs = abs(spectrum);
                        



            
            %% 5. Compute the standard deviation of the noise

            [Noise,NoiseSignal_Smoothie,NoiseSignal_SubSmooth,Noise_ExcludeSubregions_StartEndPoints] = noise_computation_new_1_3(spectrum,NoiseRegion_SP,'Smooth');   %Apodize 31 Smooth 65            
            
            
			%% 6. Search for water_peak

			% 6.0 Define searching region		
			
			% initialization
			NoWaterPeak_found = 1;
			SeekWater_Step = 0;
            SeekWater_InRealOrAbs = 'abs';

			% while peak not found increase region to look for peak and check again
			while(NoWaterPeak_found && (SeekWater_Step <= SeekWater_TotalSteps-1))

                % convert chemshift to points
				% The exact chemshift is most probable not within the CS_vec, so search for the minimum
				SeekWater_PeakRegion_LoFi_SP = find(min(abs(CS_vec_zf - SeekWater_PeakRegion_LoFi_CS)) == abs(CS_vec_zf - SeekWater_PeakRegion_LoFi_CS));
				SeekWater_PeakRegion_HiFi_SP = find(min(abs(CS_vec_zf - SeekWater_PeakRegion_HiFi_CS)) == abs(CS_vec_zf - SeekWater_PeakRegion_HiFi_CS));
                SeekWater_PeakRegion_SP = [SeekWater_PeakRegion_HiFi_SP, SeekWater_PeakRegion_LoFi_SP];

         
				% 6.1 Check for WaterPeak in real part: Criteria
                if(strcmpi(SeekWater_InRealOrAbs,'real'))                                                                          
                    [NoWaterPeak_found,FoundWaterPeak_SP,FoundWaterPeak_SS,FoundWater_Basis_SP,FoundWater_Basis_SS,FoundWater_MeanBasis_SS,SeekWater_HiFiBasisRegion_SP,SeekWater_LoFiBasisRegion_SP] = ...
                     seek_peak_1_0(spectrum_zf,SeekWater_PeakRegion_SP,'real',0,'SNR_Criterion', ...
                     {{Noise,SeekWater_SNRthreshold,SeekWater_BasisRegionWidth_SP,SeekWater_BasisRegionPeakdistance_HiFi_SP,SeekWater_BasisRegionPeakdistance_LoFi_SP}});
                else
                    [NoWaterPeak_found,FoundWaterPeak_SP,FoundWaterPeak_SS,FoundWater_Basis_SP,FoundWater_Basis_SS,FoundWater_MeanBasis_SS,SeekWater_HiFiBasisRegion_SP,SeekWater_LoFiBasisRegion_SP] = ...
                     seek_peak_1_0(spectrum_zf,SeekWater_PeakRegion_SP,'abs',0,'SNR_Criterion', ...
                     {{Noise,SeekWater_SNRthreshold,SeekWater_BasisRegionWidth_SP,SeekWater_BasisRegionPeakdistance_HiFi_SP,SeekWater_BasisRegionPeakdistance_LoFi_SP}});
                end
 
                FoundWaterPeak_logical = max(FoundWaterPeak_SS{:}) == FoundWaterPeak_SS{:};
                FoundWaterPeak_SP = FoundWaterPeak_SP{FoundWaterPeak_logical};
                %FoundWaterPeak_SS = FoundWaterPeak_SS{FoundWaterPeak_logical};            
                FoundWater_Basis_SP = FoundWater_Basis_SP{FoundWaterPeak_logical};  
                FoundWater_Basis_SS = FoundWater_Basis_SS{FoundWaterPeak_logical};
                FoundWater_MeanBasis_SS = FoundWater_MeanBasis_SS{FoundWaterPeak_logical};    
                SeekWater_HiFiBasisRegion_SP = SeekWater_HiFiBasisRegion_SP{FoundWaterPeak_logical};
                SeekWater_LoFiBasisRegion_SP = SeekWater_LoFiBasisRegion_SP{FoundWaterPeak_logical};
                
                
                % OLD FUNCTION seek_metabolite_0_2
%                 if(strcmpi(SeekWater_InRealOrAbs,'real'))                                                                               % HiFi:LoFi because HiFi chemshift --> low point
%                     [WaterPeak_found, WaterPeak_SP] = seek_metabolite_0_2(spectrum_zf_real(SeekWater_PeakRegion_HiFi_SP:SeekWater_PeakRegion_LoFi_SP),Noise,{'metpeak_height'},SeekWater_SNRthreshold);
%                 else
%                     [WaterPeak_found, WaterPeak_SP] = seek_metabolite_0_2(spectrum_zf_abs(SeekWater_PeakRegion_HiFi_SP:SeekWater_PeakRegion_LoFi_SP),Noise,{'metpeak_height'},SeekWater_SNRthreshold);
%                 end
%                   WaterPeak_SP = WaterPeak_SP + SeekWater_PeakRegion_HiFi_SP - 1;     % because to the function only the chemshift region of interest is given. If function gives: peak at SP=1 
                                                                                                 % this is in reality the point described by the formula above

                                                                                  
                                                                                             
%                 % 6.D1 DEBUG MODE: Plot spectrum where water was sought
%                 figure
%                 eval([ 'plot(CS_vec_zf(SeekWater_PeakRegion_SP(1):SeekWater_PeakRegion_SP(2)),spectrum_zf_' SeekWater_InRealOrAbs ...   % plot either real or absolute value depending on where water was sought
%                        '((SeekWater_PeakRegion_SP(1):SeekWater_PeakRegion_SP(2))),''k'',''Linewidth'',1.6)' ]);
%                 hold on
%                 eval([ 'plot(CS_vec_zf(FoundWaterPeak_SP),(0:spectrum_zf_' SeekWater_InRealOrAbs '(FoundWaterPeak_SP)/500:spectrum_zf_' SeekWater_InRealOrAbs '(FoundWaterPeak_SP)),''r'')' ]);
%                 hold off
%                 set(gca,'XDir','reverse');
% 
%                 title(sprintf('WaterPeak %s, WaterPeak found = %d, voxel x_%d y_%d z_%d',SeekWater_InRealOrAbs,~NoWaterPeak_found, x,y,z),'Interpreter','none')
%                 xlabel('Chemical Shift')
%                 ylabel('Signal')
%                 legend('Measured Spectrum','Water Peak','Location','Best')
%                 % DEBUG MODE END              
                                                                                             
			

                
                
				% 6.2: If peak not found: Set Searching Routine to search in abs value of spectrum, if this has already be done: Enlarge seek-region
                if(NoWaterPeak_found)
                    if(strcmpi(SeekWater_InRealOrAbs,'abs'))
                        SeekWater_InRealOrAbs = 'real';
                    else
                        SeekWater_PeakRegion_LoFi_CS = SeekWater_PeakRegion_LoFi_CS - SeekWater_PeakRegion_StepSize_CS;
                        SeekWater_PeakRegion_HiFi_CS = SeekWater_PeakRegion_HiFi_CS + SeekWater_PeakRegion_StepSize_CS;
                        SeekWater_Step = SeekWater_Step + 1;
                        if(SeekWater_Step <= SeekWater_TotalSteps-1)
                            SeekWater_InRealOrAbs = 'abs';              % start to search again in abs part.
                        end
                    end
                end

			end %while

            
            warning1_no_WaterPeak = NoWaterPeak_found;

      
            
			%% 7. Seek & Destroy NAA


			% 7.0 Define Searching Region


            if(~NoWaterPeak_found)
				SeekNAA_PeakRegion_LoFi_CS = IdealNAA_CS - (IdealWater_CS - CS_vec_zf(FoundWaterPeak_SP)) - SeekNAA_PeakRegionWidth_IfWaterFound/2; % IdealWater_CS - CS_vec_zf(FoundWaterPeak_SP) gives the Shift of the CS
				SeekNAA_PeakRegion_HiFi_CS = IdealNAA_CS - (IdealWater_CS - CS_vec_zf(FoundWaterPeak_SP)) + SeekNAA_PeakRegionWidth_IfWaterFound/2; % so IdealNAA_CS - Shift gives the shifted NAA_peak
			else
				SeekNAA_PeakRegion_LoFi_CS = IdealNAA_CS - SeekNAA_PeakRegionWidth_IfWaterNotFound/2;
				SeekNAA_PeakRegion_HiFi_CS = IdealNAA_CS + SeekNAA_PeakRegionWidth_IfWaterNotFound/2;
            end
	

			% convert chemshift to point
			SeekNAA_PeakRegion_LoFi_SP = find(min(abs(CS_vec_zf - SeekNAA_PeakRegion_LoFi_CS)) == abs(CS_vec_zf - SeekNAA_PeakRegion_LoFi_CS));
			SeekNAA_PeakRegion_HiFi_SP = find(min(abs(CS_vec_zf - SeekNAA_PeakRegion_HiFi_CS)) == abs(CS_vec_zf - SeekNAA_PeakRegion_HiFi_CS));
            SeekNAA_PeakRegion_SP = [SeekNAA_PeakRegion_HiFi_SP, SeekNAA_PeakRegion_LoFi_SP];
            

        
            % 7.1 Search for NAA with real part, compute SNR
            
%             [NAApeak_found,SNR,FoundNAAPeak_SP,FoundNAA_Basis_SP(1),FoundNAA_Basis_SP(2),FoundNAA_Basis_SS(1),FoundNAA_Basis_SS(2),FoundNAA_MeanBasis_SS, ...
%              SeekNAA_HiFiBasisRegion_SP,SeekNAA_LoFiBasisRegion_SP] = seek_peak_0_2( ...
%              spectrum_zf_real,SeekNAA_PeakRegion_SP,Noise,'SNR_Criterion',{[SeekNAA_SNRthreshold,SeekNAA_BasisRegionWidth_SP,SeekNAA_BasisRegionPeakdistance_HiFi_SP, SeekNAA_BasisRegionPeakdistance_LoFi_SP]});

            [NoNAApeak_found,FoundNAAPeak_SP,FoundNAAPeak_SS,FoundNAA_Basis_SP,FoundNAA_Basis_SS,FoundNAA_MeanBasis_SS,SeekNAA_HiFiBasisRegion_SP,SeekNAA_LoFiBasisRegion_SP,SeekNAA_Polyfit,SeekNAA_SubBas] = seek_peak_1_1( ...
             spectrum_zf,SeekNAA_PeakRegion_SP,SeekNAA_InRealOrAbs,1,Phantom_flag,'SNR_Criterion',{{Noise,SeekNAA_SNRthreshold,SeekNAA_BasisRegionWidth_SP,SeekNAA_BasisRegionPeakdistance_HiFi_SP, SeekNAA_BasisRegionPeakdistance_LoFi_SP}});
         
%             [NoNAApeak_found,FoundNAAPeak_SP,FoundNAAPeak_SS,FoundNAA_Basis_SP,FoundNAA_Basis_SS,FoundNAA_MeanBasis_SS,SeekNAA_HiFiBasisRegion_SP,SeekNAA_LoFiBasisRegion_SP,SeekNAA_Polyfit,SeekNAA_SubBas] = seek_peak_0_7( ...
%              spectrum_zf,SeekNAA_PeakRegion_SP,SeekNAA_InRealOrAbs,3,'SubBas',{{Noise,SeekNAA_SNRthreshold,SeekNAA_BasisRegionWidth_SP,SeekNAA_BasisRegionPeakdistance_HiFi_SP, SeekNAA_BasisRegionPeakdistance_LoFi_SP}});     
         
            FoundNAAPeak_logical = max(FoundNAAPeak_SS{:}) == FoundNAAPeak_SS{:};
            FoundNAAPeak_SP = FoundNAAPeak_SP{FoundNAAPeak_logical};
            FoundNAAPeak_SS = FoundNAAPeak_SS{FoundNAAPeak_logical};            
            FoundNAA_Basis_SP = FoundNAA_Basis_SP{FoundNAAPeak_logical};  
            FoundNAA_Basis_SS = FoundNAA_Basis_SS{FoundNAAPeak_logical};
            FoundNAA_MeanBasis_SS = FoundNAA_MeanBasis_SS{FoundNAAPeak_logical};    
            SeekNAA_HiFiBasisRegion_SP = SeekNAA_HiFiBasisRegion_SP{FoundNAAPeak_logical};
            SeekNAA_LoFiBasisRegion_SP = SeekNAA_LoFiBasisRegion_SP{FoundNAAPeak_logical};
                        
         
            NAApeak_found = ~NoNAApeak_found;
            error1_no_NAApeak = NoNAApeak_found;
            NAApeak_found_in = SeekNAA_InRealOrAbs;
        

        


            %% 8. Save SNR & Shift to Matrices
            
            if(~NoNAApeak_found)
                SNR_mat(x,y,z) = (FoundNAAPeak_SS - FoundNAA_MeanBasis_SS)/(2*Noise);
                %SNR_mat(x,y,z) = SNR;
            else
                SNR_mat(x,y,z) = NaN;
            end
            if(~NoWaterPeak_found)
            	Shift_mat(x,y,z) = CS_vec_zf(FoundWaterPeak_SP) - IdealWater_CS;
            else
            	Shift_mat(x,y,z) = NaN; 
            end



         
            
            %% 9. Display Messages 
            
            if(~strcmp(quiet,'quiet'))
                % WaterPeak
                if(~NoWaterPeak_found)
                    display(['Water found at ', num2str(round(CS_vec_zf(FoundWaterPeak_SP)*100)/100), ' ppm. Shift = ', num2str(round((CS_vec_zf(FoundWaterPeak_SP) - IdealWater_CS)*100)/100) ' ppm.'])
                else
                    display('No water found')
                end

                % NAApeak
                if(error1_no_NAApeak)
                    display('Seek ''n destroy NAA   FAILED!')
                else
                    display('Seek ''n destroy NAA   SUCCEEDED!')
                end
            end
            
            
            %% 8. Write error log file
            
%             if(exist('out_dir','var') && (warning1_no_WaterPeak || error1_no_NAApeak || error2_no_NAAbasis_HiFi || error3_no_NAAbasis_LoFi))
            if(exist('out_dir','var') && (warning1_no_WaterPeak || error1_no_NAApeak))  
                % Write Voxel in error log file
                fprintf(error_log_fid,'z_%d_y_%d_x_%d\n',z,y,x);              
                if(warning1_no_WaterPeak)
                	fprintf(error_log_fid,'W1 WARNING: No WaterPeak found. Assuming water suppression.\n');
                end
                
                if(error1_no_NAApeak)
                	fprintf(error_log_fid,'E1 ERROR: No NAA-peak found. Move to next voxel.\n');
                    
%                     if(error2_no_NAAbasis_HiFi) 
%                     	fprintf(error_log_fid,'E2 ERROR: Basis of NAA-peak at high chemical shift side not found. Move to next voxel.\n');
%                     end
%                     
%                     if(error3_no_NAAbasis_LoFi) 
%                         fprintf(error_log_fid,'E3 ERROR: Basis of NAA-peak at low chemical shift side not found. Move to next voxel.\n');
%                     end
                    
                end   
            end	
            
            
            
            

			%% 9. plot: 
            %            - noise region
            %            - Metabolite spectrum
            %            - whole spectrum zerofilling
            %            - seek water with real .OR. absolute value
            %            - seek NAA with real
            %            - spectrum with NAA region, NAApeak, BasisRegions, BasisPoints 1,2;


            % 9.0 Preparations
            
            
%             if(NAApeak_found_in_real)
%                 NAApeak_found_in = 'real';
%             else
%                 NAApeak_found_in = 'abs';
%             end           
            


            %% 9.0.1 Create all figure handles and make them visible or invisible
            if(exist('out_dir','var') && ~strcmpi(out_dir,'No'))
                NoisePeak_fig = figure('visible','off');
                NoiseSubPeak_fig = figure('visible','off');                
                Metabo_spectrum_fig = figure('visible','off');
%               whole_spectrum_zf_fig = figure('visible','off');
                SeekWater_fig = figure('visible','off'); 
                SeekNAA_fig = figure('visible','off'); 
                %succeeded_spectrum_fig = figure('visible','off'); 

                

                %% 9.0.2 Compute the NAA & Water region to plot

                Plot_NAAregion_HiFi_SP = min([SeekNAA_PeakRegion_SP(1), SeekNAA_HiFiBasisRegion_SP(1)]);            % If the PeakRegion is farther left then plot spectrum from PeakRegion
                Plot_NAAregion_LoFi_SP = max([SeekNAA_PeakRegion_SP(2), SeekNAA_LoFiBasisRegion_SP(2)]);        % if Basis Region farther left then Basis Region; Same for right point of spectrum
                Plot_WaterRegion_HiFi_SP = min([SeekWater_PeakRegion_SP(1), SeekWater_HiFiBasisRegion_SP(1)]);      
                Plot_WaterRegion_LoFi_SP = max([SeekWater_PeakRegion_SP(2), SeekWater_LoFiBasisRegion_SP(2)]);




                %% 9.0.3 DEFINE VECTORS TO PLOT

                % Noise
                %Plot_NoiseRegion_x = CS_vec(NoiseRegion_SP_vec);

                % SeekWater
                % SPECTRUM ITSELF
                Plot_SeekWater_Spec_x = CS_vec_zf(Plot_WaterRegion_HiFi_SP:Plot_WaterRegion_LoFi_SP);
                Plot_SeekWater_Spec_y = eval([ 'spectrum_zf_' SeekWater_InRealOrAbs '(Plot_WaterRegion_HiFi_SP:Plot_WaterRegion_LoFi_SP)' ]);       % plot either real or absolute value depending on where water was sought


                % Vertical Lines for Peak & Basis Regions
                Plot_SeekWater_VertLine_Min = min(Plot_SeekWater_Spec_y);
                if(Plot_SeekWater_VertLine_Min >= 0)
                   Plot_SeekWater_VertLine_Min = -max(Plot_SeekWater_Spec_y)/20;    % If the absolute value is plotted, then the vertical lines were so small
                end

                % PEAK REGION VERTICAL LINE
                Plot_SeekWater_PeakRegion_VertLine_Max = eval([ 'max(spectrum_zf_' SeekWater_InRealOrAbs '(SeekWater_PeakRegion_SP(1)),spectrum_zf_' SeekWater_InRealOrAbs '(SeekWater_PeakRegion_SP(2)));' ]);
                Plot_SeekWater_PeakRegion_VertLine_y = [Plot_SeekWater_VertLine_Min,Plot_SeekWater_PeakRegion_VertLine_Max];
                Plot_SeekWater_PeakRegion_VertLine_HiFi_x = repmat(CS_vec_zf(SeekWater_PeakRegion_SP(1)),[1 2]);                % if you use just plot(x,vector), x scalar, only dots get plotted
                Plot_SeekWater_PeakRegion_VertLine_LoFi_x = repmat(CS_vec_zf(SeekWater_PeakRegion_SP(2)),[1 2]);              % So plot([x,x],[y_start,y_end]) 

                % HiFi BASIS REGION VERTICAL LINE
                Plot_SeekWater_HiFiBasisRegion_VertLine_Max = ...
                eval([ 'max([spectrum_zf_' SeekWater_InRealOrAbs '(SeekWater_HiFiBasisRegion_SP(1)),spectrum_zf_' SeekWater_InRealOrAbs '(SeekWater_HiFiBasisRegion_SP(2))]);' ]);
                Plot_SeekWater_HiFiBasisRegion_VertLine_y = [Plot_SeekWater_VertLine_Min,Plot_SeekWater_HiFiBasisRegion_VertLine_Max];
                Plot_SeekWater_HiFiBasisRegion_HiFi_VertLine_x = repmat(CS_vec_zf(SeekWater_HiFiBasisRegion_SP(1)),[1 2]);
                Plot_SeekWater_HiFiBasisRegion_LoFi_VertLine_x = repmat(CS_vec_zf(SeekWater_HiFiBasisRegion_SP(2)), [1 2]);

                % LoFi BASIS REGION VERTICAL LINE
                Plot_SeekWater_LoFiBasisRegion_VertLine_Max = ...
                eval([ 'max([spectrum_zf_' SeekWater_InRealOrAbs '(SeekWater_LoFiBasisRegion_SP(1)),spectrum_zf_' SeekWater_InRealOrAbs '(SeekWater_LoFiBasisRegion_SP(2))]);' ]);
                Plot_SeekWater_LoFiBasisRegion_VertLine_y = [Plot_SeekWater_VertLine_Min,Plot_SeekWater_LoFiBasisRegion_VertLine_Max];
                Plot_SeekWater_LoFiBasisRegion_HiFi_VertLine_x = repmat(CS_vec_zf(SeekWater_LoFiBasisRegion_SP(1)), [1 2]);
                Plot_SeekWater_LoFiBasisRegion_LoFi_VertLine_x = repmat(CS_vec_zf(SeekWater_LoFiBasisRegion_SP(2)), [1 2]); 


                % Found Water Peak Point
                Plot_SeekWater_PeakPoint_y = eval([ 'spectrum_zf_' SeekWater_InRealOrAbs '(FoundWaterPeak_SP);' ]);





                % SeekNAA
                % SPECTRUM ITSELF
                Plot_SeekNAA_Spec_x = CS_vec_zf(Plot_NAAregion_HiFi_SP:Plot_NAAregion_LoFi_SP);
                %Plot_SeekNAA_Spec_y = eval(['spectrum_zf_' SeekNAA_InRealOrAbs '(Plot_NAAregion_HiFi_SP:Plot_NAAregion_LoFi_SP);']);
                Plot_SeekNAA_Spec_y = SeekNAA_SubBas(Plot_NAAregion_HiFi_SP:Plot_NAAregion_LoFi_SP);


                % Vertical Lines for Peak & Basis Regions
                Plot_SeekNAA_VertLine_Min = min(Plot_SeekNAA_Spec_y);
                if(Plot_SeekNAA_VertLine_Min >= 0)
                   Plot_SeekNAA_VertLine_Min = -max(Plot_SeekNAA_Spec_y)/20;
                end            


                % PEAK REGION VERTICAL LINE
                %Plot_SeekNAA_PeakRegion_VertLine_Max = eval(['max(spectrum_zf_' SeekNAA_InRealOrAbs '(SeekNAA_PeakRegion_SP(1)),spectrum_zf_' SeekNAA_InRealOrAbs '(SeekNAA_PeakRegion_SP(2)));']);
                Plot_SeekNAA_PeakRegion_VertLine_Max = max([SeekNAA_SubBas(SeekNAA_PeakRegion_SP(1)), SeekNAA_SubBas(SeekNAA_PeakRegion_SP(2))]);
                Plot_SeekNAA_PeakRegion_VertLine_y = [Plot_SeekNAA_VertLine_Min,Plot_SeekNAA_PeakRegion_VertLine_Max];
                Plot_SeekNAA_PeakRegion_HiFi_VertLine_x = repmat(CS_vec_zf(SeekNAA_PeakRegion_SP(1)),[1 2]);                % if you use just plot(x,vector), x scalar, only dots get plotted
                Plot_SeekNAA_PeakRegion_LoFi_VertLine_x = repmat(CS_vec_zf(SeekNAA_PeakRegion_SP(2)),[1 2]);              % So plot([x,x],[y_start,y_end])    

                % HiFi BASIS REGION VERTICAL LINE
                %Plot_SeekNAA_HiFiBasisRegion_VertLine_Max = eval(['max(spectrum_zf_' SeekNAA_InRealOrAbs '(SeekNAA_HiFiBasisRegion_SP(1)),spectrum_zf_' SeekNAA_InRealOrAbs '(SeekNAA_HiFiBasisRegion_SP(2)));']);
                Plot_SeekNAA_HiFiBasisRegion_VertLine_Max = max([SeekNAA_SubBas(SeekNAA_HiFiBasisRegion_SP(1)), SeekNAA_SubBas(SeekNAA_HiFiBasisRegion_SP(2))]);
                Plot_SeekNAA_HiFiBasisRegion_VertLine_y = [Plot_SeekNAA_VertLine_Min,Plot_SeekNAA_HiFiBasisRegion_VertLine_Max];
                Plot_SeekNAA_HiFiBasisRegion_HiFi_VertLine_x = repmat(CS_vec_zf(SeekNAA_HiFiBasisRegion_SP(1)),[1 2]);
                Plot_SeekNAA_HiFiBasisRegion_LoFi_VertLine_x = repmat(CS_vec_zf(SeekNAA_HiFiBasisRegion_SP(2)), [1 2]);

                % LoFi BASIS REGION VERTICAL LINE
                %Plot_SeekNAA_LoFiBasisRegion_VertLine_Max = eval(['max(spectrum_zf_' SeekNAA_InRealOrAbs '(SeekNAA_LoFiBasisRegion_SP(1)),spectrum_zf_' SeekNAA_InRealOrAbs '(SeekNAA_LoFiBasisRegion_SP(2)));']);
                Plot_SeekNAA_LoFiBasisRegion_VertLine_Max = max([SeekNAA_SubBas(SeekNAA_LoFiBasisRegion_SP(1)), SeekNAA_SubBas(SeekNAA_LoFiBasisRegion_SP(2))]);            
                Plot_SeekNAA_LoFiBasisRegion_VertLine_y = [Plot_SeekNAA_VertLine_Min,Plot_SeekNAA_LoFiBasisRegion_VertLine_Max];
                Plot_SeekNAA_LoFiBasisRegion_HiFi_VertLine_x = repmat(CS_vec_zf(SeekNAA_LoFiBasisRegion_SP(1)), [1 2]);
                Plot_SeekNAA_LoFiBasisRegion_LoFi_VertLine_x = repmat(CS_vec_zf(SeekNAA_LoFiBasisRegion_SP(2)), [1 2]); 











                %% 9.1 Noise & Peak
                set(0,'CurrentFigure',NoisePeak_fig)                                                          % set noise_figure as current figure to save data there
                
                for NoiseIndex = 1:NoOfNoiseRegions
                    subplot(NoiseSubplot_Col,NoiseSubplot_Row,NoiseIndex);
                    plot(CS_vec(NoiseRegion_SP_vec{NoiseIndex}),spectrum_real(NoiseRegion_SP_vec{NoiseIndex}))
                    hold on
                    %plot(Plot_NoiseRegion_x,polyval(Polyfit_Noise_Real,NoiseRegion_SP_vec),'r')
                    plot(CS_vec(NoiseRegion_SP_vec{NoiseIndex}),real(NoiseSignal_Smoothie{NoiseIndex}(NoiseRegion_SP_vec{NoiseIndex})),'r','Linewidth',0.3)   
                    hold off
                    set(gca,'XDir','reverse');
                    if(NoiseIndex ==ceil(NoOfNoiseRegions/2))
                        title(sprintf('Noise Peak subtracted, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
                    end
                    %title(sprintf('Noise and Peak, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
                    xlabel('Chemical Shift')
                    %ylabel('Signal')        
                end



                 %% 9.2 NoiseSubPeak
                set(0,'CurrentFigure',NoiseSubPeak_fig)                                                          % set noise_figure as current figure to save data there
                
                for NoiseIndex = 1:NoOfNoiseRegions

                    subplot(NoiseSubplot_Col,NoiseSubplot_Row,NoiseIndex);
                    plot(CS_vec(NoiseRegion_SP_vec{NoiseIndex}),real(NoiseSignal_SubSmooth{NoiseIndex}(NoiseRegion_SP_vec{NoiseIndex})))
                    hold on
                    % plot X-es in Regions that are excluded
                    for ExcludeSubregions_LoopIndex = 1:size(Noise_ExcludeSubregions_StartEndPoints{NoiseIndex},2)
                        StartPoint = Noise_ExcludeSubregions_StartEndPoints{NoiseIndex}{ExcludeSubregions_LoopIndex}(1);
                        EndPoint = Noise_ExcludeSubregions_StartEndPoints{NoiseIndex}{ExcludeSubregions_LoopIndex}(2);
                        StartCS = CS_vec(StartPoint);
                        EndCS = CS_vec(EndPoint);
                        plot([StartCS,EndCS],[min(real(NoiseSignal_SubSmooth{NoiseIndex}(StartPoint:EndPoint))),max(real(NoiseSignal_SubSmooth{NoiseIndex}(StartPoint:EndPoint)))],'r','Linewidth',1.4) % plot line increasing diagonal from left to right
                        plot([StartCS,EndCS],[max(real(NoiseSignal_SubSmooth{NoiseIndex}(StartPoint:EndPoint))),min(real(NoiseSignal_SubSmooth{NoiseIndex}(StartPoint:EndPoint)))],'r','Linewidth',1.4) % plot line decreasing diagonal from left to right
                    end  
                    hold off
                    set(gca,'XDir','reverse');
                    if(NoiseIndex == floor(NoiseSubplot_Col/2))
                        title(sprintf('Noise Peak subtracted, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
                    end
                        
                    %title(sprintf('Noise Peak subtracted, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
                    xlabel('Chemical Shift')
                    %ylabel('Signal')
                end






                %% 9.3 Metabolite Spectrum
                set(0,'CurrentFigure',Metabo_spectrum_fig)                                                                  
                plot(CS_vec(MetaboRegion_SP(1):MetaboRegion_SP(2)),spectrum_real(MetaboRegion_SP(1):MetaboRegion_SP(2)))
                set(gca,'XDir','reverse');
                title(sprintf('Metabolite Spectrum, voxel x_%d y_%d z_%d', x,y,z), 'Interpreter', 'none')
                xlabel('Chemical Shift')
                ylabel('Signal')






    %             % 9.3 Whole Spectrum zf
    %             set(0,'CurrentFigure',whole_spectrum_zf_fig)                                                                   
    %             plot(CS_vec_zf,spectrum_zf_real)
    %             set(gca,'XDir','reverse');
    % 		      title(sprintf('Spectrum With zerofilling, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
    %             xlabel('Chemical Shift')
    %             ylabel('Signal')         





                %% 9.4 WATER REGION
                set(0,'CurrentFigure',SeekWater_fig) 
                plot(Plot_SeekWater_Spec_x,Plot_SeekWater_Spec_y,'k','Linewidth',1.6)
                hold on

                % Peak & Basis Points            
                WaterLegendHandi = zeros([1 4]);
                plot(CS_vec_zf(FoundWaterPeak_SP),Plot_SeekWater_PeakPoint_y,'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k')                           % Where peak was computed
                plot(CS_vec_zf(FoundWater_Basis_SP(1)),FoundWater_Basis_SS(1),'o','MarkerSize',6,'MarkerFaceColor','c','MarkerEdgeColor','k')         % Where HiFi_Basis was computed 
                plot(CS_vec_zf(FoundWater_Basis_SP(2)),FoundWater_Basis_SS(2),'o','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k')     % Where LoFi_Basis was computed
                WaterLegendHandi(1) = plot(CS_vec_zf(FoundWaterPeak_SP),FoundWater_MeanBasis_SS,'o','MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','k');

                % Peak & Basis Regions
                % HiFi & LoFi SeekWaterPeak_region
                WaterLegendHandi(2) = plot(Plot_SeekWater_PeakRegion_VertLine_HiFi_x,Plot_SeekWater_PeakRegion_VertLine_y, '--r', 'Linewidth',1.6);             % LegendHandi means plot handle, Set handle for Legend
                plot(Plot_SeekWater_PeakRegion_VertLine_LoFi_x,Plot_SeekWater_PeakRegion_VertLine_y,'--r','Linewidth',1.6);
                % HiFi SeekWaterBasis_HiFiRegion 
                WaterLegendHandi(3) = plot(Plot_SeekWater_HiFiBasisRegion_HiFi_VertLine_x,Plot_SeekWater_HiFiBasisRegion_VertLine_y,'--c','Linewidth',1.6);     % Only the LegendHandis get displayed in Legend
                plot(Plot_SeekWater_HiFiBasisRegion_LoFi_VertLine_x,Plot_SeekWater_HiFiBasisRegion_VertLine_y,'--c','Linewidth',1.6); 
                % LoFi SeekWaterBasis_LoFiRegion
                WaterLegendHandi(4) = plot(Plot_SeekWater_LoFiBasisRegion_HiFi_VertLine_x,Plot_SeekWater_LoFiBasisRegion_VertLine_y,'--b','Linewidth',1.6);
                plot(Plot_SeekWater_LoFiBasisRegion_LoFi_VertLine_x,Plot_SeekWater_LoFiBasisRegion_VertLine_y,'--b','Linewidth',1.6); 

                hold off
                set(gca,'XDir','reverse');

                title(sprintf('WaterPeak %s, WaterPeak found = %d,voxel x_%d y_%d z_%d',SeekWater_InRealOrAbs,~NoWaterPeak_found, x,y,z),'Interpreter','none')
                xlabel('Chemical Shift')
                ylabel('Signal')
                leg = legend(WaterLegendHandi, 'Mean Basis Point','Peak Region','HiFi Basis Region','LoFi Basis Region');
                set(leg,'FontSize',6);
                set(leg,'Location','Best');






                %% 9.5 NAA REGION
                set(0,'CurrentFigure',SeekNAA_fig)              
                plot(Plot_SeekNAA_Spec_x,Plot_SeekNAA_Spec_y,'k','Linewidth',1.6);
                hold on

                % Baseline
                %plot(CS_vec_zf(Plot_NAAregion_HiFi_SP:Plot_NAAregion_LoFi_SP),polyval(SeekNAA_Polyfit,Plot_NAAregion_HiFi_SP:Plot_NAAregion_LoFi_SP),'r','Linewidth',1.6)

                % Peak & Basis Points
                NAALegendHandi = zeros([1 4]);
                plot(CS_vec_zf(FoundNAAPeak_SP),FoundNAAPeak_SS,'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k');                               % Where peak was computed
                plot(CS_vec_zf(FoundNAA_Basis_SP(1)),FoundNAA_Basis_SS(1),'o','MarkerSize',6,'MarkerFaceColor','c','MarkerEdgeColor','k');                    % Where HiFi_Basis was computed 
                plot(CS_vec_zf(FoundNAA_Basis_SP(2)),FoundNAA_Basis_SS(2),'o','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k');                         % Where LoFi_Basis was computed
                NAALegendHandi(1) = plot(CS_vec_zf(FoundNAAPeak_SP),FoundNAA_MeanBasis_SS,'o','MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','k');

                % Peak & Basis Regions
                % HiFi & LoFi SeekNAAPeak_region
                NAALegendHandi(2) = plot(Plot_SeekNAA_PeakRegion_HiFi_VertLine_x,Plot_SeekNAA_PeakRegion_VertLine_y, '--r', 'Linewidth',1.6);                                                         
                plot(Plot_SeekNAA_PeakRegion_LoFi_VertLine_x,Plot_SeekNAA_PeakRegion_VertLine_y,'--r','Linewidth',1.6);
                % HiFi SeekNAABasis_HiFiRegion 
                NAALegendHandi(3) = plot(Plot_SeekNAA_HiFiBasisRegion_HiFi_VertLine_x,Plot_SeekNAA_HiFiBasisRegion_VertLine_y,'--c','Linewidth',1.6);                                          
                plot(Plot_SeekNAA_HiFiBasisRegion_LoFi_VertLine_x,Plot_SeekNAA_HiFiBasisRegion_VertLine_y,'--c','Linewidth',1.6); 
                % LoFi SeekNAABasis_LoFiRegion
                NAALegendHandi(4) = plot(Plot_SeekNAA_LoFiBasisRegion_HiFi_VertLine_x,Plot_SeekNAA_LoFiBasisRegion_VertLine_y,'--b','Linewidth',1.6);
                plot(Plot_SeekNAA_LoFiBasisRegion_LoFi_VertLine_x,Plot_SeekNAA_LoFiBasisRegion_VertLine_y,'--b','Linewidth',1.6); 


                hold off
                set(gca,'XDir','reverse');

                title(sprintf('NAApeak %s, NAApeak found = %d,voxel x_%d y_%d z_%d',NAApeak_found_in,NAApeak_found, x,y,z),'Interpreter','none')
                xlabel('Chemical Shift')
                ylabel('Signal')
                leg = legend(NAALegendHandi, 'Mean Basis','Peak Region','HiFi Basis Region','LoFi Basis Region');
                set(leg,'FontSize',6);
                set(leg,'Location','Best');






    %             % 9.6 succeeded spectrum
    % 
    %             set(0,'CurrentFigure',succeeded_spectrum_fig)              
    %             plot(CS_vec_zf(FoundNAA_Basis_SP(1):FoundNAA_Basis_SP(2)),spectrum_zf_real(FoundNAA_Basis_SP(1):FoundNAA_Basis_SP(2)),'k','Linewidth',1.6)
    % %             hold on
    % %             plot(CS_vec_zf(FoundNAAPeak_SP),spectrum_zf_real(FoundNAAPeak_SP),'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k')
    % %             plot(CS_vec_zf(FoundNAA_Basis_SP(1)),spectrum_zf_real(FoundNAA_Basis_SP(1)),'o','MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','k')
    % %             plot(CS_vec_zf(FoundNAA_Basis_SP(2)),spectrum_zf_real(FoundNAA_Basis_SP(2)),'o','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k')
    % %             plot(CS_vec_zf(FoundNAAPeak_SP),FoundNAA_Basis_S,'o','MarkerSize',6,'MarkerFaceColor','m','MarkerEdgeColor','k')
    % %             hold off
    %             set(gca,'XDir','reverse');
    % 
    %             title(sprintf('NAApeak, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
    %             xlabel('Chemical Shift')
    %             ylabel('Signal')
    %             %leg = legend('Measured Spectrum','NAA Peak','Basispoint hf NAA','Basispoint lf NAA','Basis NAA');      % hf = high field, lf = low field
    %             %set(leg,'FontSize',6);
    %             %set(leg,'Location','Best');
    %             






                %% 10. Save all figures

                if(exist('out_dir','var'))

                    if(NAApeak_found)
                        out_dir_FailedOrWon = out_dir_succeeded;
                    else
                        out_dir_FailedOrWon = out_dir_failed;
                    end

                    % 10.1 Noise Region & Peak
                    saveas(NoisePeak_fig,sprintf('%s/x_%s_y_%s_z_%s_NoiseAndPeak', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)),'fig')
                    saveas(NoisePeak_fig,sprintf('%s/x_%s_y_%s_z_%s_NoiseAndPeak', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)),'epsc2')

                    % 10.2 Noise Region & Peak
                    saveas(NoiseSubPeak_fig,sprintf('%s/x_%s_y_%s_z_%s_NoiseSubPeak', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)),'fig')
                    saveas(NoiseSubPeak_fig,sprintf('%s/x_%s_y_%s_z_%s_NoiseSubPeak', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)),'epsc2')                

                    % 10.3 Metabolite spectrum
                    saveas(Metabo_spectrum_fig,sprintf('%s/x_%s_y_%s_z_%s_Metabo_spectrum.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                    saveas(Metabo_spectrum_fig,sprintf('%s/x_%s_y_%s_z_%s_Metabo_spectrum', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)),'epsc2')                 

    %                 % 10.3 Whole Spectrum zf
    %                 saveas(whole_spectrum_zf_fig,sprintf('%s/x_%s_y_%s_z_%s_whole_spectrum_zf.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
    %                 saveas(whole_spectrum_zf_fig,sprintf('%s/x_%s_y_%s_z_%s_whole_spectrum_zf.eps', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))                  

                    % 10.4 Water_peak
                    saveas(SeekWater_fig,sprintf('%s/x_%s_y_%s_z_%s_SeekWater.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                    saveas(SeekWater_fig,sprintf('%s/x_%s_y_%s_z_%s_SeekWater', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)),'epsc2')                  

                    % 10.5 NAA_peak
                    saveas(SeekNAA_fig,sprintf('%s/x_%s_y_%s_z_%s_SeekNAA.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                    saveas(SeekNAA_fig,sprintf('%s/x_%s_y_%s_z_%s_SeekNAA', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)),'epsc2')                  

    %                 % 10.6 Succeeeeeeded Spectrum
    %                 saveas(succeeded_spectrum_fig,sprintf('%s/x_%s_y_%s_z_%s_NAApeak_and_basis.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
    %                 saveas(succeeded_spectrum_fig,sprintf('%s/x_%s_y_%s_z_%s_NAApeak_and_basis.eps', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))

                end
            end
            
            
            

            %% 11. Loop END Preps
            
            close(NoisePeak_fig,NoiseSubPeak_fig, Metabo_spectrum_fig, SeekWater_fig, SeekNAA_fig);
            CurrentVoxel = CurrentVoxel + 1;
            
            
     
        end % x-loop
	end % y-loop
end % z-loop


%% 12. Save SNR map and Shift Map, Close error log file

SNRMap_fig = figure;
imagesc(SNR_mat)
colorbar
saveas(SNRMap_fig,sprintf('%s/SNR_map.fig', out_dir))
saveas(SNRMap_fig,sprintf('%s/SNR_map', out_dir),'epsc2')
close(SNRMap_fig)

ShiftMap_fig = figure;
imagesc(SNR_mat)
colorbar
saveas(ShiftMap_fig,sprintf('%s/Shift_map.fig', out_dir))
saveas(ShiftMap_fig,sprintf('%s/Shift_map', out_dir),'epsc2')
close(ShiftMap_fig)





fclose(error_log_fid);

pause off