function SNR_mat = compute_SNR_matrix_0_8(csi_mat_time,zerofilling,water_frequency,dwelltime,mask,out_dir,Debug_Mode)
%
% Compute the SNR for each voxel of a CSI-matrix. The matrix must have the size (ROW,COL,SLC,vecSize), so it must already be summed up over the channels and it must be phased.
% Function will write a error-log-file and plots concerning the NAA_seeking process and the Noise_seeking process if a out_dir is given

%% 0. PREPARATIONS, DEFINITIONS

% PREPARATIONS

if(~exist('Debug_Mode','var'))
    Debug_Mode = 0;
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

SeekWater_PeakRegion_Downfield_CS = 4.65-0.14;          % initial downfield (low chemical shift, low frequency) chemshift region to search for water, CS = ChemicalShift, so this is not measured in spectral points but ppm
SeekWater_PeakRegion_Upfield_CS = 4.65+0.14;            % initial upfield region
SeekWater_PeakRegion_StepSize_CS = 0.04;                % Stepsize with which the PeakRegion gets increased in each loop
SeekWater_TotalSteps = 1;                               % Total Steps for increasing the PeakRegion    
SeekWater_SNRthreshold = 55.0;                          % Threshold for SNR in order to consider peak as water peak
SeekWater_BasisRegionWidth_CS = 0.1;                    % The Basis values will be searched in a certain region. This value defines how large this region should be in ppm (chemshift)
SeekWater_BasisRegionPeakdistance_Upfield_CS = 0.2;     % Compute the "Basis"-value in the upfield region compared to peak at a Chemical Shift distance of this value
SeekWater_BasisRegionPeakdistance_Downfield_CS = 0.2;   % Same for downfield

% Compute the last three values from Chemical Shift to Specpoints (SP)
% How this is done: We compute the differences of our CS_vector to this CS; Then we take the absolute value, because we want the minimum "distance" of these differences (distance in 1-dim is abs
% value). Then we create a logical vector, containing zeros and one 1, namely at the position where this minimum is in the vector. "find" finds then the non-zero elements
SeekWater_BasisRegionWidth_SP = find(min(abs(CS_vec_zf(1) - SeekWater_BasisRegionWidth_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekWater_BasisRegionWidth_CS - CS_vec_zf));
SeekWater_BasisRegionPeakdistance_Upfield_SP = find(min(abs(CS_vec_zf(1) - SeekWater_BasisRegionPeakdistance_Upfield_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekWater_BasisRegionPeakdistance_Upfield_CS - CS_vec_zf));
SeekWater_BasisRegionPeakdistance_Downfield_SP = find(min(abs(CS_vec_zf(1) - SeekWater_BasisRegionPeakdistance_Downfield_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekWater_BasisRegionPeakdistance_Downfield_CS - CS_vec_zf));



% 0.3 Variables for seeking NAA

IdealWater_CS = 4.65;                                   % Chemical Shift of water in ideal case                
IdealNAA_CS = 2.01;                                     % Chemical Shift of NAA in ideal case
SeekNAA_PeakRegionWidth_IfWaterFound = 0.07;             % Seek for NAA in the region (IdealNAA_CS - Chemshift of WaterPeak +- this value) If WaterPeak was found
SeekNAA_PeakRegionWidth_IfWaterNotFound = 0.2;          % Same but if WaterPeak was not found
seekNAA_SNRthreshold = 4;                               % Threshold for SNR in order to consider peak as NAA peak
SeekNAA_BasisRegionWidth_CS = 0.05;                     % The Basis values will be searched in a certain region. This value defines how large this region should be in ppm (chemshift)
SeekNAA_BasisRegionPeakdistance_Upfield_CS = 0.1;       % Compute the "Basis"-value in the upfield region compared to peak at a Chemical Shift distance of this value
SeekNAA_BasisRegionPeakdistance_Downfield_CS = 0.1;     % Same for downfield

% Compute the last three values from Chemical Shift to Specpoints (SP)
SeekNAA_BasisRegionWidth_SP = find(min(abs(CS_vec_zf(1) - SeekNAA_BasisRegionWidth_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekNAA_BasisRegionWidth_CS - CS_vec_zf));
SeekNAA_BasisRegionPeakdistance_Upfield_SP = find(min(abs(CS_vec_zf(1) - SeekNAA_BasisRegionPeakdistance_Upfield_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekNAA_BasisRegionPeakdistance_Upfield_CS - CS_vec_zf));
SeekNAA_BasisRegionPeakdistance_Downfield_SP = find(min(abs(CS_vec_zf(1) - SeekNAA_BasisRegionPeakdistance_Downfield_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekNAA_BasisRegionPeakdistance_Downfield_CS - CS_vec_zf));



% 0.4 Create error log file
if(exist('out_dir','var'))
    mkdir(out_dir);
    error_log_file = sprintf('%s/compute_SNR_error_log.txt',out_dir);
    error_log_fid = fopen(error_log_file, 'w+');
    fprintf(error_log_fid,'Error Log File\nThese errors, sorted by z,y,x of matrix were found:\n\n');
end


% 0.5 SNR_mat initialization

SNR_mat = zeros(ROW,COL,SLC);


% 0.6 Variables for Noise computation

% NoiseRegion_Upfield_CS = CS_vec(1);
% NoiseRegion_Downfield_CS = CS_vec(1) - 10;
NoiseRegion_Downfield_CS = 6;   % Downfield region where to compute the standard deviation of the noise
NoiseRegion_Upfield_CS = 7.5;   % Same for upfield
Noise_polyfit_order = 2;        % The noise is often not along a flat line, but more or less linear dependent on chemical shift

% Compute Specpoints out of Chemical Shift
NoiseRegion_Downfield_SP = find(min(abs(CS_vec - NoiseRegion_Downfield_CS)) == abs(CS_vec - NoiseRegion_Downfield_CS));
NoiseRegion_Upfield_SP = find(min(abs(CS_vec - NoiseRegion_Upfield_CS)) == abs(CS_vec - NoiseRegion_Upfield_CS));




% 0.7 out_dir

if(exist('out_dir','var'))
    out_dir_failed = [out_dir, '/failed'];
    out_dir_won = [out_dir, '/won'];
    mkdir(out_dir_failed)
    mkdir(out_dir_won)
end


%% 0.D1 DEBUG MODE: How close are the CS_vec values to the SeekWater_PeakRegion_Downfield_CS and upfield?; is the zerofilling enough?

% diff_SeekWaterDownfield_to_chemshiftvec = CS_vec_zf(min(abs(CS_vec_zf - SeekWater_PeakRegion_Downfield_CS)) == abs(CS_vec_zf - SeekWater_PeakRegion_Downfield_CS)) - SeekWater_PeakRegion_Downfield_CS
% diff_SeekWaterupfield_to_chemshiftvec = CS_vec_zf(min(abs(CS_vec_zf - SeekWater_PeakRegion_Upfield_CS)) == abs(CS_vec_zf - SeekWater_PeakRegion_Upfield_CS)) - SeekWater_PeakRegion_Upfield_CS
% clear diff_SeekWaterDownfield_to_chemshiftvec diff_SeekWaterupfield_to_chemshiftvec




%% 0.D2 DEBUG MODE2: Apodize the FID to see what is in the Noise region


% csi_mat_time_apod = transpose(squeeze(csi_mat_time(14,10,SLC,:))) .* exp(-(1:vecSize)*0.0);
% csi_mat_freq_apod = fftshift(fft(csi_mat_time_apod));
% figure
% plot(CS_vec(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP),real(csi_mat_freq_apod(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP)))
% set(gca,'XDir','reverse');
% pause


%% 1. zerofilling & Time Fourier Transform

% squeeze nur in bestimmte richtung möglich?

csi_mat_time_zf = zeros(ROW,COL,SLC,zerofilling*vecSize);
csi_mat_time_zf(:,:,:,1:vecSize) = csi_mat_time;

csi_mat_freq_zf = fftshift(fft(csi_mat_time_zf, [],4),4);
csi_mat_freq = fftshift(fft(csi_mat_time, [],4),4);



%% 2. Process each voxel individually: LOOPS

for z = 1:SLC
	for y = 1:COL
        for x = 1:ROW

           

            
            
            %% 4. Voxel-specific Preparations
            
            
			if(mask(x,y) == 0)
				continue
            end
            
            display([ 'Processing Voxel z = ' num2str(z) ', y = ' num2str(y) ', x = ' num2str(x) ])            
            
            
			spectrum_zf = transpose(squeeze(csi_mat_freq_zf(x,y,z,:)));
            spectrum_zf_real = real(spectrum_zf);
            spectrum_zf_abs = abs(spectrum_zf);
            spectrum = transpose(squeeze(csi_mat_freq(x,y,z,:)));
            spectrum_real = real(spectrum);
            spectrum_imag = imag(spectrum);           
            spectrum_abs = abs(spectrum);
                        



            
            %% 5. Compute the standard deviation of the noise

            % Real Part
            Polyfit_Noise_Real = polyfit(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP,spectrum_real(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP),Noise_polyfit_order);   % In noise region there is baseline; Fit baseline
            Noise_Spectrum_Real = spectrum_real(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP) - polyval(Polyfit_Noise_Real,NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP); % Subtract that baseline
            
            % Imaginary Part
            Polyfit_Noise_Imag = polyfit(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP,spectrum_imag(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP),Noise_polyfit_order);                           
            Noise_Spectrum_Imag = spectrum_imag(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP) - polyval(Polyfit_Noise_Imag,NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP);
            
            % Compute Noise using real and imaginary part, should be the same for both.
            Noise_Spectrum = horzcat(Noise_Spectrum_Real, Noise_Spectrum_Imag);                                                                                                 % Make one long array for std-computation 
			Noise = std(Noise_Spectrum);                                                                                                                                        % Compute std of Noise without baseline
            
                                                                                                                                                                                

%             % 5.D DEBUG MODE: Plot Polyfit and Noise
%             figure
%             plot(CS_vec(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP),spectrum_real(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP))
%             hold on
%             plot(CS_vec(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP),polyval(Polyfit_Noise_Real,NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP),'r')
%             hold off
%             set(gca,'XDir','reverse');
%             
%             figure
%             plot(CS_vec(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP),spectrum_imag(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP))
%             hold on
%             plot(CS_vec(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP),polyval(Polyfit_Noise_Imag,NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP),'r')
%             hold off
%             set(gca,'XDir','reverse');  
%             
%             
%             % plot the real and imaginary Noise - Polyval
%             figure
%             plot(CS_vec(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP),Noise_Spectrum_Real)
%             set(gca,'XDir','reverse');            
%             figure
%             plot(CS_vec(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP),Noise_Spectrum_Imag)
%             set(gca,'XDir','reverse');             
%             
%             pause
   

			%% 6. Search for water_peak

			% 6.0 Define searching region		
			
			% initialization
			WaterPeak_found = 0;
			SeekWater_Step = 0;
            SeekWater_InRealOrAbs = 'real';

			% while peak not found increase region to look for peak and check again
			while(~WaterPeak_found && (SeekWater_Step <= SeekWater_TotalSteps-1))

                % convert chemshift to points
				% The exact chemshift is most probable not within the CS_vec, so search for the minimum
				SeekWater_PeakRegion_Downfield_SP = find(min(abs(CS_vec_zf - SeekWater_PeakRegion_Downfield_CS)) == abs(CS_vec_zf - SeekWater_PeakRegion_Downfield_CS));
				SeekWater_PeakRegion_Upfield_SP = find(min(abs(CS_vec_zf - SeekWater_PeakRegion_Upfield_CS)) == abs(CS_vec_zf - SeekWater_PeakRegion_Upfield_CS));

				% 6.1 Check for WaterPeak in real part: Criteria
                if(strcmpi(SeekWater_InRealOrAbs,'real'))                                                                          
                    [WaterPeak_found,SNR,WaterPeak_SP,Water_UpfieldBasis_SP,Water_DownfieldBasis_SP,Water_Basis_Signal,Water_UpfieldBasisRegion,Water_DownfieldBasisRegion] = ...
                    seek_peak_0_1(spectrum_zf_real,[SeekWater_PeakRegion_Upfield_SP,SeekWater_PeakRegion_Downfield_SP],Noise, ...
                    {'SNR_Criterion'},{[SeekWater_SNRthreshold,SeekWater_BasisRegionWidth_SP,SeekWater_BasisRegionPeakdistance_Upfield_SP,SeekWater_BasisRegionPeakdistance_Downfield_SP]});
                else
                    [WaterPeak_found,SNR,WaterPeak_SP,Water_UpfieldBasis_SP,Water_DownfieldBasis_SP,Water_Basis_Signal,Water_UpfieldBasisRegion,Water_DownfieldBasisRegion] = ...
                    seek_peak_0_1(spectrum_zf_abs,[SeekWater_PeakRegion_Upfield_SP,SeekWater_PeakRegion_Downfield_SP],Noise, ...
                    {'SNR_Criterion'},{[SeekWater_SNRthreshold,SeekWater_BasisRegionWidth_SP,SeekWater_BasisRegionPeakdistance_Upfield_SP,SeekWater_BasisRegionPeakdistance_Downfield_SP]});
                end
 
                
                % OLD FUNCTION seek_metabolite_0_2
%                 if(strcmpi(SeekWater_InRealOrAbs,'real'))                                                                               % upfield:Downfield because upfield chemshift --> low point
%                     [WaterPeak_found, WaterPeak_SP] = seek_metabolite_0_2(spectrum_zf_real(SeekWater_PeakRegion_Upfield_SP:SeekWater_PeakRegion_Downfield_SP),Noise,{'metpeak_height'},SeekWater_SNRthreshold);
%                 else
%                     [WaterPeak_found, WaterPeak_SP] = seek_metabolite_0_2(spectrum_zf_abs(SeekWater_PeakRegion_Upfield_SP:SeekWater_PeakRegion_Downfield_SP),Noise,{'metpeak_height'},SeekWater_SNRthreshold);
%                 end
%                   WaterPeak_SP = WaterPeak_SP + SeekWater_PeakRegion_Upfield_SP - 1;     % because to the function only the chemshift region of interest is given. If function gives: peak at SP=1 
                                                                                                 % this is in reality the point described by the formula above

                                                                                  
                                                                                             
                % 6.D1 DEBUG MODE: Plot spectrum where water was sought
                if(Debug_Mode)
                    figure
                    eval([ 'plot(CS_vec_zf(SeekWater_PeakRegion_Upfield_SP:SeekWater_PeakRegion_Downfield_SP),spectrum_zf_' SeekWater_InRealOrAbs ...   % plot either real or absolute value depending on where water was sought
                           '((SeekWater_PeakRegion_Upfield_SP:SeekWater_PeakRegion_Downfield_SP)),''k'',''Linewidth'',1.6)' ]);
                    hold on
                    eval([ 'plot(CS_vec_zf(WaterPeak_SP),(0:spectrum_zf_' SeekWater_InRealOrAbs '(WaterPeak_SP)/500:spectrum_zf_' SeekWater_InRealOrAbs '(WaterPeak_SP)),''r'')' ]);
                    hold off
                    set(gca,'XDir','reverse');

                    title(sprintf('WaterPeak %s, WaterPeak found = %d, voxel x_%d y_%d z_%d',SeekWater_InRealOrAbs,WaterPeak_found, x,y,z),'Interpreter','none')
                    xlabel('Chemical Shift')
                    ylabel('Signal')
                    legend('Measured Spectrum','Water Peak','Location','Best')
                end
                % DEBUG MODE END              
                                                                                             
			

                
                
				% 6.2: If peak not found: Set Searching Routine to search in abs value of spectrum, if this has already be done: Enlarge seek-region
                if(~WaterPeak_found)
                    if(strcmpi(SeekWater_InRealOrAbs,'real'))
                        SeekWater_InRealOrAbs = 'abs';
                    else
                        SeekWater_PeakRegion_Downfield_CS = SeekWater_PeakRegion_Downfield_CS - SeekWater_PeakRegion_StepSize_CS;
                        SeekWater_PeakRegion_Upfield_CS = SeekWater_PeakRegion_Upfield_CS + SeekWater_PeakRegion_StepSize_CS;
                        SeekWater_Step = SeekWater_Step + 1;
                        SeekWater_InRealOrAbs = 'real';              % start to search again in real part.
                    end
                end

			end %while

            
            warning1_no_WaterPeak = ~WaterPeak_found;



  
            
            %% 6.D2 DEBUG MODE: Plot spectrum with line where water was found
            
            if(Debug_Mode)
                pause
                figure
                plot(CS_vec_zf,spectrum_zf_real,'k','Linewidth',1.6)
                hold on
                plot(CS_vec_zf(WaterPeak_SP),(0:spectrum_zf_real(WaterPeak_SP)/500:spectrum_zf_real(WaterPeak_SP)),'r','Linewidth',1.6)     % plot line crossing spectrum where water was found.
                hold off
                set(gca,'XDir','reverse');

                title(sprintf('WaterPeak found = %d, voxel x_%d y_%d z_%d', WaterPeak_found,x,y,z),'Interpreter','none')
                xlabel('Chemical Shift')
                ylabel('Signal')
                legend('Measured Spectrum','Water Peak','Location','Best') 
            end
            % DEBUG MODE END

		
            
            
			%% 7. Seek & Destroy NAA


			% 7.0 Define Searching Region


			if(WaterPeak_found)
				SeekNAA_PeakRegion_Downfield_CS = CS_vec_zf(WaterPeak_SP) - (IdealWater_CS - IdealNAA_CS) - SeekNAA_PeakRegionWidth_IfWaterFound;
				SeekNAA_PeakRegion_Upfield_CS = CS_vec_zf(WaterPeak_SP) - (IdealWater_CS - IdealNAA_CS) + SeekNAA_PeakRegionWidth_IfWaterFound;
			else
				SeekNAA_PeakRegion_Downfield_CS = IdealNAA_CS - SeekNAA_PeakRegionWidth_IfWaterNotFound;
				SeekNAA_PeakRegion_Upfield_CS = IdealNAA_CS + SeekNAA_PeakRegionWidth_IfWaterNotFound;
            end
	

			% convert chemshift to point
			SeekNAA_PeakRegion_Downfield_SP = find(min(abs(CS_vec_zf - SeekNAA_PeakRegion_Downfield_CS)) == abs(CS_vec_zf - SeekNAA_PeakRegion_Downfield_CS));
			SeekNAA_PeakRegion_Upfield_SP = find(min(abs(CS_vec_zf - SeekNAA_PeakRegion_Upfield_CS)) == abs(CS_vec_zf - SeekNAA_PeakRegion_Upfield_CS));
            seekNAA_SeekPeakRegion = [SeekNAA_PeakRegion_Upfield_SP, SeekNAA_PeakRegion_Downfield_SP];
            

        
            % 7.1 Search for NAA with real part, compute SNR
            
            [NAApeak_found,SNR,NAApeak_SP,NAA_UpfieldBasis_SP,NAA_DownfieldBasis_SP,NAA_Basis_Signal,NAA_UpfieldBasisRegion,NAA_DownfieldBasisRegion] = seek_peak_0_1( ...
             spectrum_zf_real,seekNAA_SeekPeakRegion,Noise,'SNR_Criterion',{[seekNAA_SNRthreshold,SeekNAA_BasisRegionWidth_SP,SeekNAA_BasisRegionPeakdistance_Upfield_SP, SeekNAA_BasisRegionPeakdistance_Downfield_SP]});
            
         
            error1_no_NAApeak = ~NAApeak_found;
            NAApeak_found_in = 'real';
        
        
          
%             % 7.1D: DEBUG MODE PLOT REGION OF SPECTRUM WHERE TO SEARCH FOR NAA IN REAL PART
%             if(Debug_Mode)
%                 pause
%                 figure
%                 plot(CS_vec_zf(SeekNAA_PeakRegion_Upfield_SP:SeekNAA_PeakRegion_Downfield_SP),spectrum_zf_real((SeekNAA_PeakRegion_Upfield_SP:SeekNAA_PeakRegion_Downfield_SP)),'k','Linewidth',1.6)
%                 hold on
%                 plot(CS_vec_zf(NAApeak_SP),(0:spectrum_zf_real(NAApeak_SP)/500:spectrum_zf_real(NAApeak_SP)),'r')             % plot line crossing spectrum where water was found.
%                 hold off
%                 set(gca,'XDir','reverse');
%                 title(sprintf('NAApeak real, NAApeak found = %d, voxel x_%d y_%d z_%d', NAApeak_found,x,y,z),'Interpreter','none')
%                 xlabel('Chemical Shift')
%                 ylabel('Signal')
%                 legend('Measured Spectrum','NAA Peak','Location','Best')  
%             end
%             % DEBUG MODE MODE END





            %% Save SNR to Matrix
            
            SNR_mat(x,y,z) = SNR;


            
            
            %% 8. Write error log file
            
%             if(exist('out_dir','var') && (warning1_no_WaterPeak || error1_no_NAApeak || error2_no_NAAbasis_upfield || error3_no_NAAbasis_Downfield))
            if(exist('out_dir','var') && (warning1_no_WaterPeak || error1_no_NAApeak))  
                % Write Voxel in error log file
                fprintf(error_log_fid,'z_%d_y_%d_x_%d\n',z,y,x);              
                if(warning1_no_WaterPeak)
                	fprintf(error_log_fid,'W1 WARNING: No WaterPeak found. Assuming water suppression.\n');
                end
                
                if(error1_no_NAApeak)
                	fprintf(error_log_fid,'E1 ERROR: No NAA-peak found. Move to next voxel.\n');
                    
%                     if(error2_no_NAAbasis_upfield) 
%                     	fprintf(error_log_fid,'E2 ERROR: Basis of NAA-peak at high chemical shift side not found. Move to next voxel.\n');
%                     end
%                     
%                     if(error3_no_NAAbasis_Downfield) 
%                         fprintf(error_log_fid,'E3 ERROR: Basis of NAA-peak at low chemical shift side not found. Move to next voxel.\n');
%                     end
                    
                end   
            end	
            
            
            
            

			%% 9. plot: 
            %            - noise region
            %            - whole spectrum
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
            


            % 9.0.1 Create all figure handles and make them visible or invisible
            if(exist('out_dir','var'))
                noise_fig = figure('visible','off');
                whole_spectrum_fig = figure('visible','off');
%               whole_spectrum_zf_fig = figure('visible','off');
                seek_water_fig = figure('visible','off'); 
                seek_NAA_fig = figure('visible','off'); 
                succeeded_spectrum_fig = figure('visible','off'); 
            else
                noise_fig = figure;
                whole_spectrum_fig = figure;
%               whole_spectrum_zf_fig = figure;  
                seek_water_fig = figure;           
                seek_NAA_fig = figure;              
                succeeded_spectrum_fig = figure;
            end
                 

            

            % 9.0.2 Compute the NAA & Water region to plot

            Plot_NAAregion_Upfield_SP = min(min(seekNAA_SeekPeakRegion), NAA_UpfieldBasisRegion);
            Plot_NAAregion_Downfield_SP = max(max(seekNAA_SeekPeakRegion), NAA_DownfieldBasisRegion);

            Plot_WaterRegion_Upfield_SP = min(min(seekWater_SeekPeakRegion), Water_UpfieldBasisRegion);
            Plot_WaterRegion_Downfield_SP = max(max(seekWater_SeekPeakRegion), Water_DownfieldBasisRegion);



            % 9.0.3 Define the vectors what to plot
            
            Noise_Region_x = CS_vec(NoiseRegion_Upfield_SP:NoiseRegion_Downfield_SP);
            
            SeekWater_Spec_x = CS_vec_zf(Plot_WaterRegion_Upfield_SP:Plot_WaterRegion_Downfield_SP);
            SeekWater_Spec_y = eval([ 'spectrum_zf_' SeekWater_InRealOrAbs '(Plot_WaterRegion_Upfield_SP:Plot_WaterRegion_Downfield_SP)' ]);
            SeekWater_PeakRegion_Upfield_y = min(SeekWater_Spex_y):spectrum_zf_real(SeekWater_PeakRegion_Upfield_SP)/50:spectrum_zf_real(SeekWater_PeakRegion_Upfield_SP);



            % 9.1 Noise Region
            set(0,'CurrentFigure',noise_fig)                                                          % set noise_figure as current figure to save data there
            plot(Noise_Region_x,Noise_Spectrum_Real)
            set(gca,'XDir','reverse');
            title(sprintf('Noise region, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
            xlabel('Chemical Shift')
            ylabel('Signal')            


            % 9.2 Whole Spectrum
            set(0,'CurrentFigure',whole_spectrum_fig)                                                                  
            plot(CS_vec,spectrum_real)
            set(gca,'XDir','reverse');
            title(sprintf('Spectrum No zerofilling, voxel x_%d y_%d z_%d', x,y,z), 'Interpreter', 'none')
            xlabel('Chemical Shift')
            ylabel('Signal')


%             % 9.3 Whole Spectrum zf
%             set(0,'CurrentFigure',whole_spectrum_zf_fig)                                                                   
%             plot(CS_vec_zf,spectrum_zf_real)
%             set(gca,'XDir','reverse');
% 		      title(sprintf('Spectrum With zerofilling, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
%             xlabel('Chemical Shift')
%             ylabel('Signal')         



            % 9.4 Largest region where water was sought in real .OR. absolute value
            set(0,'CurrentFigure',seek_water_fig)              
            eval([ 'plot(CS_vec_zf(SeekWater_PeakRegion_Upfield_SP:SeekWater_PeakRegion_Downfield_SP),spectrum_zf_' SeekWater_InRealOrAbs ...   % plot either real or absolute value depending on where water was sought
                   '((SeekWater_PeakRegion_Upfield_SP:SeekWater_PeakRegion_Downfield_SP)),''k'',''Linewidth'',1.6)' ]);
            hold on
            % Peak & Basis Points
            plot(CS_vec_zf(WaterPeak_SP),spectrum_zf_real(WaterPeak_SP),'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k')                           % Where peak was computed
            plot(CS_vec_zf(Water_UpfieldBasis_SP),spectrum_zf_real(Water_UpfieldBasis_SP),'o','MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','k')         % Where Upfield_Basis was computed 
            plot(CS_vec_zf(Water_DownfieldBasis_SP),spectrum_zf_real(Water_DownfieldBasis_SP),'o','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k')     % Where Downfield_Basis was computed

            % Peak & Basis Regions
            % Upfield & Downfield SeekWaterPeak_region
            plot(CS_vec_zf(SeekWater_PeakRegion_Upfield_SP),(0:spectrum_zf_real(SeekWater_PeakRegion_Upfield_SP)/500:spectrum_zf_real(SeekWater_PeakRegion_Upfield_SP)),'r', ...
                 CS_vec_zf(SeekWater_PeakRegion_Downfield_SP),(0:spectrum_zf_real(SeekWater_PeakRegion_Downfield_SP)/500:spectrum_zf_real(SeekWater_PeakRegion_Downfield_SP)),'r','Linewidth',2.6); 
            % Upfield SeekWaterBasis_UpfieldRegion  
            plot(CS_vec_zf(Water_UpfieldBasisRegion(1)),(0:spectrum_zf_real(Water_UpfieldBasisRegion(1))/500:spectrum_zf_real(Water_UpfieldBasisRegion(1))), 'g', ...
                 CS_vec_zf(Water_UpfieldBasisRegion(2)),(0:spectrum_zf_real(Water_UpfieldBasisRegion(2))/500:spectrum_zf_real(Water_UpfieldBasisRegion(2))),'g','Linewidth',2.6);          
            % Upfield SeekWaterBasis_DownfieldRegion
            plot(CS_vec_zf(Water_DownfieldBasisRegion(1)),(0:spectrum_zf_real(Water_DownfieldBasisRegion(1))/500:spectrum_zf_real(Water_DownfieldBasisRegion(1))), 'b', ...
                 CS_vec_zf(Water_DownfieldBasisRegion(2)),(0:spectrum_zf_real(Water_DownfieldBasisRegion(2))/500:spectrum_zf_real(Water_DownfieldBasisRegion(2))),'b','Linewidth',2.6); 

            hold off
            set(gca,'XDir','reverse');

            title(sprintf('WaterPeak %s, WaterPeak found = %d,voxel x_%d y_%d z_%d',SeekWater_InRealOrAbs,WaterPeak_found, x,y,z),'Interpreter','none')
            xlabel('Chemical Shift')
            ylabel('Signal')
            leg = legend('Measured Spectrum','Water Peak','','Location','Best');
            set(leg,'FontSize',8);




            % 9.5 Largest region where NAA was sought in real value
            set(0,'CurrentFigure',seek_NAA_fig)              
            plot(CS_vec_zf(Plot_NAAregion_Upfield_SP:Plot_NAAregion_Downfield_SP),spectrum_zf_real(Plot_NAAregion_Upfield_SP:Plot_NAAregion_Downfield_SP),'k','Linewidth',1.6);
            hold on

            % Peak & Basis Points
            plot(CS_vec_zf(NAApeak_SP),spectrum_zf_real(NAApeak_SP),'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k')                           % Where peak was computed
            plot(CS_vec_zf(NAA_UpfieldBasis_SP),spectrum_zf_real(NAA_UpfieldBasis_SP),'o','MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','k')         % Where Upfield_Basis was computed 
            plot(CS_vec_zf(NAA_DownfieldBasis_SP),spectrum_zf_real(NAA_DownfieldBasis_SP),'o','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k')     % Where Downfield_Basis was computed


            % Peak & Basis Regions
            % Upfield & Downfield SeekNAAPeak_region
            plot(CS_vec_zf(SeekNAA_PeakRegion_Upfield_SP),(0:spectrum_zf_real(NAApeak_SP)/500:spectrum_zf_real(NAApeak_SP)),'r', ...
            CS_vec_zf(SeekNAA_PeakRegion_Downfield_SP),(0:spectrum_zf_real(NAApeak_SP)/500:spectrum_zf_real(NAApeak_SP)),'r','Linewidth',2.6); 
            % Upfield SeekNAABasis_UpfieldRegion  
            plot(CS_vec_zf(NAA_UpfieldBasisRegion(1)),(0:spectrum_zf_real(NAApeak_SP)/500:spectrum_zf_real(NAApeak_SP)), 'g', ...
                 CS_vec_zf(NAA_UpfieldBasisRegion(2)),(0:spectrum_zf_real(NAApeak_SP)/500:spectrum_zf_real(NAApeak_SP)),'g','Linewidth',2.6);          
            % Upfield SeekNAABasis_DownfieldRegion
            plot(CS_vec_zf(NAA_DownfieldBasisRegion(1)),(0:spectrum_zf_real(NAApeak_SP)/500:spectrum_zf_real(NAApeak_SP)), 'b', ...
                 CS_vec_zf(NAA_DownfieldBasisRegion(2)),(0:spectrum_zf_real(NAApeak_SP)/500:spectrum_zf_real(NAApeak_SP)),'b','Linewidth',2.6);  

            set(gca,'XDir','reverse');
            hold off

            title(sprintf('NAApeak %s, NAApeak found = %d,voxel x_%d y_%d z_%d',NAApeak_found_in,NAApeak_found, x,y,z),'Interpreter','none')
            xlabel('Chemical Shift')
            ylabel('Signal')
            leg = legend('Measured Spectrum','','','','SeekNAAPeak_region','SeekNAABasis_UpfieldRegion', ...
                   'SeekNAABasis_DownfieldRegion','Location','Best','Interpreter','none');
            set(leg,'FontSize',10);





            % 9.6 succeeded spectrum

            set(0,'CurrentFigure',succeeded_spectrum_fig)              
            plot(CS_vec_zf(NAA_UpfieldBasis_SP:NAA_DownfieldBasis_SP),spectrum_zf_real(NAA_UpfieldBasis_SP:NAA_DownfieldBasis_SP),'k','Linewidth',1.6)
            hold on
            plot(CS_vec_zf(NAApeak_SP),spectrum_zf_real(NAApeak_SP),'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k')
            plot(CS_vec_zf(NAA_UpfieldBasis_SP),spectrum_zf_real(NAA_UpfieldBasis_SP),'o','MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','k')
            plot(CS_vec_zf(NAA_DownfieldBasis_SP),spectrum_zf_real(NAA_DownfieldBasis_SP),'o','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k')
            plot(CS_vec_zf(NAApeak_SP),NAA_Basis_Signal,'o','MarkerSize',6,'MarkerFaceColor','m','MarkerEdgeColor','k')
            hold off
            set(gca,'XDir','reverse');

            title(sprintf('NAApeak, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
            xlabel('Chemical Shift')
            ylabel('Signal')
            legend('Measured Spectrum','NAA Peak','Basispoint hf NAA','Basispoint lf NAA','Basis NAA','Location','Best')      % hf = high field, lf = low field
            


                  
  
            
            
            %% 10. Save all figures

            if(exist('out_dir','var'))
                
                if(NAApeak_found)
                    out_dir_FailedOrWon = out_dir_won;
                else
                    out_dir_FailedOrWon = out_dir_failed;
                end
                
                % 10.1 Noise Region
                saveas(noise_fig,sprintf('%s/x_%s_y_%s_z_%s_Noise_region.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                saveas(noise_fig,sprintf('%s/x_%s_y_%s_z_%s_Noise_region.jpg', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))  
                
                % 10.2 Whole spectrum
                saveas(whole_spectrum_fig,sprintf('%s/x_%s_y_%s_z_%s_whole_spectrum.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                saveas(whole_spectrum_fig,sprintf('%s/x_%s_y_%s_z_%s_whole_spectrum.jpg', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))                 
                
%                 % 10.3 Whole Spectrum zf
%                 saveas(whole_spectrum_zf_fig,sprintf('%s/x_%s_y_%s_z_%s_whole_spectrum_zf.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
%                 saveas(whole_spectrum_zf_fig,sprintf('%s/x_%s_y_%s_z_%s_whole_spectrum_zf.jpg', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))                  
                
                % 10.4 Water_peak
                saveas(seek_water_fig,sprintf('%s/x_%s_y_%s_z_%s_seek_water.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                saveas(seek_water_fig,sprintf('%s/x_%s_y_%s_z_%s_seek_water.jpg', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))                  
                
                % 10.5 NAA_peak
                saveas(seek_NAA_fig,sprintf('%s/x_%s_y_%s_z_%s_seek_NAA.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                saveas(seek_NAA_fig,sprintf('%s/x_%s_y_%s_z_%s_seek_NAA.jpg', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))                  
                
                % 10.6 Succeeeeeeded Spectrum
                saveas(succeeded_spectrum_fig,sprintf('%s/x_%s_y_%s_z_%s_NAApeak_and_basis.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                saveas(succeeded_spectrum_fig,sprintf('%s/x_%s_y_%s_z_%s_NAApeak_and_basis.jpg', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                
            end
            
            
            
            %% 13. CLOSE ALL
            
            close all;
            
            
            
        end % x-loop
	end % y-loop
end % z-loop

pause off