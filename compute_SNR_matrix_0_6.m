function SNR_mat = compute_SNR_matrix_0_6(csi_mat_time,zerofilling,water_frequency,dwelltime,mask,out_dir,Debug_Mode)
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
CS_vec = compute_chemshift_vector_1_0(water_frequency,dwelltime/10^9,vecSize);
CS_vec_zf = compute_chemshift_vector_1_0(water_frequency,dwelltime/10^9,vecSize*zerofilling);    % dwelltime gets not increased with zero_filling: zeroes get just added at the END of vector


% 0.2 Variables for seeking water

seekwater_chemshift_Downfield = 4.62;		 % initial chemshift region to search for water
seekwater_chemshift_upfield = 4.73;
seekwater_chemshift_stepsize = 0.04;
seekwater_totalsteps = 15;
seekwater_metpeak_height_criterium_value = 55.0;   % max of waterpeak must be waterpeak_border_criterium_value times higher than at the border



% 0.3 Variables for seeking NAA

idealwater_chemshift = 4.7;
idealNAA_chemshift = 2.01;
seekNAA_chemshiftregion_water = 0.2;
seekNAA_chemshiftregion_nowater = 0.2;
seekNAA_SNR_threshold = 45.0;
%seekNAAbasis_upfield_chemshift = 0.1;
%seekNAAbasis_Downfield_chemshift = 0.07;
SeekNAA_BasisRegion_Width_CS = 0.05;
SeekNAA_BasisRegionUpfield_PeakDistance_CS = 0.1;
SeekNAA_BasisRegionDownfield_PeakDistance_CS = 0.1;

SeekNAA_BasisRegion_Width_SP = find(min(abs(CS_vec_zf(1) - SeekNAA_BasisRegion_Width_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekNAA_BasisRegion_Width_CS - CS_vec_zf));
SeekNAA_BasisRegionUpfield_PeakDistance_SP = find(min(abs(CS_vec_zf(1) - SeekNAA_BasisRegionUpfield_PeakDistance_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekNAA_BasisRegionUpfield_PeakDistance_CS - CS_vec_zf));
SeekNAA_BasisRegionDownfield_PeakDistance_SP = ...
find(min(abs(CS_vec_zf(1) - SeekNAA_BasisRegionDownfield_PeakDistance_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekNAA_BasisRegionDownfield_PeakDistance_CS - CS_vec_zf));




% 0.4 Create error log file
if(exist('out_dir','var'))
    error_log_file = sprintf('%s/compute_SNR_error_log.txt',out_dir);
    error_log_fid = fopen(error_log_file, 'w+');
    fprintf(error_log_fid,'Error Log File\nThese errors, sorted by z,y,x of matrix were found:\n\n');
end


% 0.5 SNR_mat

SNR_mat = zeros(ROW,COL,SLC);


% 0.6 Variables for Noise computation

%Noise_chemshift_upfield = CS_vec(1);
%Noise_chemshift_Downfield = CS_vec(1) - 10;
% which SP corresponds to a certain CS? We compute the differences of our CS_vector to this CS; Then we take the absolute value, because we want the minimum "distance" of these differences (distance in 1-dim is abs
% value). Then we create a logical vector, containing zeros and one 1, namely at the position where this minimum is in the vector. "find" finds then the non-zero elements
Noise_chemshift_Downfield = 6; 
Noise_chemshift_upfield = 7.5;
Noise_SP_Downfield = find(min(abs(CS_vec - Noise_chemshift_Downfield)) == abs(CS_vec - Noise_chemshift_Downfield));
Noise_SP_upfield = find(min(abs(CS_vec - Noise_chemshift_upfield)) == abs(CS_vec - Noise_chemshift_upfield));
Noise_polyfit_order = 3; 



% 0.7 out_dir

if(exist('out_dir','var'))
    out_dir_failed = [out_dir, '/failed'];
    out_dir_won = [out_dir, '/won'];
    mkdir(out_dir_failed)
    mkdir(out_dir_won)
end


%% 0.D1 DEBUG MODE: How close are the CS_vec values to the seekwater_chemshift_Downfield and upfield?; is the zerofilling enough?

% diff_seekwaterDownfield_to_chemshiftvec = CS_vec_zf(min(abs(CS_vec_zf - seekwater_chemshift_Downfield)) == abs(CS_vec_zf - seekwater_chemshift_Downfield)) - seekwater_chemshift_Downfield
% diff_seekwaterupfield_to_chemshiftvec = CS_vec_zf(min(abs(CS_vec_zf - seekwater_chemshift_upfield)) == abs(CS_vec_zf - seekwater_chemshift_upfield)) - seekwater_chemshift_upfield
% clear diff_seekwaterDownfield_to_chemshiftvec diff_seekwaterupfield_to_chemshiftvec




%% 0.D2 DEBUG MODE2: Apodize the FID to see what is in the Noise region


% csi_mat_time_apod = transpose(squeeze(csi_mat_time(14,10,SLC,:))) .* exp(-(1:vecSize)*0.01);
% csi_mat_freq_apod = fftshift(fft(csi_mat_time_apod));
% figure
% plot(CS_vec(Noise_SP_upfield:Noise_SP_Downfield),csi_mat_freq_apod(Noise_SP_upfield:Noise_SP_Downfield))
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

           
          display([ 'z = ' num2str(z) ', y = ' num2str(y) ', x = ' num2str(x) ])
            
            
            %% 4. Voxel-specific Preparations
            
            
			if(mask(x,y) == 0)
				continue
			end
			spectrum_zf = transpose(squeeze(csi_mat_freq_zf(x,y,z,:)));
            spectrum_zf_real = real(spectrum_zf);
            spectrum_zf_abs = abs(spectrum_zf);
            spectrum = transpose(squeeze(csi_mat_freq(x,y,z,:)));
            spectrum_real = real(spectrum);
            spectrum_abs = abs(spectrum);
                        



            
            %% 5. Compute the standard deviation of the noise

            Polyfit_Noise = polyfit(Noise_SP_upfield:Noise_SP_Downfield,spectrum_real(Noise_SP_upfield:Noise_SP_Downfield),Noise_polyfit_order);  % In noise region there is baseline; Fit baseline
            Noise_spectrum = spectrum_real(Noise_SP_upfield:Noise_SP_Downfield) - polyval(Polyfit_Noise,Noise_SP_upfield:Noise_SP_Downfield);     % Subtract that baseline
			Noise = std(Noise_spectrum);                                                                                                                                        % Compute std of Noise without baseline
                                                                                                                                                                                

   

			%% 6. Search for water_peak

			% 6.0 Define searching region		
			
			% initialization
			waterpeak_found = 0;
			seekwater_step = 0;
            seek_water_in_real_abs = 'real';

			% while peak not found increase region to look for peak and check again
			while(~waterpeak_found && (seekwater_step <= seekwater_totalsteps-1))

                % convert chemshift to points
				% The exact chemshift is most probable not within the CS_vec, so search for the minimum
				seekwater_SP_Downfield = find(min(abs(CS_vec_zf - seekwater_chemshift_Downfield)) == abs(CS_vec_zf - seekwater_chemshift_Downfield));
				seekwater_SP_upfield = find(min(abs(CS_vec_zf - seekwater_chemshift_upfield)) == abs(CS_vec_zf - seekwater_chemshift_upfield));


				% 6.1 Check for waterpeak in real part: Criteria
                if(strcmpi(seek_water_in_real_abs,'real'))                                                                               % upfield:Downfield because upfield chemshift --> low point
                    [waterpeak_found, waterpeak_SP] = seek_metabolite_0_2(spectrum_zf_real(seekwater_SP_upfield:seekwater_SP_Downfield),Noise,{'metpeak_height'},seekwater_metpeak_height_criterium_value);
                else
                    [waterpeak_found, waterpeak_SP] = seek_metabolite_0_2(spectrum_zf_abs(seekwater_SP_upfield:seekwater_SP_Downfield),Noise,{'metpeak_height'},seekwater_metpeak_height_criterium_value);
                end
                
				waterpeak_SP = waterpeak_SP + seekwater_SP_upfield - 1;     % because to the function only the chemshift region of interest is given. If function gives: peak at SP=1 
                                                                                                 % this is in reality the point described by the formula above

                                                                                  
                                                                                             
                % 6.D1 DEBUG MODE: Plot spectrum where water was sought
                if(Debug_Mode)
                    figure
                    eval([ 'plot(CS_vec_zf(seekwater_SP_upfield:seekwater_SP_Downfield),spectrum_zf_' seek_water_in_real_abs ...   % plot either real or absolute value depending on where water was sought
                           '((seekwater_SP_upfield:seekwater_SP_Downfield)),''k'',''Linewidth'',1.6)' ]);
                    hold on
                    eval([ 'plot(CS_vec_zf(waterpeak_SP),(0:spectrum_zf_' seek_water_in_real_abs '(waterpeak_SP)/500:spectrum_zf_' seek_water_in_real_abs '(waterpeak_SP)),''r'')' ]);
                    hold off
                    set(gca,'XDir','reverse');

                    title(sprintf('Waterpeak %s, Waterpeak found = %d, voxel x_%d y_%d z_%d',seek_water_in_real_abs,waterpeak_found, x,y,z),'Interpreter','none')
                    xlabel('Chemical Shift')
                    ylabel('Signal')
                    legend('Measured Spectrum','Water Peak','Location','Best')
                end
                % DEBUG MODE END              
                                                                                             
			

                
                
				% 6.2: If peak not found: Search in abs value of spectrum, if this has already be done: Enlarge seek-region
                if(~waterpeak_found)
                    if(strcmpi(seek_water_in_real_abs,'real'))
                        seek_water_in_real_abs = 'abs';
                    else
                        seekwater_chemshift_Downfield = seekwater_chemshift_Downfield - seekwater_chemshift_stepsize;
                        seekwater_chemshift_upfield = seekwater_chemshift_upfield + seekwater_chemshift_stepsize;
                        seekwater_step = seekwater_step + 1;
                        seek_water_in_real_abs = 'real';              % start to search again in real part.
                    end
                end

			end %while

            
            warning1_no_waterpeak = ~waterpeak_found;



  
            
            %% 6.D2 DEBUG MODE: Plot spectrum with line where water was found
            
            if(Debug_Mode)
                pause
                figure
                plot(CS_vec_zf,spectrum_zf_real,'k','Linewidth',1.6)
                hold on
                plot(CS_vec_zf(waterpeak_SP),(0:spectrum_zf_real(waterpeak_SP)/500:spectrum_zf_real(waterpeak_SP)),'r','Linewidth',1.6)     % plot line crossing spectrum where water was found.
                hold off
                set(gca,'XDir','reverse');

                title(sprintf('Waterpeak found = %d, voxel x_%d y_%d z_%d', waterpeak_found,x,y,z),'Interpreter','none')
                xlabel('Chemical Shift')
                ylabel('Signal')
                legend('Measured Spectrum','Water Peak','Location','Best') 
            end
            % DEBUG MODE END

		
            
            
			%% 7. Seek & Destroy NAA


			% 7.0 Define Searching Region


			if(waterpeak_found)
				seekNAA_chemshift_Downfield = CS_vec_zf(waterpeak_SP) - (idealwater_chemshift - idealNAA_chemshift) - seekNAA_chemshiftregion_water;
				seekNAA_chemshift_upfield = CS_vec_zf(waterpeak_SP) - (idealwater_chemshift - idealNAA_chemshift) + seekNAA_chemshiftregion_water;
			else
				seekNAA_chemshift_Downfield = idealNAA_chemshift - seekNAA_chemshiftregion_nowater;
				seekNAA_chemshift_upfield = idealNAA_chemshift + seekNAA_chemshiftregion_nowater;
            end
	

			% convert chemshift to point
			seekNAA_SP_Downfield = find(min(abs(CS_vec_zf - seekNAA_chemshift_Downfield)) == abs(CS_vec_zf - seekNAA_chemshift_Downfield));
			seekNAA_SP_upfield = find(min(abs(CS_vec_zf - seekNAA_chemshift_upfield)) == abs(CS_vec_zf - seekNAA_chemshift_upfield));
            seekNAA_SeekPeakRegion = [seekNAA_SP_upfield, seekNAA_SP_Downfield];
            

            
        
            % 7.1 Search for NAA with real part, compute SNR
            
            [NAApeak_found,SNR,NAApeak_SP,NAA_LeftBasis_SP,NAA_RightBasis_SP,LeftBasisRegion,RightBasisRegion] = seek_peak_0_1( ...
             spectrum_zf_real,seekNAA_SeekPeakRegion,Noise,'SNR_Criterion',{[seekNAA_SNRthreshold,SeekNAA_BasisRegion_Width_SP,SeekNAA_BasisRegionUpfield_PeakDistance_SP, SeekNAA_BasisRegionDownfield_PeakDistance_SP]});
            
        
        
        
%             % 7.1 Search for NAA with real part
% 			[NAApeak_found, NAApeak_SP] = seek_metabolite_0_2(spectrum_zf_real(seekNAA_SP_upfield:seekNAA_SP_Downfield), Noise,{'metpeak_height'},seekNAA_metpeak_height_criterium_value);
%             NAApeak_SP = NAApeak_SP + seekNAA_SP_upfield - 1;
%             NAApeak_found_in_real = NAApeak_found;
%             
%             
%             % 7.1D: DEBUG MODE PLOT REGION OF SPECTRUM WHERE TO SEARCH FOR NAA IN REAL PART
%             if(Debug_Mode)
%                 pause
%                 figure
%                 plot(CS_vec_zf(seekNAA_SP_upfield:seekNAA_SP_Downfield),spectrum_zf_real((seekNAA_SP_upfield:seekNAA_SP_Downfield)),'k','Linewidth',1.6)
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
%             
%             
%             
%             % 7.2 Search for NAA with abs value
%             if(~NAApeak_found)
%                 % Try again with absolute value
%                 [NAApeak_found, NAApeak_SP] = seek_metabolite_0_2(spectrum_zf_abs(seekNAA_SP_upfield:seekNAA_SP_Downfield),Noise, {'metpeak_height'},seekNAA_metpeak_height_criterium_value);
%                 NAApeak_SP = NAApeak_SP + seekNAA_SP_upfield - 1; 
%             end
%                      
%                 
%                 
%             % 7.2D: DEBUG MODE: PLOT REGION OF SPECTRUM WHERE TO SEARCH FOR NAA IN ABS VALUE
%             if(Debug_Mode)
%                 pause
%                 figure
%                 plot(CS_vec_zf(seekNAA_SP_upfield:seekNAA_SP_Downfield),spectrum_zf_abs((seekNAA_SP_upfield:seekNAA_SP_Downfield)),'k','Linewidth',1.6)
%                 hold on
%                 plot(CS_vec_zf(NAApeak_SP),(0:spectrum_zf_abs(NAApeak_SP)/500:spectrum_zf_abs(NAApeak_SP)),'r')             % plot line crossing spectrum where water was found.
%                 hold off
%                 set(gca,'XDir','reverse');
% 
%                 title(sprintf('NAApeak abs, NAApeak found = %d, voxel x_%d y_%d z_%d', NAApeak_found,x,y,z),'Interpreter','none')
%                 xlabel('Chemical Shift')
%                 ylabel('Signal')
%                 legend('Measured Spectrum','NAA Peak','Location','Best')
%             end
%             % PLOT MODE END
%             
%                 
%                 
%             error1_no_NAApeak = ~NAApeak_found;
% %                 if(~NAApeak_found)
% %                     error1_no_NAApeak = 1;
% %                     continue
% %                 end
% 
% 
% 
% %pause
% 
% 
% 
% 
% 			%% 8. Search for "Basis" of NAA
% 
%             
% 			[NAA_basis_upfield_found,NAA_basis_Downfield_found, NAA_basis_upfield_SP, NAA_basis_Downfield_SP] = ...
%             seek_basis_of_metabolite_0_1(spectrum_zf_real(seekNAA_SP_upfield:seekNAA_SP_Downfield),CS_vec_zf,'fixed_distance_to_peak', ...
%             [seekNAAbasis_upfield_chemshift,seekNAAbasis_Downfield_chemshift],NAApeak_SP);
%         
%             error2_no_NAAbasis_upfield = ~NAA_basis_upfield_found;
%             error3_no_NAAbasis_Downfield = ~NAA_basis_Downfield_found;
%         
% 
%      
% % 			if(~NAA_basis_upfield_found)
% %                 error2_no_NAAbasis_upfield = 1;                
% % 				continue
% % 			end
% % 			if(~NAA_basis_Downfield_found)
% %                 error3_no_NAAbasis_Downfield = 1;                              
% % 				continue
% % 			end


%			%% 9. Compute the Basis at the chemical shift of the peak; Compute Signal
%            
%            
% 			% Put a line through the points (NAA_basis_upfield_SP, spectrum_zf_real(...)) and (NAA_basis_Downfield_SP,spectrum_zf_real(...)) and compute the value of this line at NAA_peak_SP; 
%             % NOT IMPLEMENTED YET
%             % Now implemented: Average
%                   
% 			NAA_Basis_Signal = mean([spectrum_zf_real(NAA_basis_upfield_SP),spectrum_zf_real(NAA_basis_Downfield_SP)]);
% 			Signal = spectrum_zf_real(NAApeak_SP) - NAA_Basis_Signal;
% 
% 
% 
%             %% 10. Compute SNR
% 
% 			SNR_mat(x,y,z) = Signal/Noise;



            
            
            %% 11. Write error log file
            
            if(exist('out_dir','var') && (warning1_no_waterpeak || error1_no_NAApeak || error2_no_NAAbasis_upfield || error3_no_NAAbasis_Downfield))
                % Write Voxel in error log file
                fprintf(error_log_fid,'z_%d_y_%d_x_%d\n',z,y,x);              
                if(warning1_no_waterpeak)
                	fprintf(error_log_fid,'W1 WARNING: No waterpeak found. Assuming water suppression.\n');
                end
                
                if(error1_no_NAApeak)
                	fprintf(error_log_fid,'E1 ERROR: No NAA-peak found. Move to next voxel.\n');
                    
                    if(error2_no_NAAbasis_upfield) 
                    	fprintf(error_log_fid,'E2 ERROR: Basis of NAA-peak at high chemical shift side not found. Move to next voxel.\n');
                    end
                    
                    if(error3_no_NAAbasis_Downfield) 
                        fprintf(error_log_fid,'E3 ERROR: Basis of NAA-peak at low chemical shift side not found. Move to next voxel.\n');
                    end
                    
                end   
            end	
            
            
            
            

			%% 12. plot: 
            %            - noise region
            %            - whole spectrum
            %            - whole spectrum zerofilling
            %            - seek water with real .OR. absolute value
            %            - seek NAA with real .OR. absolute value
            %            - spectrum with NAA region, NAApeak, Basis1,2;


            % 12.0 Preparations
            
            
            if(NAApeak_found_in_real)
                NAApeak_found_in = 'real';
            else
                NAApeak_found_in = 'abs';
            end           
            
            % Create all figure handles and make them visible or invisible
            if(exist('out_dir','var'))
                noise_fig = figure('visible','off');
                whole_spectrum_fig = figure('visible','off');
                whole_spectrum_zf_fig = figure('visible','off');
                seek_water_fig = figure('visible','off'); 
                seek_NAA_fig = figure('visible','off'); 
                succeeded_spectrum_fig = figure('visible','off'); 
                 
            else
                noise_fig = figure;
                whole_spectrum_fig = figure;
                whole_spectrum_zf_fig = figure;  
                seek_water_fig = figure;           
                seek_NAA_fig = figure;              
                succeeded_spectrum_fig = figure;                   
            end
            
            
            
            
            % 12.1 Noise Region
            set(0,'CurrentFigure',noise_fig)                                                          % set noise_figure as current figure to save data there
			plot(CS_vec(Noise_SP_upfield:Noise_SP_Downfield),Noise_spectrum)
            set(gca,'XDir','reverse');
		    title(sprintf('Noise region, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
            xlabel('Chemical Shift')
            ylabel('Signal')            
            
            
            % 12.2 Whole Spectrum
            set(0,'CurrentFigure',whole_spectrum_fig)                                                                  
            plot(CS_vec,spectrum_real)
            set(gca,'XDir','reverse');
		    title(sprintf('Spectrum No zerofilling, voxel x_%d y_%d z_%d', x,y,z), 'Interpreter', 'none')
            xlabel('Chemical Shift')
            ylabel('Signal')
            
            
            % 12.3 Whole Spectrum zf
            set(0,'CurrentFigure',whole_spectrum_zf_fig)                                                                   
            plot(CS_vec_zf,spectrum_zf_real)
            set(gca,'XDir','reverse');
		    title(sprintf('Spectrum With zerofilling, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
            xlabel('Chemical Shift')
            ylabel('Signal')         
                        
    
            
            % 12.4 Largest region where water was sought in real .OR. absolute value
            set(0,'CurrentFigure',seek_water_fig)              
            eval([ 'plot(CS_vec_zf(seekwater_SP_upfield:seekwater_SP_Downfield),spectrum_zf_' seek_water_in_real_abs ...   % plot either real or absolute value depending on where water was sought
                   '((seekwater_SP_upfield:seekwater_SP_Downfield)),''k'',''Linewidth'',1.6)' ]);
            hold on
            eval([ 'plot(CS_vec_zf(waterpeak_SP),(0:spectrum_zf_' seek_water_in_real_abs '(waterpeak_SP)/500:spectrum_zf_' seek_water_in_real_abs '(waterpeak_SP)),''r'')' ]);            
            set(gca,'XDir','reverse');

            title(sprintf('Waterpeak %s, Waterpeak found = %d,voxel x_%d y_%d z_%d',seek_water_in_real_abs,waterpeak_found, x,y,z),'Interpreter','none')
            xlabel('Chemical Shift')
            ylabel('Signal')
            legend('Measured Spectrum','Water Peak','Location','Best')            


            
            % 12.5 Largest region where NAA was sought in real .OR. absolute value
            set(0,'CurrentFigure',seek_NAA_fig)              
            eval([ 'plot(CS_vec_zf(seekNAA_SP_upfield:seekNAA_SP_Downfield),spectrum_zf_' NAApeak_found_in ...   % plot either real or absolute value depending on where NAA was sought
                   '((seekNAA_SP_upfield:seekNAA_SP_Downfield)),''k'',''Linewidth'',1.6)' ]);
            hold on
            eval([ 'plot(CS_vec_zf(NAApeak_SP),(0:spectrum_zf_' NAApeak_found_in '(NAApeak_SP)/500:spectrum_zf_' NAApeak_found_in '(NAApeak_SP)),''r'')' ]);            
            set(gca,'XDir','reverse');

            title(sprintf('NAApeak %s, NAApeak found = %d,voxel x_%d y_%d z_%d',NAApeak_found_in,NAApeak_found, x,y,z),'Interpreter','none')
            xlabel('Chemical Shift')
            ylabel('Signal')
            legend('Measured Spectrum','NAA Peak','Location','Best')                 
            
            
            
         

            % 12.6 succeeded spectrum
            
            set(0,'CurrentFigure',succeeded_spectrum_fig)              
			plot(CS_vec_zf(NAA_basis_upfield_SP:NAA_basis_Downfield_SP),spectrum_zf_real(NAA_basis_upfield_SP:NAA_basis_Downfield_SP),'k','Linewidth',1.6)
            hold on
            plot(CS_vec_zf(NAApeak_SP),spectrum_zf_real(NAApeak_SP),'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k')
            plot(CS_vec_zf(NAA_basis_upfield_SP),spectrum_zf_real(NAA_basis_upfield_SP),'o','MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','k')
            plot(CS_vec_zf(NAA_basis_Downfield_SP),spectrum_zf_real(NAA_basis_Downfield_SP),'o','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k')
            plot(CS_vec_zf(NAApeak_SP),NAA_Basis_Signal,'o','MarkerSize',6,'MarkerFaceColor','m','MarkerEdgeColor','k')
            hold off
            set(gca,'XDir','reverse');

		    title(sprintf('NAApeak, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
            xlabel('Chemical Shift')
            ylabel('Signal')
            legend('Measured Spectrum','NAA Peak','Basispoint hf NAA','Basispoint lf NAA','Basis NAA','Location','Best')      % hf = high field, lf = low field
            


            
            
            
            
            %% 13. Save all figures

            if(exist('out_dir','var'))
                
                if(NAApeak_found)
                    out_dir_FailedOrWon = out_dir_won;
                else
                    out_dir_FailedOrWon = out_dir_failed;
                end
                
                % 13.1 Noise Region
                saveas(noise_fig,sprintf('%s/x_%s_y_%s_z_%s_Noise_region.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                saveas(noise_fig,sprintf('%s/x_%s_y_%s_z_%s_Noise_region.jpg', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))  
                
                % 13.2 Whole spectrum
                saveas(whole_spectrum_fig,sprintf('%s/x_%s_y_%s_z_%s_whole_spectrum.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                saveas(whole_spectrum_fig,sprintf('%s/x_%s_y_%s_z_%s_whole_spectrum.jpg', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))                 
                
                % 13.3 Whole Spectrum zf
                saveas(whole_spectrum_zf_fig,sprintf('%s/x_%s_y_%s_z_%s_whole_spectrum_zf.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                saveas(whole_spectrum_zf_fig,sprintf('%s/x_%s_y_%s_z_%s_whole_spectrum_zf.jpg', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))                  
                
                % 13.4 Water_peak
                saveas(seek_water_fig,sprintf('%s/x_%s_y_%s_z_%s_seek_water.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                saveas(seek_water_fig,sprintf('%s/x_%s_y_%s_z_%s_seek_water.jpg', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))                  
                
                % 13.5 NAA_peak
                saveas(seek_NAA_fig,sprintf('%s/x_%s_y_%s_z_%s_seek_NAA.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                saveas(seek_NAA_fig,sprintf('%s/x_%s_y_%s_z_%s_seek_NAA.jpg', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))                  
                
                % 13.6 Succeeeeeeded Spectrum
                saveas(succeeded_spectrum_fig,sprintf('%s/x_%s_y_%s_z_%s_NAApeak_and_basis.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                saveas(succeeded_spectrum_fig,sprintf('%s/x_%s_y_%s_z_%s_NAApeak_and_basis.jpg', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                
            end
            
            
            
            %% 14. CLOSE ALL
            
            close all;
            
            
            
        end % x-loop
	end % y-loop
end % z-loop

pause off