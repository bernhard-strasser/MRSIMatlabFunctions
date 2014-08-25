function SNR_mat = compute_SNR_matrix_0_4(CSI_mat_time,zerofilling,water_frequency,dwelltime,mask,out_dir)
%
% Compute the SNR for each voxel of a CSI-matrix. The matrix must have the size (ROW,COL,SLC,vecSize), so it must already be summed up over the channels
% Function will write a error-log-file and plots concerning the NAA_seeking process and the Noise_seeking process if a out_dir is given

%% 0. PREPARATIONS, DEFINITIONS

pause on


% 0.1 ROW, COL, SLC, vecSize

ROW = size(CSI_mat_time,1);
COL = size(CSI_mat_time,2);
SLC = size(CSI_mat_time,3);
vecSize = size(CSI_mat_time,4);


% 0.1 chemshift vector
chemshift_vector = compute_chemshift_vector_1_0(water_frequency,dwelltime/10^9,vecSize);
chemshift_vector_zf = compute_chemshift_vector_1_0(water_frequency,dwelltime/10^9,vecSize*zerofilling);    % dwelltime gets not increased with zero_filling: zeroes get just added at the END of vector


% 0.2 Variables for seeking water

seekwater_chemshift_lowfield = 4.62;		 % initial chemshift region to search for water
seekwater_chemshift_upfield = 4.73;
seekwater_chemshift_stepsize = 0.04;
seekwater_totalsteps = 15;
seekwater_metpeak_height_criterium_value = 55.0;   % max of waterpeak must be waterpeak_border_criterium_value times higher than at the border



% 0.3 Variables for seeking NAA

idealwater_chemshift = 4.7;
idealNAA_chemshift = 2.01;
seekNAA_chemshiftregion_water = 0.2;
seekNAA_chemshiftregion_nowater = 0.2;
seekNAA_metpeak_height_criterium_value = 45.0;
seekNAAbasis_upfield_chemshift = 0.1;
seekNAAbasis_lowfield_chemshift = 0.07;



% 0.4 Create error log file
if(exist('out_dir','var'))
    mkdir(sprintf('%s/SNR_computations/',out_dir));
    error_log_file = sprintf('%s/SNR_computations/compute_SNR_error_log.txt',out_dir);
    error_log_fid = fopen(error_log_file, 'w+');
    fprintf(error_log_fid,'Error Log File\nThese errors, sorted by z,y,x of matrix were found:\n\n');
end


% 0.5 SNR_mat

SNR_mat = zeros(ROW,COL,SLC);


% 0.6 Variables for Noise computation

%Noise_chemshift_upfield = chemshift_vector(1);
%Noise_chemshift_lowfield = chemshift_vector(1) - 10;
Noise_chemshift_upfield = 24;
Noise_chemshift_lowfield = 8; 
Noise_specpoint_lowfield = find(min(abs(chemshift_vector - Noise_chemshift_lowfield)) == abs(chemshift_vector - Noise_chemshift_lowfield));
Noise_specpoint_upfield = find(min(abs(chemshift_vector - Noise_chemshift_upfield)) == abs(chemshift_vector - Noise_chemshift_upfield));
Noise_polyfit_order = 3; 



%% 0.D DEBUG MODE: How close are the chemshift_vector values to the seekwater_chemshift_lowfield and upfield?; is the zerofilling enough?

% diff_seekwaterlowfield_to_chemshiftvec = chemshift_vector_zf(min(abs(chemshift_vector_zf - seekwater_chemshift_lowfield)) == abs(chemshift_vector_zf - seekwater_chemshift_lowfield)) - seekwater_chemshift_lowfield
% diff_seekwaterupfield_to_chemshiftvec = chemshift_vector_zf(min(abs(chemshift_vector_zf - seekwater_chemshift_upfield)) == abs(chemshift_vector_zf - seekwater_chemshift_upfield)) - seekwater_chemshift_upfield
% clear diff_seekwaterlowfield_to_chemshiftvec diff_seekwaterupfield_to_chemshiftvec



%% 1. OVERSAMPLING & Time Fourier Transform

% squeeze nur in bestimmte richtung möglich?

CSI_mat_time_zf = zeros(ROW,COL,SLC,zerofilling*vecSize);
CSI_mat_time_zf(:,:,:,1:vecSize) = CSI_mat_time;

CSI_mat_freq_zf = fftshift(fft(CSI_mat_time_zf, [],4),4);
CSI_mat_freq = fftshift(fft(CSI_mat_time, [],4),4);



%% 2. Process each voxel individually: LOOPS

for z = 1:SLC
	for y = 1:COL
		for x = 1:ROW

           
            
            
            
            %% 4. Voxel-specific Preparations
            
            
			if(mask(x,y) == 0)
				continue
			end
			spectrum_zf = transpose(squeeze(CSI_mat_freq_zf(x,y,z,:)));
            spectrum_zf_real = real(spectrum_zf);
            spectrum_zf_abs = abs(spectrum_zf);
            spectrum = transpose(squeeze(CSI_mat_freq(x,y,z,:)));
            spectrum_real = real(spectrum);
            spectrum_abs = abs(spectrum);
            
            if(exist('out_dir','var'))
                % Write Voxel in error log file
                fprintf(error_log_fid,'z_%d_y_%d_x_%d\n',z,y,x);
            end


            
            
            %% 4.D DEBUG MODE: PLOT SPECTRUM
            
            figure
            plot(chemshift_vector,abs(spectrum))
            set(gca,'XDir','reverse');
		    title(sprintf('Spectrum No OverSampling, voxel x_%d y_%d z_%d', x,y,z), 'Interpreter', 'none')
            xlabel('Chemical Shift')
            ylabel('Signal')
            
            figure
            plot(chemshift_vector_zf,spectrum_zf_abs)
            set(gca,'XDir','reverse');
		    title(sprintf('Spectrum With OverSampling, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
            xlabel('Chemical Shift')
            ylabel('Signal')         
            

            %pause
            close all;
         

            
            
            %% 5. Compute the standard deviation of the noise

            Polyfit_Noise = polyfit(Noise_specpoint_upfield:Noise_specpoint_lowfield,spectrum_real(Noise_specpoint_upfield:Noise_specpoint_lowfield),Noise_polyfit_order);  % In noise region there is baseline; Fit baseline
            Noise_spectrum = spectrum_real(Noise_specpoint_upfield:Noise_specpoint_lowfield) - polyval(Polyfit_Noise,Noise_specpoint_upfield:Noise_specpoint_lowfield);     % Subtract that baseline
			Noise = std(Noise_spectrum);                                                                                                                                        % Compute std of Noise without baseline
                                                                                                                                                                                
            
            
            
            %% 5.D DEBUG MODE: PLOT NOISE
            
            figure
			plot(chemshift_vector(Noise_specpoint_upfield:Noise_specpoint_lowfield),Noise_spectrum)
            set(gca,'XDir','reverse');
		    title(sprintf('Noise region, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
            xlabel('Chemical Shift')
            ylabel('Signal')
            %pause
            close(gcf);
            
            
            

			%% 6. Search for water_peak

			% 6.0 Define searching region		
			
			% initialization
			waterpeak_found = 0;
			seekwater_step = 0;

			% while peak not found increase region to look for peak and check again
			while(~waterpeak_found && (seekwater_step <= seekwater_totalsteps-1))

                % convert chemshift to points
				% The exact chemshift is most probable not within the chemshift_vector, so search for the minimum
				seekwater_specpoint_lowfield = find(min(abs(chemshift_vector_zf - seekwater_chemshift_lowfield)) == abs(chemshift_vector_zf - seekwater_chemshift_lowfield));
				seekwater_specpoint_upfield = find(min(abs(chemshift_vector_zf - seekwater_chemshift_upfield)) == abs(chemshift_vector_zf - seekwater_chemshift_upfield));
				%[seekwater_specpoint_lowfield, seekwater_specpoint_upfield] = computer_inverse_chemshift_vecotr_1_0([seekwater_chemshift_lowfield,seekwater_chemshift_upfield],water_frequency,dwelltime);


				% 6.1 Check for waterpeak: Criteria
                                                                                                   % upfield:lowfield because upfield chemshift --> low point
				[waterpeak_found, waterpeak_specpoint] = seek_metabolite_0_2(spectrum_zf_real(seekwater_specpoint_upfield:seekwater_specpoint_lowfield),Noise, {'metpeak_height'},seekwater_metpeak_height_criterium_value);
				waterpeak_specpoint = waterpeak_specpoint + seekwater_specpoint_upfield - 1;     % because to the function only the chemshift region of interest is given. If function gives: peak at specpoint=1 
                                                                                                   % this is in reality the point described by the formula above


                                                                        
                                                                                             
                                                                                             
                % 6.D1 DEBUG MODE: Plot spectrum where water was sought
                figure
                plot(chemshift_vector_zf(seekwater_specpoint_upfield:seekwater_specpoint_lowfield),spectrum_zf_real((seekwater_specpoint_upfield:seekwater_specpoint_lowfield)),'k','Linewidth',1.6)
                hold on
                plot(chemshift_vector_zf(waterpeak_specpoint),(0:spectrum_zf_real(waterpeak_specpoint)/500:spectrum_zf_real(waterpeak_specpoint)),'r')             % plot line crossing spectrum where water was found.
                hold off
                set(gca,'XDir','reverse');

                title(sprintf('Waterpeak, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
                xlabel('Chemical Shift')
                ylabel('Signal')
                legend('Measured Spectrum','Water Peak','Location','Best')
                
                
                %pause
                close(gcf)                                                                                   
                % DEBUD MODE END                                                                         
                                                                                             
			
                
                
                
                
				% 6.2: Enlarge region if peak not found
                if(~waterpeak_found)
					seekwater_chemshift_lowfield = seekwater_chemshift_lowfield - seekwater_chemshift_stepsize;
					seekwater_chemshift_upfield = seekwater_chemshift_upfield + seekwater_chemshift_stepsize;
					seekwater_step = seekwater_step + 1;
                end

			end %while


            if(waterpeak_found && exist('out_dir','var'))
				fprintf(error_log_fid,'W1 WARNING: No waterpeak found. Assuming water suppression.');
            end



  
            
            %% 6.D2 DEBUG MODE: Plot spectrum with line where water was found
            
            if(waterpeak_found)
                figure
                plot(chemshift_vector_zf,spectrum_zf_real,'k','Linewidth',1.6)
                hold on
                plot(chemshift_vector_zf(waterpeak_specpoint),(0:spectrum_zf_real(waterpeak_specpoint)/500:spectrum_zf_real(waterpeak_specpoint)),'r','Linewidth',1.6)     % plot line crossing spectrum where water was found.
                hold off
                set(gca,'XDir','reverse');

                title(sprintf('Waterpeak, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
                xlabel('Chemical Shift')
                ylabel('Signal')
                legend('Measured Spectrum','Water Peak','Location','Best') 
                
                %pause
                close(gcf)
            end

		
            
            
			%% 7. Seek & Destroy NAA


			% 7.0 Define Searching Region


			if(waterpeak_found)
				seekNAA_chemshift_lowfield = chemshift_vector_zf(waterpeak_specpoint) - (idealwater_chemshift - idealNAA_chemshift) - seekNAA_chemshiftregion_water;
				seekNAA_chemshift_upfield = chemshift_vector_zf(waterpeak_specpoint) - (idealwater_chemshift - idealNAA_chemshift) + seekNAA_chemshiftregion_water;
			else
				seekNAA_chemshift_lowfield = idealNAA_chemshift - seekNAA_chemshiftregion_nowater;
				seekNAA_chemshift_upfield = idealNAA_chemshift + seekNAA_chemshiftregion_nowater;
			end
	

			% convert chemshift to point
			seekNAA_specpoint_lowfield = find(min(abs(chemshift_vector_zf - seekNAA_chemshift_lowfield)) == abs(chemshift_vector_zf - seekNAA_chemshift_lowfield));
			seekNAA_specpoint_upfield = find(min(abs(chemshift_vector_zf - seekNAA_chemshift_upfield)) == abs(chemshift_vector_zf - seekNAA_chemshift_upfield));
            
                      
            
            % 7.1 Search for NAA with real part
			[NAApeak_found, NAApeak_specpoint] = seek_metabolite_0_2(spectrum_zf_real(seekNAA_specpoint_upfield:seekNAA_specpoint_lowfield), Noise,{'metpeak_height'},seekNAA_metpeak_height_criterium_value);
            NAApeak_specpoint = NAApeak_specpoint + seekNAA_specpoint_upfield - 1;

            
            
            % 7.1D DEBUG MODE: PLOT REGION OF SPECTRUM WHERE TO SEARCH FOR NAA
                figure
                plot(chemshift_vector_zf(seekNAA_specpoint_upfield:seekNAA_specpoint_lowfield),spectrum_zf_real((seekNAA_specpoint_upfield:seekNAA_specpoint_lowfield)),'k','Linewidth',1.6)
                hold on
                plot(chemshift_vector_zf(NAApeak_specpoint),(0:spectrum_zf_real(NAApeak_specpoint)/500:spectrum_zf_real(NAApeak_specpoint)),'r')             % plot line crossing spectrum where water was found.
                hold off
                set(gca,'XDir','reverse');

                title(sprintf('NAApeak, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
                xlabel('Chemical Shift')
                ylabel('Signal')
                legend('Measured Spectrum','NAA Peak','Location','Best')
                
                
                %pause
                close(gcf)  
            % DEBUG MODE END
            
            
            
            % 7.2 Search for NAA with abs value
            if(~NAApeak_found)
                % Try again with absolute value
                [NAApeak_found, NAApeak_specpoint] = seek_metabolite_0_2(spectrum_zf_abs(seekNAA_specpoint_upfield:seekNAA_specpoint_lowfield),Noise, {'metpeak_height'},seekNAA_metpeak_height_criterium_value);
                NAApeak_specpoint = NAApeak_specpoint + seekNAA_specpoint_upfield - 1;    
                     
                
                
            % 7.2D DEBUG MODE: PLOT REGION OF SPECTRUM WHERE TO SEARCH FOR NAA
                figure
                plot(chemshift_vector_zf(seekNAA_specpoint_upfield:seekNAA_specpoint_lowfield),spectrum_zf_abs((seekNAA_specpoint_upfield:seekNAA_specpoint_lowfield)),'k','Linewidth',1.6)
                hold on
                plot(chemshift_vector_zf(NAApeak_specpoint),(0:spectrum_zf_abs(NAApeak_specpoint)/500:spectrum_zf_abs(NAApeak_specpoint)),'r')             % plot line crossing spectrum where water was found.
                hold off
                set(gca,'XDir','reverse');

                title(sprintf('NAApeak, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
                xlabel('Chemical Shift')
                ylabel('Signal')
                legend('Measured Spectrum','NAA Peak','Location','Best')
                
                
                %pause
                close(gcf)  
            % DEBUG MODE END
                
                
                if(~NAApeak_found)
                    if(exist('out_dir','var'))
                        fprintf(error_log_fid,'E1 ERROR: No NAA-peak found. Move to next voxel.');
                    end
                    continue
                end
            end



%pause




			%% 8. Search for "Basis" of NAA



			[NAA_basis_upfield_found,NAA_basis_lowfield_found, NAA_basis_upfield_specpoint, NAA_basis_lowfield_specpoint] = ...
            seek_basis_of_metabolite_0_1(spectrum_zf_real(seekNAA_specpoint_upfield:seekNAA_specpoint_lowfield),chemshift_vector_zf,'fixed_distance_to_peak', ...
            [seekNAAbasis_upfield_chemshift,seekNAAbasis_lowfield_chemshift],NAApeak_specpoint);

            

			if(~NAA_basis_upfield_found)
                if(exist('out_dir','var'))
                    fprintf(error_log_fid,'E2 ERROR: Basis of NAA-peak at high chemical shift side not found. Move to next voxel.');
                end
				continue
			end
			if(~NAA_basis_lowfield_found)
                if(exist('out_dir','var'))
                    fprintf(error_log_fid,'E3 ERROR: Basis of NAA-peak at low chemical shift side not found. Move to next voxel.');
                end
				continue
			end




			%% 9. Compute the Basis at the chemical shift of the peak; Compute Signal
            
            
			% Put a line through the points (NAA_basis_upfield_specpoint, spectrum_zf_real(...)) and (NAA_basis_lowfield_specpoint,spectrum_zf_real(...)) and compute the value of this line at NAA_peak_specpoint; 
            % NOT IMPLEMENTED YET
            % Now implemented: Average
                  
			NAA_Basis_Signal = mean([spectrum_zf_real(NAA_basis_upfield_specpoint),spectrum_zf_real(NAA_basis_lowfield_specpoint)]);
			Signal = spectrum_zf_real(NAApeak_specpoint) - NAA_Basis_Signal;



            %% 10. Compute SNR

			SNR_mat(x,y,z) = Signal/Noise;




			%% 11. save a plot of spectrum in the NAA region with NAApeak, Basis1,2, line connecting Basis1,2; Save plot of noise region


		    figure
			plot(chemshift_vector_zf(NAA_basis_upfield_specpoint:NAA_basis_lowfield_specpoint),spectrum_zf_real(NAA_basis_upfield_specpoint:NAA_basis_lowfield_specpoint),'k','Linewidth',1.6)
            hold on
            plot(chemshift_vector_zf(NAApeak_specpoint),spectrum_zf_real(NAApeak_specpoint),'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k')
            plot(chemshift_vector_zf(NAA_basis_upfield_specpoint),spectrum_zf_real(NAA_basis_upfield_specpoint),'o','MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','k')
            plot(chemshift_vector_zf(NAA_basis_lowfield_specpoint),spectrum_zf_real(NAA_basis_lowfield_specpoint),'o','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k')
            plot(chemshift_vector_zf(NAApeak_specpoint),NAA_Basis_Signal,'o','MarkerSize',6,'MarkerFaceColor','m','MarkerEdgeColor','k')
            hold off
            set(gca,'XDir','reverse');

		    title(sprintf('NAApeak, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
            xlabel('Chemical Shift')
            ylabel('Signal')
            legend('Measured Spectrum','NAA Peak','Basispoint hf NAA','Basispoint lf NAA','Basis NAA','Location','Best')      % hf = high field, lf = low field
            
			saveas(gcf,sprintf('%s/SNR_computations/x_%s_y_%s_z_%s_NAApeak.fig', out_dir,num2str(x),num2str(y),num2str(z)))
			saveas(gcf,sprintf('%s/SNR_computations/x_%s_y_%s_z_%s_NAApeak.jpg', out_dir,num2str(x),num2str(y),num2str(z)))

            
            
		    figure
			plot(chemshift_vector(Noise_specpoint_upfield:Noise_specpoint_lowfield),Noise_spectrum)
            set(gca,'XDir','reverse');
		    title(sprintf('Noise region, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
            xlabel('Chemical Shift')
            ylabel('Signal')
            
			saveas(gcf,sprintf('%s/SNR_computations/x_%s_y_%s_z_%s_Noise_region.fig', out_dir,num2str(x),num2str(y),num2str(z)))
			saveas(gcf,sprintf('%s/SNR_computations/x_%s_y_%s_z_%s_Noise_region.jpg', out_dir,num2str(x),num2str(y),num2str(z)))
            
            % DEBUG MODE
            %pause
            % DEBUG MODE END
    		close all


		end % x-loop
	end % y-loop
end % z-loop

pause off