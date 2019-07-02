function SNR_mat = compute_SNR_matrix(CSI_mat_time,oversampling,water_frequency,dwelltime,mask,out_dir)
%
% Function to compute the SNR for each voxel of a CSI-matrix. The matrix must have the size (ROW,COL,SLC,vecSize), so it must already be summed up over the channels
%

%% 0. PREPARATIONS, DEFINITIONS


% 0.1 ROW, COL, SLC, vecSize

ROW = size(CSI_mat_time,1);
COL = size(CSI_mat_time,2);
SLC = size(CSI_mat_time,3);
vecSize = size(CSI_mat_time,4);


% 0.1 chemshift vector
chemshift_vector = compute_chemshift_vector_1_0(water_frequency,dwelltime/oversampling,vecSize*oversampling);


% 0.2 Variables for seeking water

seekwater_chemshift_min = 4.55;		 % initial chemshift region to search for water
seekwater_chemshift_max = 4.8;
seekwater_chemshift_stepsize = 0.05;
seekwater_totalsteps = 3;
seekwater_metpeak_height_criterium_value = 30;   % max of waterpeak must be waterpeak_border_criterium_value times higher than at the border


% 0.3 Variables for seeking NAA

idealwater_chemshift = 4.7;
idealNAA_chemshift = 2.01;
seekNAA_chemshiftregion_water = 0.05;
seekNAA_chemshiftregion_nowater = 0.4;
seekNAA_metpeak_height_criterium_value = 10;

% 0.4 Create error log file
error_log_file = sprintf('%s_maps/SNR_computations/compute_SNR_error_log.txt',out_dir);
error_log_fid = fopen(error_log_file, 'w+');
fprintf(error_log_fid,'Error Log File of %s\nThese errors, sorted by z,y,x of matrix were found:\n\n',Patname)


% 0.5 SNR_mat

SNR_mat = zeros(ROW,COL,SLC);


% 0.6 Region for Noise computation

Noise_chemshift_low = chemshift_vector(end);
Noise_chemshift_high = chemshift_vector(end) + 2;
Noise_specpoint_low = find(min(abs(chemshift_vector - Noise_chemshift_low)) == abs(chemshift_vector - Noise_chemshift_low));
Noise_specpoint_high = find(min(abs(chemshift_vector - Noise_chemshift_high)) == abs(chemshift_vector - Noise_chemshift_high));



%% 0.D DEBUG MODE: How close are the chemshift_vector values to the seekwater_chemshift_min and max?; is the oversampling enough?
diff_seekwatermin_to_chemshiftvec = chemshift_vector(find(min(abs(chemshift_vector - seekwater_chemshift_min)) == abs(chemshift_vector - seekwater_chemshift_min))) - seekwater_chemshift_min
diff_seekwatermax_to_chemshiftvec = chemshift_vector(find(min(abs(chemshift_vector - seekwater_chemshift_max)) == abs(chemshift_vector - seekwater_chemshift_max))) - seekwater_chemshift_max
clear diff_seekwatermin_to_chemshiftvec diff_seekwatermax_to_chemshiftvec



%% 1. OVERSAMPLING

% squeeze nur in bestimmte richtung möglich?

CSI_mat_time_os = zeros(ROW,COL,SLC,oversampling*vecSize);
CSI_mat_time_os(:,:,:,1:vecSize) = CSI_mat_time;




%% 2. Time Fourier Transform

CSI_mat_freq_os = fftshift(fft(CSI_mat_time_os, [],4),4);



%% 3. Process each voxel individually: LOOPS

for z = 1:SLC
	for y = 1:COL
		for x = 1:ROW

			if(mask(x,y) == 0)
				continue
			end
			spectrum = CSI_mat_freq_os(x,y,z,:);



			%% 4. Write Voxel in error log file

			fprintf(error_log_fid,'z_%d_y_%d_x_%d\n',z,y,x);





			%% 5. Search for water_peak in each voxel

			% 5.0 Define searching region		
			
			% initialization
			waterpeak_found = 0;
			seekwater_step = 0;

			% while peak not found increase region to look for peak and check again
			while(waterpeak_found == 0 && seekwater_step < seekwater_totalsteps-1)

				% convert chemshift to points
				% The exact chemshift is most probable not within the chemshift_vector, so search for the minimum
				seekwater_specpoint_min = find(min(abs(chemshift_vector - seekwater_chemshift_min)) == abs(chemshift_vector - seekwater_chemshift_min))
				seekwater_specpoint_max = find(min(abs(chemshift_vector - seekwater_chemshift_max)) == abs(chemshift_vector - seekwater_chemshift_max))
				%[seekwater_specpoint_min, seekwater_specpoint_max] = computer_inverse_chemshift_vecotr_1_0([seekwater_chemshift_min,seekwater_chemshift_max],water_frequency,dwelltime);


				% 5.1 Check for waterpeak: Criteria

				[waterpeak_found, waterpeak_specpoint] = seek_metabolite_0_1(spectrum(seekwater_specpoint_max:seekwater_specpoint_min), {'metpeak_height'},[seekwater_metpeak_height_criterium_value]);  % max:min because high chemshift --> low point
				waterpeak_specpoint = waterpeak_specpoint + seekwater_specpoint_max - 1; % because to the function only the chemshift region of interest is given. If function says: peak at specpoint=1 this is in reality the point described by the formula above


			
				% 5.2: Enlarge region if peak not found
				if(waterpeak_found == 0)
					seekwater_chemshift_min = seekwater_chemshift_min - seekwater_chemshift_stepsize;
					seekwater_chemshift_max = seekwater_chemshift_max + seekwater_chemshift_stepsize;
					seekwater_step = seekwater_step + 1;
				end

			end %while


			if(waterpeak_found == 0)
				fprintf(error_log_fid,'W1 WARNING: No waterpeak found. Assuming water suppression.');
			end




		
			%% 6. Seek & Destroy NAA


			% 6.0 Define Searching Region


			if(waterpeak_found)
				seekNAA_chemshift_min = chemshift_vector(waterpeak_specpoint) - (idealwater_chemshift - idealNAA_chemshift) - seekNAA_chemshiftregion_water;
				seekNAA_chemshift_max = chemshift_vector(waterpeak_specpoint) - (idealwater_chemshift - idealNAA_chemshift) + seekNAA_chemshiftregion_water;
			else
				seekNAA_chemshift_min = idealNAA_chemshift - seekNAA_chemshiftregion_nowater;
				seekNAA_chemshift_max = idealNAA_chemshift + seekNAA_chemshiftregion_nowater;
			end
	

			% convert chemshift to point
			seekNAA_specpoint_min = find(min(abs(chemshift_vector - seekNAA_chemshift_min)) == abs(chemshift_vector - seekNAA_chemshift_min))
			seekNAA_specpoint_max = find(min(abs(chemshift_vector - seekNAA_chemshift_max)) == abs(chemshift_vector - seekNAA_chemshift_max))




			[NAApeak_found, NAApeak_specpoint] = seek_metabolite_0_1(spectrum(seekNAA_specpoint_max:seekNAA_specpoint_min), {'metpeak_height'},[seekNAA_metpeak_height_criterium_value]);


			if(~NAApeak_found)
				fprintf(error_log_fid,'E1 ERROR: No NAA-peak found. Move to next voxel.');
				continue
			end





			%% 7. Search for "Basis" of NAA


			[NAA_basis_highchemshift_found, NAA_basis_highchemshift_specpoint] = seek_basis_of_metabolite_0_1(spectrum(seekNAA_specpoint_max:seekNAA_specpoint_min),NAApeak_specpoint,'',valuess,'high_chemshift');
			[NAA_basis_lowchemshift_found, NAA_basis_lowchemshift_specpoint] = seek_basis_of_metabolite_0_1(spectrum(seekNAA_specpoint_max:seekNAA_specpoint_min),NAApeak_specpoint,'',valuess,'low_chemshift');


			if(~NAA_basis_highchemshift_found)
				fprintf(error_log_fid,'E2 ERROR: Basis of NAA-peak at high chemical shift side not found. Move to next voxel.');
				continue
			end
			if(~NAA_basis_lowchemshift_found)
				fprintf(error_log_fid,'E3 ERROR: Basis of NAA-peak at low chemical shift side not found. Move to next voxel.');
				continue
			end




			%% 8. Compute the Basis at the chemical shift of the peak & Signal

			% Put a line through the points (NAA_basis_highchemshift_specpoint, spectrum(...)) and (NAA_basis_lowchemshift_specpoint,spectrum(...)) and compute the value of this line at NAA_peak_specpoint

			NAA_Basis = spectrum(


			Signal = spectrum(NAApeak_specpoint) - NAA_Basis;



			%% 9. Compute the standard deviation of the noise and the SNR


			Noise = std(spectrum(Noise_specpoint_high:Noise_specpoint_low));

			SNR_mat(x,y,z) = Signal/Noise;




			%% 10. save a plot of spectrum in the NAA region with NAApeak, Basis1,2, line connecting Basis1,2; Save plot of noise region




		end % x-loop
	end % y-loop
end % z-loop



%% 11. Write the matrix to a file

SNR_file_fid = fopen(sprintf('%s_maps/matlab_SNR_map.raw',out_dir),'w+');
fwrite(SNR_file_fid,SNR_mat,'float');


fclose(SNR_file_fid);
fclose(error_log_fid);
