function SNR_mat = compute_SNR_matrix(Patname,CSI_mat_time,oversampling,water_frequency,dwelltime)
%
% Function to compute the SNR for each voxel of a CSI-matrix. The matrix must have the size (ROW,COL,SLC,vecSize), so it must already be summed up over the channels
%

%% 0. PREPARATIONS, DEFINITIONS

seekwater_enlargeregion_stepsize = 0.05;
seekwater_enlargeregion_totalsteps = 3;

seekwater_border_criterium_value = 30;   % max of waterpeak must be waterpeak_border_criterium_value times higher than at the border

% create error log file
error_log_file = sprintf('./tmp/compute_SNR_error_log_%s.txt',Patname);
error_log_fid = fopen(error_log_file, 'w+');
fprintf(error_log_fid,'Error Log File of %s\nThese errors, sorted by z,y,x of matrix were found:\n\n',Patname)


%% 1. OVERSAMPLING

% squeeze nur in bestimmte richtung möglich?

CSI_mat_time_os = zeros(size(CSI_mat_time,1),size(CSI_mat_time,2),size(CSI_mat_time,3),oversampling*size(CSI_mat_time,4));
CSI_mat_time_os(:,:,:,1:size(CSI_mat_time,4)) = CSI_mat_time;




%% 2. Time Fourier Transform

CSI_mat_freq_os = fftshift(fft(CSI_mat_time_os, [],4),4);



%% 3. Process each voxel individually: LOOPS

for z = 1:size(CSI_mat_time,3)
	for y = 1:size(CSI_mat_time,2)
		for x = 1:size(CSI_mat_time,1)




%% 4. Write Voxel in error log file

			fprintf(error_log_fid,'z_%d_y_%d_x_%d\n',z,y,x);




			%% 5. Search for water_peak in each voxel


			% 5.0 Define searching region
			
			
			% initialization
			seekwater_chemshift_min = 4.55;
			seekwater_chemshift_max = 4.8;
			waterpeak_found = 0;
			water_chemshift_step = 0;

			% while peak not found increase region to look for peak and check again
			while(waterpeak_found == 0 && water_chemshift_step < water_chemshift_totalsteps-1)

				% convert chemshift to points
				[specpoint_water_min, specpoint_water_max] = computer_inverse_chemshift_vecotr_1_0([seekwater_chemshift_min,seekwater_chemshift_max],water_frequency,dwelltime);



				% 5.1 Check for waterpeak: Criteria
				waterpeak_higher_than_leftborder = logical(max(CSI_mat_freq_os(x,y,z,specpoint_water_max:specpoint_water_min)) > 30*CSI_mat_freq_os(x,y,z,specpoint_water_max));						% max:min because high chemshift --> low point
				waterpeak_higher_than_rightborder = logical(max(CSI_mat_freq_os(x,y,z,specpoint_water_max:specpoint_water_min)) > 30*CSI_mat_freq_os(x,y,z,specpoint_water_min));
				% More criteria needed ???
				waterpeak_found = waterpeak_higher_than_leftborder && waterpeak_higher_than_rightborder;

				% if peak not found --> enlarge region to search for peak
				if(waterpeak_found == 0)
					seekwater_chemshift_min = seekwater_chemshift_min - seekwater_enlargeregion_stepsize;
					seekwater_chemshift_max = seekwater_chemshift_max + seekwater_enlargeregion_stepsize;
				end

			end %while


			if(waterpeak_found == 0)
				fprintf(error_log_fid,'W1 WARNING: No waterpeak found. Assuming water suppression.');
			end








fclose(error_log_fid);
