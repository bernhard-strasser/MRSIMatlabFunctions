function read_and_plot_phasemaps_1_1(file_struct,plot_mode,pha_abs,Hamming,LCM_add_phase)


% file_struct = struct('name', 'CSI', 'Image', 'Flip', 'LCM-phamap', 'Weighting', 'water_frequency', 'dwelltime', 'fredir_shift')


if(~exist('plot_mode','var') || ~strcmpi(plot_mode,'bandwidth') && ~strcmpi(plot_mode,'Compare_CSI_Imaging'))
    display( [char(10) 'What mode of the program do you want to use?' char(10) '(1) Compare Imaging of different bandwidths' char(10) '(2) Compare CSI and Imaging phases' char(10)])
    plot_mode = input('plot_mode = ');
    if(plot_mode == 1)
        plot_mode = 'bandwidth';
    elseif(plot_mode == 2)
        plot_mode = 'Compare_CSI_Imaging';
    else
        display( [char(10) 'Wrong input. Computer will explode. Drink tea and leave work.' char(10)])
    end
end


if(exist('pha_abs', 'var') && strcmpi(pha_abs,'abs'))
    pha_abs = 'abs';
else
    pha_abs = 'angle';
end



%% 1. if you want to compare water_images of different bandwidths
if(strcmpi(plot_mode, 'bandwidth'))
    
    
    size(file_struct,2)
    
    % 1.0 DEFINTIONS
    [headersize,k_fredir_center,fredir_measured,phadir_measured,total_channel_no] = read_image_dat_meas_header_1_0(file_struct(1).Image);
    ROW = phadir_measured; 
    COL = phadir_measured; % no bug here, we always have square images
    SLC = 1;
    Image_Normal = zeros(size(file_struct,2),total_channel_no,ROW,COL,SLC);
    Image_Flip = zeros(size(file_struct,2),total_channel_no,ROW,COL,SLC);
    Image_Corr = zeros(size(file_struct,2),total_channel_no,ROW,COL,SLC);
    
    
    
    
    % 1.1 READ DATA
    for read_index = 1:size(file_struct,2)
        
        Image_Normal(read_index,:,:,:,:) = read_image_dat_1_6_Flip_test(file_struct(read_index).Image,128,128,0,0,'AP',0);

        if(~strcmp(file_struct(read_index).Flip,'NO'))
            Image_Flip(read_index,:,:,:,:) = read_image_dat_1_6_Flip_test(file_struct(read_index).Flip,128,128,0,0,'AP',0);
            %Image_Flip(read_index,:,:,:,:) = flipdim(flipdim(Image_Flip(read_index,:,:,:,:),3),4);
            %Image_Flip(read_index,:,:,:,:) = circshift(Image_Flip(read_index,:,:,:,:), [0 0 -1 1]);
            
            Image_Corr(read_index,:,:,:,:) = (Image_Normal(read_index,:,:,:,:) + Image_Flip(read_index,:,:,:,:).*abs(Image_Normal(read_index,:,:,:,:)) ./ abs(Image_Normal(read_index,:,:,:,:)))/2;
        end
    end
     

    
    
    % 1.2 PLOT DATA
    plot_channel = 1;
    while(plot_channel ~= 666)
        
        display( [char(10) 'Enter channel you want to plot. Enter 666 for continueing with further program' char(10)] )
        plot_channel = input('channel = ');
        
        if(plot_channel == 666)
            continue
        end
        
        
        for plot_index = 1:size(file_struct,2)

            % Plot Image Normal, Flip, Corr
            figure;
            imagesc(rad2deg(angle(squeeze(Image_Normal(plot_index,plot_channel,:,:,1)))),[-180 180])
            title(sprintf('ImageNormal %s channel %i', file_struct(plot_index).name, plot_channel))
           
            if(~strcmp(file_struct(plot_index).Flip,'NO'))
                figure;
                imagesc(rad2deg(angle(squeeze(Image_Flip(plot_index,plot_channel,:,:,1)))),[-180 180])
                title(sprintf('ImageFlip %s channel %i', file_struct(plot_index).name, plot_channel))          
                
                figure;
                imagesc(rad2deg(angle(squeeze(Image_Corr(plot_index,plot_channel,:,:,1)))),[-180 180])
                title(sprintf('ImageCorr %s channel %i', file_struct(plot_index).name, plot_channel))
            end
            
            
            
            
            
            % Plot Subtraction maps Normal, Flip, Corr
            figure;
            imagesc(rad2deg(angle(squeeze(Image_Normal(1,plot_channel,:,:,1)))) - rad2deg(angle(squeeze(Image_Normal(plot_index,plot_channel,:,:,1)))),[-15 15])
            title(sprintf('ImageNormal Sub %s - %s channel %i', file_struct(1).name,file_struct(plot_index).name, plot_channel))    
        
            if(~strcmp(file_struct(plot_index).Flip,'NO'))

                figure;
                imagesc(rad2deg(angle(squeeze(Image_Flip(1,plot_channel,:,:,1)))) - rad2deg(angle(squeeze(Image_Flip(plot_index,plot_channel,:,:,1)))),[-15 15])
                title(sprintf('ImageFlip Sub %s - %s channel %i', file_struct(1).name,file_struct(plot_index).name, plot_channel))      

                figure;
                imagesc(rad2deg(angle(squeeze(Image_Corr(1,plot_channel,:,:,1)))) - rad2deg(angle(squeeze(Image_Corr(plot_index,plot_channel,:,:,1)))),[-15 15])
                title(sprintf('ImageCorr Sub %s - %s channel %i', file_struct(1).name,file_struct(plot_index).name, plot_channel))  
            end

            waitforbuttonpress
            close all;

            
            
            
        end
  
    end
    
    
    
    
    
    
    
    
    
    
    
%% 2. if you want to compare csi with imaging phasemaps
elseif(strcmpi(plot_mode, 'Compare_CSI_Imaging'))
%%

 







    for read_index = 1:size(file_struct,2)  % CSI DATA = LARGE --> FOR COMPUTATIONS FOR EACH FILE INDIVIDUALLY
   
        
        
        
        
        % 2.1 DEFINITIONS
        [ROW,COL,SLC,vecSize,total_channel_no,total_k_points,headersize] = read_csi_dat_meas_header_1_0(file_struct(read_index).CSI);      
        water_frequency = file_struct(read_index).water_frequency;
        dwelltime = file_struct(read_index).dwelltime;

        % 2.2 READ IN DATA & COMPUTE PHASE AND PHASE DIFFERENCES
        
        % CSI, hamming it
        CSI = read_csi_dat_1_3(file_struct(read_index).CSI);
        if(strcmpi(Hamming,'hamming'))
            CSI = Hamming_filtering_1_1(CSI,[2 3]);
        end
		CSI_phase = rad2deg(angle(squeeze(CSI(:,:,:,1,1))));
        
        
        % Image NORMAL
        Image_Normal = read_image_dat_1_6_Flip_test(file_struct(read_index).Image,ROW,COL,0,file_struct(read_index).fredir_shift,'AP',0);
        Image_Normal = circshift(Image_Normal, [0 0 0]);
		Im_Normal_Phase = rad2deg(angle(squeeze(Image_Normal(:,:,:,1))));
		Sub_CSI_ImNormal = CSI_phase - Im_Normal_Phase; 
		Sub_CSI_ImNormal = unwrap_phase_1_0(Sub_CSI_ImNormal);
        
        % Image FLIP, COMPUTE CORR
        if(~strcmp(file_struct(read_index).Flip,'NO'))
            Image_Flip = read_image_dat_1_6_Flip_test(file_struct(read_index).Flip,ROW,COL,0,file_struct(read_index).fredir_shift,'AP',1);
            %Image_Flip = flipdim(flipdim(Image_Flip,2),3);
            Image_Flip = circshift(Image_Flip, [0 0 0]);
            Image_Corr = ( Image_Normal + Image_Flip .* abs(Image_Normal)./abs(Image_Flip))/2;
        
			Im_Flip_Phase = rad2deg(angle(squeeze(Image_Flip(:,:,:,1))));
			Im_Corr_Phase = rad2deg(angle(squeeze(Image_Corr(:,:,:,1))));

			Sub_CSI_ImFlip = CSI_phase - Im_Flip_Phase;
			Sub_CSI_ImCorr = CSI_phase - Im_Corr_Phase;

			Sub_CSI_ImFlip = unwrap_phase_1_0(Sub_CSI_ImFlip);
			Sub_CSI_ImCorr = unwrap_phase_1_0(Sub_CSI_ImCorr);

		end 
        
        
        
        % LCM phasemap
        if(~strcmpi(file_struct(read_index).LCM_phasemap,'NO'))
            LCM_phasemap_fid = fopen(file_struct(read_index).LCM_phasemap,'r');
            LCM_phasemap = fread(LCM_phasemap_fid, [64,64], 'double');
            LCM_phasemap = wrapTo180(LCM_phasemap + LCM_add_phase);
            fclose(LCM_phasemap_fid);
			Sub_CSI_LCM = unwrap_phase_1_ß(CSI_phase - LCM_phasemap);
        end
        
        
        % WEIGHTING
        if(~strcmp(file_struct(read_index).Weighting,'NO'))
            Image_Weighting = abs(read_image_dat_1_5(file_struct(read_index).Weighting));
        else
            Image_Weighting = sum(abs(Image_Normal),1);
        end        
    
    
        
        

        
        
        % PART I: Phasemaps
        % Plot: CSI-Phasemp, Image_phasemap, Image_Flip_Phasemap, Image_Corr_phasemap, Sub_CSI_to_Image, Sub_CSI_to_Flip, Sub_CSI_to_Corr



        plot_channel = 1;
        while(plot_channel ~= 666)
            
            display( [char(10) 'Enter channel you want to plot. Enter 666 for continueing with further program' char(10)] )
            plot_channel = input('channel = ');
            close all;
            
            if(plot_channel == 666)
                continue
            end           	



            % PLOT CSI NORMAL FLIP CORR 
            figure;
            eval([ 'imagesc(-rad2deg(' pha_abs '(squeeze(CSI(plot_channel,:,:,1,1)))))' ]);
            title(sprintf('CSI %s channel %i', file_struct(read_index).name, plot_channel))

            figure;
            eval([ 'imagesc(rad2deg(' pha_abs '(squeeze(Image_Normal(plot_channel,:,:,1)))))' ]);
            title(sprintf('ImageNormal %s channel %i', file_struct(read_index).name, plot_channel))

            if(~strcmp(file_struct(read_index).Flip,'NO'))
                figure;
                eval([ 'imagesc(rad2deg(' pha_abs '(squeeze(Image_Flip(plot_channel,:,:,1)))))' ]);
                title(sprintf('ImageFlip %s channel %i', file_struct(read_index).name, plot_channel))            
            
                figure;
                eval([ 'imagesc(rad2deg(' pha_abs '(squeeze(Image_Corr(plot_channel,:,:,1,1)))))' ]);
                title(sprintf('ImageCorr %s channel %i', file_struct(read_index).name, plot_channel))

            end
            
            
            
            % PLOT LCM PHAMAP
            
            if(~strcmp(file_struct(read_index).LCM_phasemap,'NO'))       
            
                figure;
                imagesc(LCM_phasemap(:,:))
                title(sprintf('LCM phasemap %s channel %i', file_struct(read_index).name, plot_channel))

            end     
            
            
            
            
            % PLOT SUB CSI - NORMAL, FLIP, CORR
            figure;
            imagesc(rad2deg(-angle(squeeze(CSI(plot_channel,:,:,1,1)))) - rad2deg(angle(squeeze(Image_Normal(plot_channel,:,:,1)))),[-60 60])
            title(sprintf('Sub CSI - ImageNormal %s channel %i', file_struct(read_index).name, plot_channel))  

            if(~strcmp(file_struct(read_index).Flip,'NO'))
                figure;
                imagesc(rad2deg(-angle(squeeze(CSI(plot_channel,:,:,1,1)))) - rad2deg(angle(squeeze(Image_Flip(plot_channel,:,:,1)))),[-60 60])
                title(sprintf('Sub CSI - ImageFlip %s channel %i', file_struct(read_index).name, plot_channel)) 

                figure;
                imagesc(rad2deg(-angle(squeeze(CSI(plot_channel,:,:,1,1)))) - rad2deg(angle(squeeze(Image_Corr(plot_channel,:,:,1)))),[-60 60])
                title(sprintf('Sub CSI - ImageCorr %s channel %i', file_struct(read_index).name, plot_channel))         
            end
            
            
            
            % PLOT LCM - CORR, CSI
            
            if(~strcmp(file_struct(read_index).LCM_phasemap,'NO'))      
            
                if(~strcmp(file_struct(read_index).Flip,'NO'))
                    figure;
                    imagesc(LCM_phasemap - rad2deg(angle(squeeze(Image_Corr(plot_channel,:,:,1)))),[-60 60])
                    title(sprintf('SUB LCM - ImageCorr %s channel %i', file_struct(read_index).name, plot_channel))
                else
                    figure;
                    imagesc(LCM_phasemap - rad2deg(angle(squeeze(Image_Normal(plot_channel,:,:,1)))),[-60 60])
                    title(sprintf('SUB LCM - ImageNormal %s channel %i', file_struct(read_index).name, plot_channel)) 
                end
                
                
                figure;
                imagesc(-LCM_phasemap - rad2deg(angle(squeeze(CSI(plot_channel,:,:,1,1)))),[-60 60])
                title(sprintf('SUB LCM - CSI %s channel %i', file_struct(read_index).name, plot_channel))                
                

            end             
            
            
            
            
        end
    
        








		if(~strcmpi(mask_left_down_corner,'NO'))
			% PART II: Statistical computations of subtractions.

			% Let user define a mask where the mean and standard deviation of the difference maps are computed.


			if(~(exist('mask_left_down_corner', 'var') && exist('mask_left_up_corner', 'var') && exist('mask_right_down_corner','var') && exist('mask_right_up_corner','var')))
				display('Please give me the corners of the mask for statistical computations in voxel units')
			    mask_left_border = input('mask_left_border = ');
			    mask_right_border = input('mask_left_border = ');
			    mask_down_border = input('mask_down_border = ');
			    mask_up_border = input('mask_up_border = ');
		    end

			
			mask = zeros(size(CSI_phase,2),size(CSI_phase,3));
			mask(mask_up_border:mask_down_border,mask_left_border:mask_right_border) = 1;




    
        
        
        
        
       
        
        % PART III: Spectra
        
        % 2.II.1: Preparations
        close all;
        if(strcmp(file_struct(read_index).Flip,'NO'))   % If just normal image is inputted, then just use this for phasing
            Phasing_index = 1;
        else
            Phasing_index = 3;
        end
        
        
        plot_x = 1; plot_y = 1;
        % 2.II.2: LOOP FOR VOXEL TO PLOT AND USER INPUT
        while(plot_x ~= 666 && plot_y ~= 666)      % let user define which voxel to show
            
            display( [char(10) 'Enter x, y you want to plot. Enter x = 666 or y = 666 for continueing with further program' char(10)] )
            plot_x = input('x = ');
            plot_y = input('y = ');
            plot_x_y = [plot_x,plot_y];
            
            if(plot_x == 666)
                continue
            end 
        
            
            
            
            % 2.II.3: LOOP SO THAT NORMAL IMAGE, FLIP AND CORR ARE USED FOR PHASING
            for Use_Normal_Flip_Corr_To_Phase = 1:Phasing_index

                if(Use_Normal_Flip_Corr_To_Phase == 1)
                    weighting_n_phasing_factor = Image_Normal(:,plot_x_y(1),plot_x_y(2),1);
                elseif(Use_Normal_Flip_Corr_To_Phase == 2)
                    weighting_n_phasing_factor = Image_Flip(:,plot_x_y(1),plot_x_y(2),1);
                elseif(Use_Normal_Flip_Corr_To_Phase == 3)
                    weighting_n_phasing_factor = Image_Corr(:,plot_x_y(1),plot_x_y(2),1);
                end


                
                
                % 2.II.4: COMPUTE WEIGHTING AND PHASING FACTOR
                weighting_sum_square_mag = sum(abs(Image_Normal(:,plot_x_y(1),plot_x_y(2),1)).^2,1);
                weighting_normalize = Image_Weighting(1,plot_x_y(1),plot_x_y(2)) ./ weighting_sum_square_mag;
                weighting_n_phasing_factor = repmat(weighting_normalize, [total_channel_no 1 1 1]) .* weighting_n_phasing_factor; 


                % 2.II.4: WEIGHT AND PHASE THE CSI-DATA
                CSI_phased = repmat(weighting_n_phasing_factor,[1 1 1 1 vecSize]) .* CSI(:,plot_x_y(1),plot_x_y(2),1,:);
                CSI_summed = sum(CSI_phased,1);



                % 2.II.5: FFT
                CSI_unphased = fftshift(fft(CSI(:,plot_x_y(1),plot_x_y(2),1,:),vecSize,5),5);
                CSI_phased = fftshift(fft(CSI_phased(:,plot_x_y(1),plot_x_y(2),1,:),vecSize,5),5);
                CSI_summed = fftshift(fft(CSI_summed(plot_x_y(1),plot_x_y(2),1,:),vecSize,4),4);


                
                
                
                
                % 2.II.6: plot unphased & phased for all channel, and summed 
                chemshift_vec = compute_chemshift_vector_1_0(water_frequency,dwelltime,vecSize);
                for plot_channel = 1:total_channel_no

                    figure;
                    plot(chemshift_vec,CSI_unphased(plot_channel,1,1,1,:))
                    set(gca,'XDir','reverse');
                    title(sprintf('CSI unphased %s, channel %d, [x,y]=[%d,%d]', file_struct(read_index).name, plot_channel,plot_x_y(1),plot_x_y(2)))

                    figure;
                    plot(chemshift_vec,CSI_phased(plot_channel,1,1,1,:))
                    set(gca,'XDir','reverse');
                    title(sprintf('CSI phased %s, channel %d, [x,y]=[%d,%d]', file_struct(read_index).name, plot_channel,plot_x_y(1),plot_x_y(2)))               

                end

                figure;
                plot(chemshift_vec,CSI_summed(1,1,1,:))
                set(gca,'XDir','reverse');
                title(sprintf('CSI summed %s, [x,y]=[%d,%d]', file_struct(read_index).name,plot_x_y(1),plot_x_y(2)))   
  
                
            end % Normal_Flip_Corr Phasing
    
    
        end % x_y_plot
        
  
        
        
        
        
        
        
clear all;
    end % CSI_file_to_compute
    
    

    
    % save all open figures to save_fig_dir
    
    
    
    
end % if(plotmode)

close all
