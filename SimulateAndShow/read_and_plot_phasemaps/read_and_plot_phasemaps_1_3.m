function read_and_plot_phasemaps_1_3(file_struct,plot_mode,Hamming,LCM_add_phase)


% file_struct = struct('name', 'CSI', 'Image', 'Flip', 'LCM', 'Weighting', 'water_frequency', 'dwelltime', 'fredir_shift', 'mask_position', 'out_dir')



%% 0. Definitions, Preparations


% FLAGS

Flip_flag = 0;
LCM_flag = 0;
Weighting_flag = 0;
out_dir_flag = 0;
if(~strcmpi(file_struct(1).Flip,'NO'))
    Flip_flag = 1;
end
if(~strcmpi(file_struct(1).LCM,'NO'))
    LCM_flag = 1;
end
if(~strcmpi(file_struct(1).Weighting,'NO'))
    Weighting_flag = 1;
end
if(~strcmpi(file_struct(1).out_dir,'NO'))
    out_dir_flag = 1;
end


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
    if(out_dir_flag)
        out_dir = file_struct(read_index).out_dir;
    end
    
    
    
    % 1.1 READ DATA
    for read_index = 1:size(file_struct,2)
        
        Image_Normal(read_index,:,:,:,:) = read_image_dat_1_6_Flip_test(file_struct(read_index).Image,128,128,0,0,'AP',0);

        if(Flip_flag)
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
            
            if(out_dir_flag)
                saveas(gcf,sprintf('%s/ImNormal_%s_channel%02d.fig', out_dir,file_struct(plot_index).name,plot_channel))
                saveas(gcf,sprintf('%s/ImNormal_%s_channel%02d.jpg', out_dir,file_struct(plot_index).name,plot_channel))
            end
          
            if(Flip_flag)
                fig1 = figure;
                imagesc(rad2deg(angle(squeeze(Image_Flip(plot_index,plot_channel,:,:,1)))),[-180 180])
                title(sprintf('ImageFlip %s channel %i', file_struct(plot_index).name, plot_channel))  
                
                
                fig2 = figure;
                imagesc(rad2deg(angle(squeeze(Image_Corr(plot_index,plot_channel,:,:,1)))),[-180 180])
                title(sprintf('ImageCorr %s channel %i', file_struct(plot_index).name, plot_channel))
                
                if(out_dir_flag)
                    saveas(fig1,sprintf('%s/ImFlip_%s_channel%02d.fig', out_dir,file_struct(plot_index).name,plot_channel))
                    saveas(fig1,sprintf('%s/ImFlip_%s_channel%02d.jpg', out_dir,file_struct(plot_index).name,plot_channel))
                    saveas(fig2,sprintf('%s/ImCorr_%s_channel%02d.fig', out_dir,file_struct(plot_index).name,plot_channel))
                    saveas(fig2,sprintf('%s/ImCorr_%s_channel%02d.jpg', out_dir,file_struct(plot_index).name,plot_channel))
                end
                
                
            end
            
            
            
            
            
            % Plot Subtraction maps Normal, Flip, Corr
            figure;
            imagesc(rad2deg(angle(squeeze(Image_Normal(1,plot_channel,:,:,1)))) - rad2deg(angle(squeeze(Image_Normal(plot_index,plot_channel,:,:,1)))),[-15 15])
            title(sprintf('ImageNormal Sub %s - %s channel %i', file_struct(1).name,file_struct(plot_index).name, plot_channel))  
            if(out_dir_flag)
            	saveas(gcf,sprintf('%s/SubNormal_%s_minus_%s_channel%02d.fig', out_dir,file_struct(1).name,file_struct(plot_index).name,plot_channel))
                saveas(gcf,sprintf('%s/SubNormal_%s_minus_%s_channel%02d.jpg', out_dir,file_struct(1).name,file_struct(plot_index).name,plot_channel))
            end          
            
            
            
        
            if(Flip_flag)

                fig1 = figure;
                imagesc(rad2deg(angle(squeeze(Image_Flip(1,plot_channel,:,:,1)))) - rad2deg(angle(squeeze(Image_Flip(plot_index,plot_channel,:,:,1)))),[-15 15])
                title(sprintf('ImageFlip Sub %s - %s channel %i', file_struct(1).name,file_struct(plot_index).name, plot_channel))
                
                

                fig2 = figure;
                imagesc(rad2deg(angle(squeeze(Image_Corr(1,plot_channel,:,:,1)))) - rad2deg(angle(squeeze(Image_Corr(plot_index,plot_channel,:,:,1)))),[-15 15])
                title(sprintf('ImageCorr Sub %s - %s channel %i', file_struct(1).name,file_struct(plot_index).name, plot_channel))
                if(out_dir_flag)
                    saveas(fig1,sprintf('%s/SubFlip_%s_minus_%s_channel%02d.fig', out_dir,file_struct(1).name,file_struct(plot_index).name,plot_channel))
                    saveas(fig1,sprintf('%s/SubFlip_%s_minus_%s_channel%02d.jpg', out_dir,file_struct(1).name,file_struct(plot_index).name,plot_channel))
                    saveas(fig2,sprintf('%s/SubCorr_%s_minus_%s_channel%02d.fig', out_dir,file_struct(1).name,file_struct(plot_index).name,plot_channel))
                    saveas(fig2,sprintf('%s/SubCorr_%s_minus_%s_channel%02d.jpg', out_dir,file_struct(1).name,file_struct(plot_index).name,plot_channel))
                end                   
                
                
            end

            waitforbuttonpress
            close all;

            
            
            
        end
  
    end
    
    
    
    
    
    
    
    
    
    
    
%% 2. if you want to compare csi with imaging phasemaps
elseif(strcmpi(plot_mode, 'Compare_CSI_Imaging'))
%%

 







    for read_index = 1:size(file_struct,2)  % CSI DATA = LARGE --> FOR COMPUTATIONS FOR EACH FILE INDIVIDUALLY
   
        
        if(out_dir_flag)
            out_dir = file_struct(read_index).out_dir
        end
        sprintf('Processing %s',file_struct(read_index).name)  
        
        
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
        CSI = conj(CSI);
        CSI_phase = -rad2deg(angle(squeeze(CSI(:,:,:,1,1))));
        
        
        % Image NORMAL
        Image_Normal = read_image_dat_1_6_Flip_test(file_struct(read_index).Image,ROW,COL,0,file_struct(read_index).fredir_shift,'AP',0);
        Image_Normal = circshift(Image_Normal, [0 1 -1]);
        Im_Normal_Phase = rad2deg(angle(squeeze(Image_Normal(:,:,:,1))));
		Sub_CSI_ImNormal = CSI_phase - Im_Normal_Phase; 
		Sub_CSI_ImNormal = unwrap_phase_1_0(Sub_CSI_ImNormal);

        
        % Image FLIP, COMPUTE CORR
        if(Flip_flag)
            Image_Flip = read_image_dat_1_6_Flip_test(file_struct(read_index).Flip,ROW,COL,0,file_struct(read_index).fredir_shift,'AP',1);
            %Image_Flip = flipdim(flipdim(Image_Flip,2),3);
            Image_Flip = circshift(Image_Flip, [0 1 -1]);
            Image_Corr = ( Image_Normal + Image_Flip .* abs(Image_Normal)./abs(Image_Flip))/2;
            
            Im_Flip_Phase = rad2deg(angle(squeeze(Image_Flip(:,:,:,1))));
			Im_Corr_Phase = rad2deg(angle(squeeze(Image_Corr(:,:,:,1))));

			Sub_CSI_ImFlip = CSI_phase - Im_Flip_Phase;
			Sub_CSI_ImCorr = CSI_phase - Im_Corr_Phase;

			Sub_CSI_ImFlip = unwrap_phase_1_0(Sub_CSI_ImFlip);
			Sub_CSI_ImCorr = unwrap_phase_1_0(Sub_CSI_ImCorr);
        end 
        
        
        
        % LCM phasemap
        if(LCM_flag)
            
            LCM_phasemap_files = file_struct(read_index).LCM_phasemap;
            LCM_phasemap = zeros(size(LCM_phasemap_files,2),64,64);
            
            for LCM_phasemap_index = 1:size(LCM_phasemap_files,2)
                LCM_phasemap_fid = fopen(LCM_phasemap_files(LCM_phasemap_index,'r'));
                LCM_phasemap(LCM_phasemap_index,:,:) = fread(LCM_phasemap_fid, [64,64], 'double');
                fclose(LCM_phasemap_fid);
            end
            
            LCM_phasemap = wrapTo180(LCM_phasemap + LCM_add_phase);
            Sub_LCM_CSI = unwrap_phase_1_0(-CSI_phase + LCM_phasemap);
            Sub_LCM_ImNormal = unwrap_phase_1_0(LCM_phasemap - Im_Normal_Phase);
            
            if(Flip_flag)
                Sub_LCM_ImFlip = unwrap_phase_1_0(LCM_phasemap - Im_Flip_Phase);
                Sub_LCM_ImCorr = unwrap_phase_1_0(LCM_phasemap - Im_Corr_Phase);
            end
            
        end
        
        
        % WEIGHTING
        if(Weighting_flag)
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
            fig1 = figure;
            imagesc(abs(squeeze(CSI(plot_channel,:,:,1,1))));
            title(sprintf('CSI mag %s channel %i', file_struct(read_index).name, plot_channel))
            fig2 = figure;
            imagesc(rad2deg(angle(squeeze(CSI(plot_channel,:,:,1,1)))));
            title(sprintf('CSI pha %s channel %i', file_struct(read_index).name, plot_channel))            

            fig3 = figure;
            imagesc(abs(squeeze(Image_Normal(plot_channel,:,:,1))));
            title(sprintf('ImageNormal mag %s channel %i', file_struct(read_index).name, plot_channel))
            
            fig4 = figure;
            imagesc(rad2deg(angle(squeeze(Image_Normal(plot_channel,:,:,1)))));
            title(sprintf('ImageNormal pha %s channel %i', file_struct(read_index).name, plot_channel))            
           
            if(out_dir_flag)
            	saveas(fig1,sprintf('%s/%s_CSI_mag_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                saveas(fig1,sprintf('%s/%s_CSI_mag_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel))
            	saveas(fig2,sprintf('%s/%s_CSI_pha_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                saveas(fig2,sprintf('%s/%s_CSI_pha_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel))                

            	saveas(fig3,sprintf('%s/%s_ImNormal_mag_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                saveas(fig3,sprintf('%s/%s_ImNormal_mag_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel))
            	saveas(fig4,sprintf('%s/%s_ImNormal_pha_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                saveas(fig4,sprintf('%s/%s_ImNormal_pha_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel))                 
            end  
            
            

            if(Flip_flag)
                fig1 = figure;
                imagesc(abs(squeeze(Image_Flip(plot_channel,:,:,1))));
                title(sprintf('ImageFlip mag %s channel %i', file_struct(read_index).name, plot_channel)) 
                fig2 = figure;
                imagesc(rad2deg(angle((squeeze(Image_Flip(plot_channel,:,:,1))))));
                title(sprintf('ImageFlip pha %s channel %i', file_struct(read_index).name, plot_channel))                
            
                fig3 = figure;
                imagesc(rad2deg(angle(squeeze(Image_Corr(plot_channel,:,:,1,1)))));
                title(sprintf('ImageCorr %s channel %i', file_struct(read_index).name, plot_channel))
                
                if(out_dir_flag)
                    saveas(fig1,sprintf('%s/%s_ImFlip_mag_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                    saveas(fig1,sprintf('%s/%s_ImFlip_mag_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel))
                    saveas(fig2,sprintf('%s/%s_ImFlip_pha_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                    saveas(fig2,sprintf('%s/%s_ImFlip_pha_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel))                

                    saveas(fig3,sprintf('%s/%s_ImCorr_pha_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                    saveas(fig3,sprintf('%s/%s_ImCorr_pha_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel))                 
                end                
                
                

            end
            
            
            
            % PLOT LCM PHAMAP
            
            if(LCM_flag)       
            
                figure;
                imagesc(LCM_phasemap(plot_channel,:,:))
                title(sprintf('LCM phasemap %s channel %i', file_struct(read_index).name, plot_channel))
                
                if(out_dir_flag)
                    saveas(gcf,sprintf('%s/%s_LCM_pha_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                    saveas(gcf,sprintf('%s/%s_LCM_pha_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel))                    
                end

            end     
            
            
            
            
            % PLOT SUB CSI - NORMAL, FLIP, CORR
            fig1 = figure;
            imagesc(squeeze(Sub_CSI_ImNormal(plot_channel,:,:)),[-60 60])
            title(sprintf('Sub CSI - ImageNormal %s channel %i', file_struct(read_index).name, plot_channel))  

            if(out_dir_flag)
                saveas(fig1,sprintf('%s/%s_Sub_CSI_Normal_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                saveas(fig1,sprintf('%s/%s_Sub_CSI_Normal_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel))
            end
            
            
            if(Flip_flag)
                
                fig2 = figure;
                imagesc(Sub_CSI_ImFlip(plot_channel,:,:),[-60 60])
                title(sprintf('Sub CSI - ImageFlip %s channel %i', file_struct(read_index).name, plot_channel)) 

                fig3 = figure;
                imagesc(Sub_CSI_ImCorr(plot_channel,:,:),[-60 60])
                title(sprintf('Sub CSI - ImageCorr %s channel %i', file_struct(read_index).name, plot_channel))    
                
                if(out_dir_flag)                    
                	saveas(fig2,sprintf('%s/%s_Sub_CSI_Flip_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                    saveas(fig2,sprintf('%s/%s_Sub_CSI_Flip_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel))  
                    
                	saveas(fig3,sprintf('%s/%s_Sub_CSI_Corr_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                    saveas(fig3,sprintf('%s/%s_Sub_CSI_Corr_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel))                    
                end
                
                
                
            end
            
            
            
            % PLOT LCM - CORR, CSI
            
            if(LCM_flag)      
            
                fig1 = figure;
                imagesc(Sub_LCM_ImNormal(plot_channel,:,:),[-60 60])
                title(sprintf('SUB LCM - ImageNormal %s channel %i', file_struct(read_index).name, plot_channel))                 
                

                if(Flip_flag)
                    fig2 = figure;
                    imagesc(Sub_LCM_ImFlip(plot_channel,:,:),[-60 60])
                    title(sprintf('SUB LCM - ImageFlip %s channel %i', file_struct(read_index).name, plot_channel))                    
                    
                    fig3 = figure;
                    imagesc(Sub_LCM_ImCorr(plot_channel,:,:),[-60 60])
                    title(sprintf('SUB LCM - ImageCorr %s channel %i', file_struct(read_index).name, plot_channel))
                 	
                end
                
                
                fig4 = figure;
                imagesc(Sub_LCM_CSI(plot_channel,:,:),[-60 60])
                title(sprintf('SUB LCM - CSI %s channel %i', file_struct(read_index).name, plot_channel))                
                

                if(out_dir_flag)                    
                	saveas(fig1,sprintf('%s/%s_Sub_LCM_Image_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                    saveas(fig1,sprintf('%s/%s_Sub_LCM_Image_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel))  
                    
                	saveas(fig4,sprintf('%s/%s_Sub_LCM_CSI_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                    saveas(fig4,sprintf('%s/%s_Sub_LCM_CSI_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel)) 
                    
                    if(Flip_flag)               
                        saveas(fig2,sprintf('%s/%s_Sub_LCM_ImFlip_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                        saveas(fig2,sprintf('%s/%s_Sub_LCM_ImFlip_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel))  

                        saveas(fig3,sprintf('%s/%s_Sub_LCM_ImCorr_channel%02d.fig', out_dir,file_struct(read_index).name,plot_channel))
                        saveas(fig3,sprintf('%s/%s_Sub_LCM_ImCorr_channel%02d.jpg', out_dir,file_struct(read_index).name,plot_channel))    
                    end                   
                    
                    
                    
                end               
                
                
            end             
            
            
            
            
        end
    
        
    
        
        
        % PART II: Statistical computations of subtractions.
        
        if(~strcmpi(file_struct(read_index).mask_position,'NOSTATISTICS'))


            if(strcmpi(file_struct(read_index).mask_position,'NO'))
                % Let user define a mask where the mean and standard deviation of the difference maps are computed.
                imagesc(squeeze(Sub_CSI_ImNormal(plot_channel,:,:)))
                display('Please give me the borders of the mask for statistical computations in voxel units')
                mask_position(1) = input('mask_left_border = ');
                mask_position(2) = input('mask_right_border = ');
                mask_position(3) = input('mask_down_border = ');
                mask_position(4) = input('mask_up_border = ');
            else
                mask_position = file_struct(read_index).mask_position;
            end


            mask = zeros(size(CSI_phase,2),size(CSI_phase,3));
            mask(mask_position(4):mask_position(3),mask_position(1):mask_position(2)) = 1;

            
            Sub_CSI_ImNormal_masked = Sub_CSI_ImNormal_masked;
            Sub_CSI_ImNormal_masked(mask) = -180;
            imagesc(squeeze(Sub_CSI_ImNormal_masked(plot_channel,:,:)),[-60 60])
            
            
            
            % do statistics
            
            Sub_CSI_ImNormal_mean = mean(Sub_CSI_ImNormal(plot_channel,mask))
            Sub_CSI_ImNormal_std = std(Sub_CSI_ImNormal(plot_channel,mask))
                        
            if(Flip_flag)
                Sub_CSI_ImFlip_mean = mean(Sub_CSI_ImFlip(plot_channel,mask))
                Sub_CSI_ImFlip_std = std(Sub_CSI_ImFlip(plot_channel,mask))
                Sub_CSI_ImCorr_mean = mean(Sub_CSI_ImCorr(plot_channel,mask))
                Sub_CSI_ImCorr_std = std(Sub_CSI_ImCorr(plot_channel,mask))
            end

            
            if(LCM_flag)
                Sub_LCM_ImNormal_mean = mean(Sub_LCM_ImNormal(plot_channel,mask))
                Sub_LCM_ImNormal_std = std(Sub_LCM_ImNormal(plot_channel,mask))
                Sub_LCM_CSI_mean = mean(Sub_LCM_CSI(plot_channel,mask))
                Sub_LCM_CSI_std = std(Sub_LCM_CSI(plot_channel,mask))
                
                if(Flip_flag)
                   Sub_LCM_ImFlip_mean = mean(Sub_LCM_ImFlip(plot_channel,mask))
                   Sub_LCM_ImFlip_std = std(Sub_LCM_ImCorr(plot_channel,mask))
                   Sub_LCM_ImCorr_mean = mean(Sub_LCM_ImCorr(plot_channel,mask))
                   Sub_LCM_ImCorr_std = std(Sub_LCM_ImCorr(plot_channel,mask))                 
                end
                
            end
                
            if(out_dir_flag)
                Step1_phasevalidation_txt_fid = fopen(sprintf('%s/Step1_phasevalidation_%s.txt',out_dir,file_struct(read_index).name),'w+');
                fprintf(Step1_phasevalidation_txt_fid,'channel %s\n                   mean       stddev\n\n',plot_channel)          %printf with a field width of 9 and 2 blank characters
                fprintf(Step1_phasevalidation_txt_fid,'Sub_CSI_ImNormal  %9.3f  %9.3f\n',Sub_CSI_ImNormal_mean,Sub_CSI_ImNormal_std)
               
                if(Flip_flag)
                    fprintf(Step1_phasevalidation_txt_fid,'Sub_CSI_ImFlip    %9.3f  %9.3f\n',Sub_CSI_ImFlip_mean,Sub_CSI_ImFlip_std)
                    fprintf(Step1_phasevalidation_txt_fid,'Sub_CSI_ImCorr    %9.3f  %9.3f\n',Sub_CSI_ImCorr_mean,Sub_CSI_ImCorr_std)                    
                end
                
                if(LCM_flag)
                    fprintf(Step1_phasevalidation_txt_fid,'Sub_LCM_ImCSI     %9.3f  %9.3f\n',Sub_LCM_CSI_mean,Sub_LCM_CSI_std)
                    fprintf(Step1_phasevalidation_txt_fid,'Sub_CSI_ImNormal  %9.3f  %9.3f\n',Sub_LCM_ImNormal_mean,Sub_LCM_ImNormal_std)                    
                    
                    if(Flip_flag)
                        fprintf(Step1_phasevalidation_txt_fid,'Sub_LCM_ImFlip    %9.3f  %9.3f\n',Sub_LCM_ImFlip_mean,Sub_LCM_ImFlip_std)
                        fprintf(Step1_phasevalidation_txt_fid,'Sub_CSI_ImCorr    %9.3f  %9.3f\n',Sub_LCM_ImCorr_mean,Sub_LCM_ImCorr_std)
                    end
                    
                end
               
            end
            
            

        end
        
        
        
        
        
        
        
       
        
        % PART II: Spectra
        
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
        
  
        
        
        
        
        
        
        
    end % CSI_file_to_compute
    
    

    
    % save all open figures to save_fig_dir
    
    
    
    
end % if(plotmode)

close all