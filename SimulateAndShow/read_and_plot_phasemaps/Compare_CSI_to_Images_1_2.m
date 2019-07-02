function Compare_CSI_to_Images_1_2(file_struct,Hamming,AllProbandsComparison_outdir)


for read_index = 1:size(file_struct,2)  % CSI DATA = LARGE --> COMPUTATIONS FOR EACH FILE INDIVIDUALLY
    
    
   %% 0. DEFINITIONS, PREPARATIONS
   
    % 0.1 FLAGS
    Flip_flag = 0;
    LCM_flag = 0;
    Weighting_flag = 0;
    out_dir_flag = 0;
    if(~strcmpi(file_struct(read_index).Flip,'NO'))
        Flip_flag = 1;
    end
    if(~strcmpi(file_struct(read_index).LCM,'NO'))
        LCM_flag = 1;
    end
    if(~strcmpi(file_struct(read_index).Weighting,'NO'))
        Weighting_flag = 1;
    end
    if(~strcmpi(file_struct(read_index).out_dir,'NO'))
        out_dir_flag = 1;
    end

    Phantom_flag = file_struct(read_index).Phantom_flag

    % 0.2 mkdir 
    if(out_dir_flag)
        out_dir = file_struct(read_index).out_dir;
        mkdir(out_dir)
    end
    display([char(10),'Processing ',file_struct(read_index).name, char(10)])  


    % 0.3 DEFINITIONS
    [ROW,COL,SLC,vecSize,total_channel_no,total_k_points,headersize] = read_csi_dat_meas_header_1_0(file_struct(read_index).CSI);
    LCM_add_phase = file_struct(read_index).LCM_add_phase;
    LCM_add_phase(strcmpi('NO',LCM_add_phase)) = 0;
    
    %% 1. READ IN DATA

    % 1.1 CSI
    CSI = read_csi_dat_1_3(file_struct(read_index).CSI);
    
    
    % 1.2 Image NORMAL
    Image_Normal = read_image_dat_1_6(file_struct(read_index).Image,ROW,COL,0,file_struct(read_index).fredir_shift,'AP',0);


    % 1.3 Image FLIP, COMPUTE CORR
    if(Flip_flag)
        Image_Flip = read_image_dat_1_6(file_struct(read_index).Flip,ROW,COL,0,file_struct(read_index).fredir_shift,'AP',1);
    end 

    
    % 1.4 LCM phasemap
    if(LCM_flag)

        LCM_phasemap_files = file_struct(read_index).LCM;
        LCM_phasemap = zeros(size(LCM_phasemap_files,2),ROW,COL);

        for LCM_phasemap_index = 1:size(LCM_phasemap_files,2)
            if(~strcmpi(LCM_phasemap_files{LCM_phasemap_index}, 'NO'))
                LCM_phasemap_fid = fopen(LCM_phasemap_files{LCM_phasemap_index},'r');
                LCM_phasemap(LCM_phasemap_index,:,:) = fread(LCM_phasemap_fid, [ROW,COL], 'float');
                fclose(LCM_phasemap_fid);
            end
        end
    end


    % 1.5 WEIGHTING
    if(Weighting_flag)
        Image_Weighting = abs(read_image_dat_1_5(file_struct(read_index).Weighting));
    else
        Image_Weighting = sum(abs(Image_Normal),1);
    end        






    
    
    
    %% 2. MAKE CORRECTIONS, COMPUTE PHASES & MAGNITUDES

    % 2.1 CSI
    if(strcmpi(Hamming,'hamming'))
        CSI = Hamming_filtering_1_1(CSI,[2 3]);
    end
    CSI = conj(CSI);
    CSI_phase = rad2deg(angle(reshape(squeeze(CSI(:,:,:,1,1)), [total_channel_no,ROW,COL])));   %reshape because if just 1 channel --> this dimension gets squeezed also!
    CSI_mag = abs(reshape(squeeze(CSI(:,:,:,1,1)), [total_channel_no,ROW,COL]));

    
    if(Phantom_flag)
        SubPhase = deg2rad(63.2898);
    else
        SubPhase = deg2rad(30.0035);
    end   
    
    % 2.2 ImNormal
    Image_Normal = circshift(Image_Normal, [0 1 0]);   % 0 1 -1
    Image_Normal = Image_Normal * exp(-1i*SubPhase);
    ImNormal_Phase = rad2deg(angle(reshape(squeeze(Image_Normal(:,:,:,1)),[total_channel_no,ROW,COL])));
    ImNormal_mag = abs(reshape(squeeze(Image_Normal(:,:,:,1)),[total_channel_no,ROW,COL]));
        
    
    % 2.3 ImFlip & ImCorr
    if(Flip_flag)
        Image_Flip = circshift(Image_Flip, [0 0 -1]);
        Image_Flip = Image_Flip * exp(-1i*SubPhase);
        %Image_Corr = ( Image_Normal + Image_Flip )/2;       
        Image_Corr = ( Image_Normal + Image_Flip .* abs(Image_Normal)./abs(Image_Flip))/2;
        ImFlip_Phase = rad2deg(angle(reshape(squeeze(Image_Flip(:,:,:,1)),[total_channel_no,ROW,COL])));
        ImFlip_mag = abs(reshape(squeeze(Image_Flip(:,:,:,1)),[total_channel_no,ROW,COL]));
        ImCorr_Phase = rad2deg(angle(reshape(squeeze(Image_Corr(:,:,:,1)),[total_channel_no,ROW,COL])));
    end 
    


    
    LCM_add_phase = 0;
    
    
    % 2.4 LCM phasemap
    if(LCM_flag)
        LCM_phasemap = wrapTo180(LCM_phasemap + LCM_add_phase);  
    end
    
    
    
    
    
    
    %% 3. Subtraction mapos
    
    % 3.1 CSI AS REFERENCE
    
    Sub_CSI_ImNormal = CSI_phase - ImNormal_Phase; 
    Sub_CSI_ImNormal = wrapTo180(Sub_CSI_ImNormal);
    
    if(Flip_flag)
        Sub_CSI_ImFlip = CSI_phase - ImFlip_Phase;
        Sub_CSI_ImCorr = CSI_phase - ImCorr_Phase;
        Sub_CSI_ImFlip = wrapTo180(Sub_CSI_ImFlip);
        Sub_CSI_ImCorr = wrapTo180(Sub_CSI_ImCorr);
    end
    

    
    % 3.2 LCM AS REFERENCE
    
    if(LCM_flag)
        Sub_LCM_CSI = wrapTo180(LCM_phasemap - CSI_phase);
        Sub_LCM_ImNormal = wrapTo180(LCM_phasemap - ImNormal_Phase);
        
        if(Flip_flag)
            Sub_LCM_ImFlip = wrapTo180(LCM_phasemap - ImFlip_Phase);
            Sub_LCM_ImCorr = wrapTo180(LCM_phasemap - ImCorr_Phase);
        end
    end
    
    
    

%% 4. Channel LOOP

    plot_channel = 1;
    while(plot_channel ~= 666)

    display( [char(10) 'Enter channel you want to plot. Enter 666 for continueing with further program' char(10)] )
    plot_channel = input('channel = ');
    close all;

    if(plot_channel == 666)
        continue
    end           



        %% 6.1 PLOT MAG AND PHASEMAPS


        % 6.1.1 magnitudes

        fig_CSI_mag = figure;
        imagesc(squeeze(CSI_mag(plot_channel,:,:)));
        title(sprintf('CSI mag %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none')
        colorbar

        fig_ImNormal_mag = figure;
        imagesc(squeeze(ImNormal_mag(plot_channel,:,:)));
        title(sprintf('ImageNormal mag %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none')  
        colorbar

        if(Flip_flag)
            fig_ImFlip_mag = figure;
            imagesc(squeeze(ImFlip_mag(plot_channel,:,:)));
            title(sprintf('ImageFlip mag %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none') 
            colorbar
        end



        % 6.1.2 phasemaps

        fig_CSI_pha = figure;
        imagesc(squeeze(CSI_phase(plot_channel,:,:)));
        title(sprintf('CSI pha %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none')  
        colorbar

        fig_ImNormal_pha = figure;
        imagesc(squeeze(ImNormal_Phase(plot_channel,:,:)));
        title(sprintf('ImageNormal pha %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none')
        colorbar

        if(Flip_flag)
            fig_ImFlip_pha = figure;
            imagesc(squeeze(ImFlip_Phase(plot_channel,:,:)));
            title(sprintf('ImageFlip pha %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none')
            colorbar                

            fig_ImCorr_pha = figure;
            imagesc(squeeze(ImCorr_Phase(plot_channel,:,:)));
            title(sprintf('ImageCorr %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none')
            colorbar            
        end

        if(LCM_flag)
            fig_LCM = figure;
            imagesc(squeeze(LCM_phasemap(plot_channel,:,:)))
            title(sprintf('LCM phasemap %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none')
            colorbar            
        end



        %% 6.2 PLOT SUBTRACTION MAPS



%         % 6.2.1 CSI AS REFERENCE
        fig_Sub_CSI_ImNormal = figure;
        imagesc(squeeze(Sub_CSI_ImNormal(plot_channel,:,:)),[-60 60])
        title(sprintf('Sub CSI - ImageNormal %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none')
        colorbar         


        if(Flip_flag)
            fig_Sub_CSI_ImFlip = figure;
            imagesc(squeeze(Sub_CSI_ImFlip(plot_channel,:,:)),[-60 60])
            title(sprintf('Sub CSI - ImageFlip %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none')
            colorbar 

            fig_Sub_CSI_ImCorr = figure;
            imagesc(squeeze(Sub_CSI_ImCorr(plot_channel,:,:)),[-60 60])
            title(sprintf('Sub CSI - ImageCorr %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none')
            colorbar           
        end




        % 6.2.2 LCM AS REFERENCE
        if(LCM_flag)
            % CSI
            fig_Sub_LCM_CSI = figure;
            imagesc(squeeze(Sub_LCM_CSI(plot_channel,:,:)),[-60 60])
            title(sprintf('SUB LCM - CSI %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none')
            colorbar  

            % ImNormal
            fig_Sub_LCM_ImNormal = figure;
            imagesc(squeeze(Sub_LCM_ImNormal(plot_channel,:,:)),[-60 60])
            title(sprintf('SUB LCM - ImageNormal %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none')
            colorbar 

            if(Flip_flag)
                % ImFlip
                fig_Sub_LCM_ImFlip = figure;
                imagesc(squeeze(Sub_LCM_ImFlip(plot_channel,:,:)),[-60 60])
                title(sprintf('SUB LCM - ImageFlip %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none')
                colorbar                    

                % ImCorr
                fig_Sub_LCM_ImCorr = figure;
                imagesc(squeeze(Sub_LCM_ImCorr(plot_channel,:,:)),[-60 60])
                title(sprintf('SUB LCM - ImageCorr %s channel %i', file_struct(read_index).name, plot_channel), 'Interpreter', 'none')
                colorbar
            end                 
        end












    %% 6.3 Statistical computations



        if(~strcmpi(file_struct(read_index).mask_position,'NOSTATISTICS'))


            % 6.3.1 mask input
            if(strcmpi(file_struct(read_index).mask_position,'NO'))
                % Let user define a mask where the mean and standard deviation of the difference maps are computed.
                imagesc(squeeze(Sub_CSI_ImNormal(plot_channel,:,:)))
                display('Please give me the borders of the mask for statistical computations in voxel units')
                mask_position(1) = input('mask_left_border = ');
                mask_position(2) = input('mask_right_border = ');
                mask_position(3) = input('mask_up_border = ');
                mask_position(4) = input('mask_down_border = ');
            else
                mask_position = file_struct(read_index).mask_position;
            end

            mask = zeros(size(CSI_phase,2),size(CSI_phase,3));
            mask(mask_position(3):mask_position(4),mask_position(1):mask_position(2)) = 1;
            mask = logical(mask);
            size(mask)




            % 6.3.2 plot mask
            Sub_CSI_ImNormal_masked = Sub_CSI_ImNormal;
            Sub_CSI_ImNormal_masked(:,mask) = -181;
            fig_mask = figure;
            imagesc(squeeze(Sub_CSI_ImNormal_masked(plot_channel,:,:)),[-180 180])
            colorbar




            % 6.3.3 compute mean, std

            % Initialization
            Sub_CSI_ImFlip_mean=999; Sub_CSI_ImFlip_std=999; Sub_CSI_ImCorr_mean=999; Sub_CSI_ImCorr_std=999; Sub_LCM_ImNormal_mean=999; Sub_LCM_ImNormal_std=999; Sub_LCM_CSI_mean=999; Sub_LCM_CSI_std=999;
            Sub_LCM_ImFlip_mean=999; Sub_LCM_ImFlip_std=999; Sub_LCM_ImCorr_mean=999; Sub_LCM_ImCorr_std=999;



            Sub_CSI_ImNormal_mean = mean(Sub_CSI_ImNormal(plot_channel,mask))
            Sub_CSI_ImNormal_std = std(Sub_CSI_ImNormal(plot_channel,mask))
            Sub_CSI_ImNormal_rss = sqrt(sum(Sub_CSI_ImNormal(plot_channel,mask).^2))


            if(Flip_flag)
                Sub_CSI_ImFlip_mean = mean(Sub_CSI_ImFlip(plot_channel,mask))
                Sub_CSI_ImFlip_std = std(Sub_CSI_ImFlip(plot_channel,mask))
                Sub_CSI_ImFlip_rss = sqrt(sum(Sub_CSI_ImFlip(plot_channel,mask).^2))              
                Sub_CSI_ImCorr_mean = mean(Sub_CSI_ImCorr(plot_channel,mask))
                Sub_CSI_ImCorr_std = std(Sub_CSI_ImCorr(plot_channel,mask))
                Sub_CSI_ImCorr_rss = sqrt(sum(Sub_CSI_ImCorr(plot_channel,mask).^2))                            
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
                end %Flip_flag
            end % LCM_Flag
        end %if(NOSTATISTICS)






         %% 6.4 write (to) files

        if(out_dir_flag)
                        
            % 6.4.1 magnitudes
            saveas(fig_CSI_mag,sprintf('%s/channel%02d_%s_CSI_mag.fig', out_dir,plot_channel,file_struct(read_index).name))
            saveas(fig_CSI_mag,sprintf('%s/channel%02d_%s_CSI_mag.jpg', out_dir,plot_channel,file_struct(read_index).name))
            saveas(fig_ImNormal_mag,sprintf('%s/channel%02d_%s_ImNormal_mag.fig', out_dir,plot_channel,file_struct(read_index).name))
            saveas(fig_ImNormal_mag,sprintf('%s/channel%02d_%s_ImNormal_mag.jpg', out_dir,plot_channel,file_struct(read_index).name)) 

            if(Flip_flag)
                saveas(fig_ImFlip_mag,sprintf('%s/channel%02d_%s_ImFlip_mag.fig', out_dir,plot_channel,file_struct(read_index).name))
                saveas(fig_ImFlip_mag,sprintf('%s/channel%02d_%s_ImFlip_mag.jpg', out_dir,plot_channel,file_struct(read_index).name)) 
            end




            % 6.4.2 phasemaps
            % CSI
            saveas(fig_CSI_pha,sprintf('%s/channel%02d_%s_CSI_pha.fig', out_dir,plot_channel,file_struct(read_index).name))
            saveas(fig_CSI_pha,sprintf('%s/channel%02d_%s_CSI_pha.jpg', out_dir,plot_channel,file_struct(read_index).name))  
            % ImNormal
            saveas(fig_ImNormal_pha,sprintf('%s/channel%02d_%s_ImNormal_pha.fig', out_dir,plot_channel,file_struct(read_index).name))
            saveas(fig_ImNormal_pha,sprintf('%s/channel%02d_%s_ImNormal_pha.jpg', out_dir,plot_channel,file_struct(read_index).name)) 

            if(Flip_flag)
                % ImFlip
                saveas(fig_ImFlip_pha,sprintf('%s/channel%02d_%s_ImFlip_pha.fig', out_dir,plot_channel,file_struct(read_index).name))
                saveas(fig_ImFlip_pha,sprintf('%s/channel%02d_%s_ImFlip_pha.jpg', out_dir,plot_channel,file_struct(read_index).name)) 
                % ImCorr
                saveas(fig_ImCorr_pha,sprintf('%s/channel%02d_%s_ImCorr_pha.fig', out_dir,plot_channel,file_struct(read_index).name))
                saveas(fig_ImCorr_pha,sprintf('%s/channel%02d_%s_ImCorr_pha.jpg', out_dir,plot_channel,file_struct(read_index).name)) 
            end

            if(LCM_flag)
                % LCM
                saveas(fig_LCM,sprintf('%s/channel%02d_%s_LCM_pha.fig', out_dir,plot_channel,file_struct(read_index).name))
                saveas(fig_LCM,sprintf('%s/channel%02d_%s_LCM_pha.jpg', out_dir,plot_channel,file_struct(read_index).name))
            end





            
            % 6.4.3 Sub CSI as reference
            % ImNormal
            saveas(fig_Sub_CSI_ImNormal,sprintf('%s/channel%02d_%s_Sub_CSI_Normal.fig', out_dir,plot_channel,file_struct(read_index).name))
            saveas(fig_Sub_CSI_ImNormal,sprintf('%s/channel%02d_%s_Sub_CSI_Normal.jpg', out_dir,plot_channel,file_struct(read_index).name))

            if(Flip_flag)
                % ImFlip
                saveas(fig_Sub_CSI_ImFlip,sprintf('%s/channel%02d_%s_Sub_CSI_Flip.fig', out_dir,plot_channel,file_struct(read_index).name))
                saveas(fig_Sub_CSI_ImFlip,sprintf('%s/channel%02d_%s_Sub_CSI_Flip.jpg', out_dir,plot_channel,file_struct(read_index).name))  
                % ImCorr
                saveas(fig_Sub_CSI_ImCorr,sprintf('%s/channel%02d_%s_Sub_CSI_Corr.fig', out_dir,plot_channel,file_struct(read_index).name))
                saveas(fig_Sub_CSI_ImCorr,sprintf('%s/channel%02d_%s_Sub_CSI_Corr.jpg', out_dir,plot_channel,file_struct(read_index).name))                    
            end





            % 6.4.4 Sub LCM as reference
            % CSI
            if(LCM_flag)
                saveas(fig_Sub_LCM_CSI,sprintf('%s/channel%02d_%s_Sub_LCM_CSI.fig', out_dir,plot_channel,file_struct(read_index).name))
                saveas(fig_Sub_LCM_CSI,sprintf('%s/channel%02d_%s_Sub_LCM_CSI.jpg', out_dir,plot_channel,file_struct(read_index).name))  
                % ImNormal
                saveas(fig_Sub_LCM_ImNormal,sprintf('%s/channel%02d_%s_Sub_LCM_ImNormal.fig', out_dir,plot_channel,file_struct(read_index).name))
                saveas(fig_Sub_LCM_ImNormal,sprintf('%s/channel%02d_%s_Sub_LCM_ImNormal.jpg', out_dir,plot_channel,file_struct(read_index).name))  

                % ImFlip & ImCorr
                if(Flip_flag) 
                    saveas(fig_Sub_LCM_ImFlip,sprintf('%s/channel%02d_%s_Sub_LCM_ImFlip.fig', out_dir,plot_channel,file_struct(read_index).name))
                    saveas(fig_Sub_LCM_ImFlip,sprintf('%s/channel%02d_%s_Sub_LCM_ImFlip.jpg', out_dir,plot_channel,file_struct(read_index).name))  
                    saveas(fig_Sub_LCM_ImCorr,sprintf('%s/channel%02d_%s_Sub_LCM_ImCorr.fig', out_dir,plot_channel,file_struct(read_index).name))
                    saveas(fig_Sub_LCM_ImCorr,sprintf('%s/channel%02d_%s_Sub_LCM_ImCorr.jpg', out_dir,plot_channel,file_struct(read_index).name))    
                end %Flip_flag
            end %LCM_flag


            % mask
            saveas(fig_mask,sprintf('%s/%s_mask.fig', out_dir,file_struct(read_index).name))
            saveas(fig_mask,sprintf('%s/%s_mask.jpg', out_dir,file_struct(read_index).name))         



            % 6.4.5 write statistics to file
            if(~strcmpi(file_struct(read_index).mask_position,'NOSTATISTICS'))

                phasevalidation_txt_fid = fopen(sprintf('%s/channel%02d_phasevalidation_%s.txt',out_dir,plot_channel,file_struct(read_index).name),'a');

                fprintf(phasevalidation_txt_fid,'channel %s\n                   mean         stddev\n',num2str(plot_channel));          %printf with a field width of 9 and 2 blank characters
                fprintf(phasevalidation_txt_fid,'Sub_CSI_ImNormal  %9.4f  %9.4f\n',Sub_CSI_ImNormal_mean,Sub_CSI_ImNormal_std);

                fprintf(phasevalidation_txt_fid,'Sub_CSI_ImFlip    %9.4f  %9.4f\n',Sub_CSI_ImFlip_mean,Sub_CSI_ImFlip_std);
                fprintf(phasevalidation_txt_fid,'Sub_CSI_ImCorr    %9.4f  %9.4f\n',Sub_CSI_ImCorr_mean,Sub_CSI_ImCorr_std);                    

                fprintf(phasevalidation_txt_fid,'Sub_LCM_CSI       %9.4f  %9.4f\n',Sub_LCM_CSI_mean,Sub_LCM_CSI_std);
                fprintf(phasevalidation_txt_fid,'Sub_LCM_ImNormal  %9.4f  %9.4f\n',Sub_LCM_ImNormal_mean,Sub_LCM_ImNormal_std);                    

                fprintf(phasevalidation_txt_fid,'Sub_LCM_ImFlip    %9.4f  %9.4f\n',Sub_LCM_ImFlip_mean,Sub_LCM_ImFlip_std);
                fprintf(phasevalidation_txt_fid,'Sub_LCM_ImCorr    %9.4f  %9.4f\n\n',Sub_LCM_ImCorr_mean,Sub_LCM_ImCorr_std);

                % If the statistics file do not exist create it and the first line
                % LCM file
                if(~logical(exist(sprintf('%s/statistics_AllProbandsComparison_LCM.txt',AllProbandsComparison_outdir),'file')))
                    AllProbandsComparison_LCM_fid = fopen(sprintf('%s/statistics_AllProbandsComparison_LCM.txt',AllProbandsComparison_outdir),'w');
                    fprintf(AllProbandsComparison_LCM_fid,'Proband_Channel        \tmean_Sub_LCM-ImNor\tstdd_Sub_LCM-ImNor\tmean_Sub_LCM-ImFli\tstdd_Sub_LCM-ImFli\tmean_Sub_LCM-ImCor\tstdd_Sub_LCM-ImCor\n');        
                else
                    AllProbandsComparison_LCM_fid = fopen(sprintf('%s/statistics_AllProbandsComparison_LCM.txt',AllProbandsComparison_outdir),'a');
                end
                % write difference of CSI - Image and LCM - Image to file
                fprintf(AllProbandsComparison_LCM_fid,'%10.10s_channel%02d\t%18.4f\t%18.4f\t%18.4f\t%18.4f\t%18.4f\t%18.4f\n', ...
                file_struct(read_index).name,plot_channel,Sub_LCM_ImNormal_mean,Sub_LCM_ImNormal_std,Sub_LCM_ImFlip_mean,Sub_LCM_ImFlip_std,Sub_LCM_ImCorr_mean,Sub_LCM_ImCorr_std);                                
               
            
                
                % CSI file
                if(~logical(exist(sprintf('%s/statistics_AllProbandsComparison_CSI.txt',AllProbandsComparison_outdir),'file')))
                    AllProbandsComparison_CSI_fid = fopen(sprintf('%s/statistics_AllProbandsComparison_CSI.txt',AllProbandsComparison_outdir),'w');
                    fprintf(AllProbandsComparison_CSI_fid,'Proband_Channel        \tmean_Sub_CSI-ImNor\tstdd_Sub_CSI-ImNor\tmean_Sub_CSI-ImFli\tstdd_Sub_CSI-ImFli\tmean_Sub_CSI-ImCor\tstdd_Sub_CSI-ImCor\n');        
                else
                    AllProbandsComparison_CSI_fid = fopen(sprintf('%s/statistics_AllProbandsComparison_CSI.txt',AllProbandsComparison_outdir),'a');
                end
                % write difference of CSI - Image and CSI - Image to file
                fprintf(AllProbandsComparison_CSI_fid,'%10.10s_channel%02d\t%18.4f\t%18.4f\t%18.4f\t%18.4f\t%18.4f\t%18.4f\n', ...
                file_struct(read_index).name,plot_channel,Sub_CSI_ImNormal_mean,Sub_CSI_ImNormal_std,Sub_CSI_ImFlip_mean,Sub_CSI_ImFlip_std,Sub_CSI_ImCorr_mean,Sub_CSI_ImCorr_std);               

            
                

                fclose(phasevalidation_txt_fid);
                fclose(AllProbandsComparison_LCM_fid);
                fclose(AllProbandsComparison_CSI_fid);               
            end
 
        end %out_dir_flag 
            
    end %while(plot_channel)

end %read_index   
        
