function read_and_plot_phasemaps_1_4(file_struct,plot_mode,Hamming,LCM_add_phase,AllProbandsComparison_outdir)
%
% file_struct = struct('name', 'CSI', 'Image', 'Flip', 'LCM', 'Weighting', 'water_frequency', 'dwelltime', 'fredir_shift', 'mask_position', 'out_dir')
%
% mask_position = [left_border,right_border,up_border,down_border]
% mask_position = 'NOSTATISTICS' for plots only

%% 0. Definitions, Preparations



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
    
    
    plot_bandwidth_1_0(file_struct)

    
  
    
%% 2. if you want to compare csi with imaging phasemaps and LCM phasemaps
elseif(strcmpi(plot_mode, 'Compare_CSI_Imaging'))
%%

 
Compare_CSI_to_Images_1_1(file_struct,Hamming,LCM_add_phase,AllProbandsComparison_outdir) 
        
        
        
        
        
elseif(strcmpi(plot_mode,'Spectra'))
        
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
                    title(sprintf('CSI unphased %s, channel %d, [x,y]=[%d,%d]', file_struct(read_index).name, plot_channel,plot_x_y(1),plot_x_y(2)), 'Interpreter', 'none')

                    figure;
                    plot(chemshift_vec,CSI_phased(plot_channel,1,1,1,:))
                    set(gca,'XDir','reverse');
                    title(sprintf('CSI phased %s, channel %d, [x,y]=[%d,%d]', file_struct(read_index).name, plot_channel,plot_x_y(1),plot_x_y(2)), 'Interpreter', 'none')               

                end

                figure;
                plot(chemshift_vec,CSI_summed(1,1,1,:))
                set(gca,'XDir','reverse');
                title(sprintf('CSI summed %s, [x,y]=[%d,%d]', file_struct(read_index).name,plot_x_y(1),plot_x_y(2)), 'Interpreter', 'none')   
  
                
            end % Normal_Flip_Corr Phasing
    
    
        end % x_y_plot
  
    
end % if(plotmode)

close all