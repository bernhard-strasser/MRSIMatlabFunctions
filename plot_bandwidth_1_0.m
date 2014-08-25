function plot_bandwidth_1_0(file_struct)


%% 0. DEFINTIONS


out_dir_flag = 0;
if(~strcmpi(file_struct(1).out_dir,'NO'))
    out_dir_flag = 1;
end

[headersize,k_fredir_center,fredir_measured,phadir_measured,total_channel_no] = read_image_dat_meas_header_1_0(file_struct(1).Image);
ROW = phadir_measured; 
COL = phadir_measured; % no bug here, we always have square images
SLC = 1;
Image_Normal = zeros(size(file_struct,2),total_channel_no,ROW,COL,SLC);
Image_Flip = zeros(size(file_struct,2),total_channel_no,ROW,COL,SLC);
Image_Corr = zeros(size(file_struct,2),total_channel_no,ROW,COL,SLC);
if(out_dir_flag)
    out_dir = file_struct(1).out_dir;
    mkdir(out_dir)
end




%% 1. READ DATA


for read_index = 1:size(file_struct,2)

    Image_Normal(read_index,:,:,:,:) = read_image_dat_1_6_Flip_test(file_struct(read_index).Image,128,128,0,0,'AP',0);

    if(Flip_flag)
        Image_Flip(read_index,:,:,:,:) = read_image_dat_1_6(file_struct(read_index).Flip,128,128,0,0,'AP',0);
        Image_Corr(read_index,:,:,:,:) = (Image_Normal(read_index,:,:,:,:) + Image_Flip(read_index,:,:,:,:).*abs(Image_Normal(read_index,:,:,:,:)) ./ abs(Image_Normal(read_index,:,:,:,:)))/2;
    end
end





%% 2. PLOT DATA


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
        title(sprintf('ImageNormal %s channel %i', file_struct(plot_index).name, plot_channel), 'Interpreter', 'none')

        if(out_dir_flag)
            saveas(gcf,sprintf('%s/ImNormal_%s_channel%02d.fig', out_dir,file_struct(plot_index).name,plot_channel))
            saveas(gcf,sprintf('%s/ImNormal_%s_channel%02d.jpg', out_dir,file_struct(plot_index).name,plot_channel))
        end

        if(Flip_flag)
            fig1 = figure;
            imagesc(rad2deg(angle(squeeze(Image_Flip(plot_index,plot_channel,:,:,1)))),[-180 180])
            title(sprintf('ImageFlip %s channel %i', file_struct(plot_index).name, plot_channel), 'Interpreter', 'none')  


            fig2 = figure;
            imagesc(rad2deg(angle(squeeze(Image_Corr(plot_index,plot_channel,:,:,1)))),[-180 180])
            title(sprintf('ImageCorr %s channel %i', file_struct(plot_index).name, plot_channel), 'Interpreter', 'none')

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
        title(sprintf('ImageNormal Sub %s - %s channel %i', file_struct(1).name,file_struct(plot_index).name, plot_channel), 'Interpreter', 'none')  
        if(out_dir_flag)
            saveas(gcf,sprintf('%s/SubNormal_%s_minus_%s_channel%02d.fig', out_dir,file_struct(1).name,file_struct(plot_index).name,plot_channel))
            saveas(gcf,sprintf('%s/SubNormal_%s_minus_%s_channel%02d.jpg', out_dir,file_struct(1).name,file_struct(plot_index).name,plot_channel))
        end          




        if(Flip_flag)

            fig1 = figure;
            imagesc(rad2deg(angle(squeeze(Image_Flip(1,plot_channel,:,:,1)))) - rad2deg(angle(squeeze(Image_Flip(plot_index,plot_channel,:,:,1)))),[-15 15])
            title(sprintf('ImageFlip Sub %s - %s channel %i', file_struct(1).name,file_struct(plot_index).name, plot_channel), 'Interpreter', 'none')



            fig2 = figure;
            imagesc(rad2deg(angle(squeeze(Image_Corr(1,plot_channel,:,:,1)))) - rad2deg(angle(squeeze(Image_Corr(plot_index,plot_channel,:,:,1)))),[-15 15])
            title(sprintf('ImageCorr Sub %s - %s channel %i', file_struct(1).name,file_struct(plot_index).name, plot_channel), 'Interpreter', 'none')
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