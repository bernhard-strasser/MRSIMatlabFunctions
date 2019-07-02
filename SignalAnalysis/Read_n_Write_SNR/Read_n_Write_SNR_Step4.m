function Read_n_Write_SNR_Step4(file_struct,SNR_textfile)

size(file_struct)
%for read_index = 1:size(file_struct,2)
 
for read_index = 1:size(file_struct,2)  
   %% 0. DEFINITIONS, PREPARATIONS
   


    %NameWritten_flag = false;

    % 0.2 mkdir 
    if(~exist(file_struct(read_index).OutDir,'dir'))
        out_dir = file_struct(read_index).OutDir;
        mkdir(out_dir)
    end
    display([char(10),'Processing ',file_struct(read_index).name, char(10)])  

    ROW = file_struct(read_index).MatrixSize(1);
    COL = file_struct(read_index).MatrixSize(2);
    SLC = 1;
    
    
    
    
    %% 1. READ IN DATA

    % 1.1 masks
    
    if(ischar(file_struct(read_index).Mask_Siem))
        mask = read_RawFiles_1_0(file_struct(read_index).Mask_Siem,'float',ROW,COL,1);
    else
        mask = file_struct(read_index).Mask_Siem;
    end
    
    if(ischar(file_struct(read_index).Mask_OurMeth)) 
        mask_OurMeth = read_RawFiles_1_0(file_struct(read_index).Mask_OurMeth,'float',ROW,COL,1);    
    else
        mask_OurMeth = file_struct(read_index).Mask_Siem;
    end    
    
    % 1.2 LCM SNR maps
    
    LCM_Siem_SNR = read_RawFiles_1_0(file_struct(read_index).LCM_Siem,'float',ROW,COL,1);
    LCM_OurMeth_SNR = read_RawFiles_1_0(file_struct(read_index).LCM_OurMeth,'float',ROW,COL,1);  
    
    
    % 1.3 MATLAB SNR maps
    
    MATLAB_Siem_SNR = read_RawFiles_1_0(file_struct(read_index).MATLAB_Siem,'float',ROW,COL,1);
    MATLAB_OurMeth_SNR = read_RawFiles_1_0(file_struct(read_index).MATLAB_OurMeth,'float',ROW,COL,1);    
            

%     if(strcmp(file_struct(read_index).Mask_Siem,file_struct(read_index).Mask_OurMeth))
        LCM_Ratio_SNR = LCM_OurMeth_SNR ./ LCM_Siem_SNR;
        MATLAB_Ratio_SNR = MATLAB_OurMeth_SNR ./ MATLAB_Siem_SNR;
%     else
%         LCM_Ratio_SNR = NaN * ones([ROW,COL]);
%         MATLAB_Ratio_SNR = NaN * ones([ROW,COL]);
%         LCM_Siem_SNR(~mask) = NaN;
%         MATLAB_Siem_SNR(~mask) = NaN;
%         LCM_OurMeth_SNR(~mask_OurMeth) = NaN;
%         MATLAB_OurMeth_SNR(~mask_OurMeth) = NaN;
%     end

    
    

    %% 2. PLOT MAPS


    % 6.1.1 masks

    fig_mask = figure('Visible','Off');
    imagesc(mask,[0,1]);
    title(sprintf('mask %s', file_struct(read_index).name), 'Interpreter', 'none')

    if(~strcmp(file_struct(read_index).Mask_Siem,file_struct(read_index).Mask_OurMeth))
        fig_mask_OurMeth = figure('Visible','Off');
        imagesc(mask_OurMeth,[0,1]);
        title(sprintf('mask OurMeth %s', file_struct(read_index).name), 'Interpreter', 'none')
    end


    if(eq(ROW,32))
        lim_LCM_OurMethSiem = 85;
        lim_MATLAB_OurMethSiem = 1500;   
        lim_LCM_Ratio = 9;
        lim_MATLAB_Ratio = 9;
    else
        lim_LCM_OurMethSiem = 25;
        lim_MATLAB_OurMethSiem = 110;  
        lim_LCM_Ratio = 7;
        lim_MATLAB_Ratio = 7;            
    end




    fig_LCM_Siem = figure('Visible','Off');
    imagesc(LCM_Siem_SNR, [0 lim_LCM_OurMethSiem]);
    title(sprintf('LCM Siem SNR map %s', file_struct(read_index).name), 'Interpreter', 'none')  
    colorbar('FontSize',30)

    fig_LCM_OurMeth = figure('Visible','Off');
    imagesc(LCM_OurMeth_SNR, [0 lim_LCM_OurMethSiem]);
    title(sprintf('LCM OurMeth SNR map %s', file_struct(read_index).name), 'Interpreter', 'none')  
    colorbar('FontSize',30)

    fig_MATLAB_Siem = figure('Visible','Off');
    imagesc(MATLAB_Siem_SNR,[0 lim_MATLAB_OurMethSiem]);
    title(sprintf('MATLAB Siem SNR map %s', file_struct(read_index).name), 'Interpreter', 'none')  
    colorbar('FontSize',30)

    fig_MATLAB_OurMeth = figure('Visible','Off');
    imagesc(MATLAB_OurMeth_SNR, [0 lim_MATLAB_OurMethSiem]);
    title(sprintf('MATLAB OurMeth SNR map %s', file_struct(read_index).name), 'Interpreter', 'none')  
    colorbar('FontSize',30)

    fig_LCM_Ratio = figure('Visible','Off');
    imagesc(LCM_Ratio_SNR, [0 lim_LCM_Ratio]);
    title(sprintf('MATLAB Ratio OurMeth/Siem SNR map %s', file_struct(read_index).name), 'Interpreter', 'none')  
    colorbar('FontSize',30)


    fig_MATLAB_Ratio = figure('Visible','Off');
    imagesc(MATLAB_Ratio_SNR, [0 lim_MATLAB_Ratio]);
    title(sprintf('MATLAB Ratio OurMeth/Siem SNR map %s', file_struct(read_index).name), 'Interpreter', 'none')  
    colorbar('FontSize',30)  


    saveas(fig_mask,sprintf('%s/mask',file_struct(read_index).OutDir),'epsc2')
    saveas(fig_mask,sprintf('%s/mask',file_struct(read_index).OutDir),'fig')
    saveas(fig_mask,sprintf('%s/mask',file_struct(read_index).OutDir),'jpg')   
    if(~strcmp(file_struct(read_index).Mask_Siem,file_struct(read_index).Mask_OurMeth))
        saveas(fig_mask_OurMeth,sprintf('%s/mask_OurMeth',file_struct(read_index).OutDir),'epsc2')
        saveas(fig_mask_OurMeth,sprintf('%s/mask_OurMeth',file_struct(read_index).OutDir),'fig')
        saveas(fig_mask_OurMeth,sprintf('%s/mask_OurMeth',file_struct(read_index).OutDir),'jpg')        
    end    



    saveas(fig_LCM_Siem,sprintf('%s/SNR_LCM_Siem',file_struct(read_index).OutDir),'epsc2')
    saveas(fig_LCM_Siem,sprintf('%s/SNR_LCM_Siem',file_struct(read_index).OutDir),'fig')
    saveas(fig_LCM_Siem,sprintf('%s/SNR_LCM_Siem',file_struct(read_index).OutDir),'jpg')    

    saveas(fig_LCM_OurMeth,sprintf('%s/SNR_LCM_OurMeth',file_struct(read_index).OutDir),'epsc2')
    saveas(fig_LCM_OurMeth,sprintf('%s/SNR_LCM_OurMeth',file_struct(read_index).OutDir),'fig')
    saveas(fig_LCM_OurMeth,sprintf('%s/SNR_LCM_OurMeth',file_struct(read_index).OutDir),'jpg')
    
    saveas(fig_LCM_Ratio,sprintf('%s/SNR_LCM_Ratio',file_struct(read_index).OutDir),'epsc2')
    saveas(fig_LCM_Ratio,sprintf('%s/SNR_LCM_Ratio',file_struct(read_index).OutDir),'fig')        
    saveas(fig_LCM_Ratio,sprintf('%s/SNR_LCM_Ratio',file_struct(read_index).OutDir),'jpg')        


    saveas(fig_MATLAB_Siem,sprintf('%s/SNR_MATLAB_Siem',file_struct(read_index).OutDir),'epsc2')
    saveas(fig_MATLAB_Siem,sprintf('%s/SNR_MATLAB_Siem',file_struct(read_index).OutDir),'fig')
    saveas(fig_MATLAB_Siem,sprintf('%s/SNR_MATLAB_Siem',file_struct(read_index).OutDir),'jpg')
    
    saveas(fig_MATLAB_OurMeth,sprintf('%s/SNR_MATLAB_OurMeth',file_struct(read_index).OutDir),'epsc2')
    saveas(fig_MATLAB_OurMeth,sprintf('%s/SNR_MATLAB_OurMeth',file_struct(read_index).OutDir),'fig')
    saveas(fig_MATLAB_OurMeth,sprintf('%s/SNR_MATLAB_OurMeth',file_struct(read_index).OutDir),'jpg')
    
    saveas(fig_MATLAB_Ratio,sprintf('%s/SNR_MATLAB_Ratio',file_struct(read_index).OutDir),'epsc2')
    saveas(fig_MATLAB_Ratio,sprintf('%s/SNR_MATLAB_Ratio',file_struct(read_index).OutDir),'fig')
    saveas(fig_MATLAB_Ratio,sprintf('%s/SNR_MATLAB_Ratio',file_struct(read_index).OutDir),'jpg')        
        
    close all;
        
              

            %% 6.4.5 write statistics to file


                % If the statistics file do not exist create it and the first line
                % LCM file
                if(~exist(SNR_textfile,'file'))
                    SNR_textfile_fid = fopen(SNR_textfile,'w');
                    fprintf(SNR_textfile_fid,'SNR of %s',file_struct(read_index).name);      
                else
                    SNR_textfile_fid = fopen(SNR_textfile,'a');
                    fprintf(SNR_textfile_fid,'\n\n\nSNR of %s',file_struct(read_index).name);
                end
                
                fprintf(SNR_textfile_fid,'\nVoxel (ROW,COL,SLC)\t\tLCM_Siem\t\tLCM_OurMeth\t\tLCM_Ratio\t\tMAT_Siem\t\tMAT_OurMeth\t\tMAT_Ratio\n');                  
                
                
                for x = 1:ROW
                    for y = 1:COL
                        for z = 1:SLC
                            if(mask(x,y,z) == 0 && mask_OurMeth(x,y,z) == 0)
                                continue;
                            end
                            
                            fprintf(SNR_textfile_fid,'(%d,%d,%d)\t\t\t\t%7.2f\t\t%7.2f\t\t%7.3f\t\t%7.2f\t\t%7.2f\t\t%7.3f\n', ...
                            x,y,z,LCM_Siem_SNR(x,y),LCM_OurMeth_SNR(x,y),LCM_Ratio_SNR(x,y),MATLAB_Siem_SNR(x,y),MATLAB_OurMeth_SNR(x,y),MATLAB_Ratio_SNR(x,y));  
                        
                        end
                    end
                end
            
                

                fclose(SNR_textfile_fid);


 
         
        


    
end
        
