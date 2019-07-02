function Read_n_Write_SNR(file_struct,SNR_textfile)

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
    mask = read_RawFiles_1_0(file_struct(read_index).Mask_VC,'float32',ROW,COL,1);
    mask_AC = read_RawFiles_1_0(file_struct(read_index).Mask_AC,'float32',ROW,COL,1);    
    
    
    % 1.2 LCM SNR maps
    
    LCM_VC_SNR = read_RawFiles_1_0(file_struct(read_index).LCM_VC,'float32',ROW,COL,1);
    LCM_AC_SNR = read_RawFiles_1_0(file_struct(read_index).LCM_AC,'float32',ROW,COL,1);  
    
    
    % 1.3 MATLAB SNR maps
    
    MATLAB_VC_SNR = read_RawFiles_1_0(file_struct(read_index).MATLAB_VC,'float32',ROW,COL,1);
    MATLAB_AC_SNR = read_RawFiles_1_0(file_struct(read_index).MATLAB_AC,'float32',ROW,COL,1);    
            

    if(strcmp(file_struct(read_index).Mask_VC,file_struct(read_index).Mask_AC))
        LCM_Ratio_SNR = LCM_AC_SNR ./ LCM_VC_SNR;
        MATLAB_Ratio_SNR = MATLAB_AC_SNR ./ MATLAB_VC_SNR;
    else
        LCM_Ratio_SNR = NaN * ones([ROW,COL]);
        MATLAB_Ratio_SNR = NaN * ones([ROW,COL]);
        LCM_VC_SNR(~mask) = NaN;
        MATLAB_VC_SNR(~mask) = NaN;
        LCM_AC_SNR(~mask_AC) = NaN;
        MATLAB_AC_SNR(~mask_AC) = NaN;
    end

    
    

    %% 2. PLOT MAPS


    % 6.1.1 masks

    fig_mask = figure('Visible','Off');
    imagesc(mask,[0,1]);
    title(sprintf('mask %s', file_struct(read_index).name), 'Interpreter', 'none')

    if(~strcmp(file_struct(read_index).Mask_VC,file_struct(read_index).Mask_AC))
        fig_mask_AC = figure('Visible','Off');
        imagesc(mask_AC,[0,1]);
        title(sprintf('mask AC %s', file_struct(read_index).name), 'Interpreter', 'none')
    end


    if(eq(ROW,32))
        lim_LCM_ACVC = 85;
        lim_MATLAB_ACVC = 1500;   
        lim_LCM_Ratio = 9;
        lim_MATLAB_Ratio = 9;
    else
        lim_LCM_ACVC = 25;
        lim_MATLAB_ACVC = 110;  
        lim_LCM_Ratio = 7;
        lim_MATLAB_Ratio = 7;            
    end




    fig_LCM_VC = figure('Visible','Off');
    imagesc(LCM_VC_SNR, [0 lim_LCM_ACVC]);
    title(sprintf('LCM VC SNR map %s', file_struct(read_index).name), 'Interpreter', 'none')  
    colorbar('FontSize',30)

    fig_LCM_AC = figure('Visible','Off');
    imagesc(LCM_AC_SNR, [0 lim_LCM_ACVC]);
    title(sprintf('LCM AC SNR map %s', file_struct(read_index).name), 'Interpreter', 'none')  
    colorbar('FontSize',30)

    fig_MATLAB_VC = figure('Visible','Off');
    imagesc(MATLAB_VC_SNR,[0 lim_MATLAB_ACVC]);
    title(sprintf('MATLAB VC SNR map %s', file_struct(read_index).name), 'Interpreter', 'none')  
    colorbar('FontSize',30)

    fig_MATLAB_AC = figure('Visible','Off');
    imagesc(MATLAB_AC_SNR, [0 lim_MATLAB_ACVC]);
    title(sprintf('MATLAB AC SNR map %s', file_struct(read_index).name), 'Interpreter', 'none')  
    colorbar('FontSize',30)

    fig_LCM_Ratio = figure('Visible','Off');
    imagesc(LCM_Ratio_SNR, [0 lim_LCM_Ratio]);
    title(sprintf('MATLAB Ratio AC/VC SNR map %s', file_struct(read_index).name), 'Interpreter', 'none')  
    colorbar('FontSize',30)


    fig_MATLAB_Ratio = figure('Visible','Off');
    imagesc(MATLAB_Ratio_SNR, [0 lim_MATLAB_Ratio]);
    title(sprintf('MATLAB Ratio AC/VC SNR map %s', file_struct(read_index).name), 'Interpreter', 'none')  
    colorbar('FontSize',30)  


    saveas(fig_mask,sprintf('%s/mask',file_struct(read_index).OutDir),'epsc2')
    saveas(fig_mask,sprintf('%s/mask',file_struct(read_index).OutDir),'fig')
    saveas(fig_mask,sprintf('%s/mask',file_struct(read_index).OutDir),'jpg')    
    if(~strcmp(file_struct(read_index).Mask_VC,file_struct(read_index).Mask_AC))
        saveas(fig_mask_AC,sprintf('%s/mask_AC',file_struct(read_index).OutDir),'epsc2')
        saveas(fig_mask_AC,sprintf('%s/mask_AC',file_struct(read_index).OutDir),'fig')
        saveas(fig_mask_AC,sprintf('%s/mask_AC',file_struct(read_index).OutDir),'jpg')        
    end    



    saveas(fig_LCM_VC,sprintf('%s/SNR_LCM_VC',file_struct(read_index).OutDir),'epsc2')
    saveas(fig_LCM_VC,sprintf('%s/SNR_LCM_VC',file_struct(read_index).OutDir),'fig')
    saveas(fig_LCM_VC,sprintf('%s/SNR_LCM_VC',file_struct(read_index).OutDir),'jpg')   

    saveas(fig_LCM_AC,sprintf('%s/SNR_LCM_AC',file_struct(read_index).OutDir),'epsc2')
    saveas(fig_LCM_AC,sprintf('%s/SNR_LCM_AC',file_struct(read_index).OutDir),'fig')
    saveas(fig_LCM_AC,sprintf('%s/SNR_LCM_AC',file_struct(read_index).OutDir),'jpg')    

    saveas(fig_LCM_Ratio,sprintf('%s/SNR_LCM_Ratio',file_struct(read_index).OutDir),'epsc2')
    saveas(fig_LCM_Ratio,sprintf('%s/SNR_LCM_Ratio',file_struct(read_index).OutDir),'fig')        
    saveas(fig_LCM_Ratio,sprintf('%s/SNR_LCM_Ratio',file_struct(read_index).OutDir),'jpg')        


    saveas(fig_MATLAB_VC,sprintf('%s/SNR_MATLAB_VC',file_struct(read_index).OutDir),'epsc2')
    saveas(fig_MATLAB_VC,sprintf('%s/SNR_MATLAB_VC',file_struct(read_index).OutDir),'fig')
    saveas(fig_MATLAB_VC,sprintf('%s/SNR_MATLAB_VC',file_struct(read_index).OutDir),'jpg')
    
    
    saveas(fig_MATLAB_AC,sprintf('%s/SNR_MATLAB_AC',file_struct(read_index).OutDir),'epsc2')
    saveas(fig_MATLAB_AC,sprintf('%s/SNR_MATLAB_AC',file_struct(read_index).OutDir),'fig')
    saveas(fig_MATLAB_AC,sprintf('%s/SNR_MATLAB_AC',file_struct(read_index).OutDir),'jpg')
    
    
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
                
                fprintf(SNR_textfile_fid,'\nVoxel (ROW,COL,SLC)\t\tLCM_VC\t\tLCM_AC\t\tLCM_Ratio\t\tMAT_VC\t\tMAT_AC\t\tMAT_Ratio\n');                  
                
                
                for x = 1:ROW
                    for y = 1:COL
                        for z = 1:SLC
                            if(mask(x,y,z) == 0 && mask_AC(x,y,z) == 0)
                                continue;
                            end
                            
                            fprintf(SNR_textfile_fid,'(%d,%d,%d)\t\t\t\t%7.2f\t\t%7.2f\t\t%7.3f\t\t%7.2f\t\t%7.2f\t\t%7.3f\n', ...
                            x,y,z,LCM_VC_SNR(x,y),LCM_AC_SNR(x,y),LCM_Ratio_SNR(x,y),MATLAB_VC_SNR(x,y),MATLAB_AC_SNR(x,y),MATLAB_Ratio_SNR(x,y));  
                        
                        end
                    end
                end
            
                

                fclose(SNR_textfile_fid);


 
         
        


    
end
        
