function Parameters = read_csi_dicom_ascconv_1_0(csi_dicom_file)
%% 0. Preparations



%% 1. Find & SAVE ASCCONV

sLine = fgets(csi_dicom_file);

while()
    sLine = fgets(csi_dicom_file);
    
    if(~begin_found)
        if(sLine < 0)
            display('Error: No Ascconv found. Aborting.')
            return
        else
            continue
        end
        
    else
        
        
        
    end



