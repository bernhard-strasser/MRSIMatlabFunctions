function Parameters = read_csi_dicom_ascconv_0_2(csi_dicom_file)
%% 0. Preparations



%% 1. Find & SAVE ASCCONV

begin_found = 0;
%end_found = 0;
ascconv = [];
sLine = 0;
%sLine = fgets(csi_dicom_file);

while(sLine > -1)
    
    sLine = fgets(csi_dicom_file);
    
    if(~begin_found)                                % If begin of ascconv not yet found
        
        
        if(strcmp('### ASCCONV BEGIN ###',sLine))
            begin_found = true;
        else
            continue                                % If begin not yet found and current line is also not begin --> read next line
        end
        
        
    else                                            % If begin of ascconv was found
        
        if(strcmp('### ASCCONV END ###',sLine))     % If the end was found --> stop while loop
            return
        else
            ascconv = [ascconv; {sLine}];           % If the begin was already found and the current line is not the end --> read in line and save it
        end
        
    end
  

        
        
        
end



