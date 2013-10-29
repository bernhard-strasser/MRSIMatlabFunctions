function [image_complex, image_complex_kspace] = read_image_1_1(image_file,fredir_desired, phadir_desired,kspace_corr_flag,fredir_shift,phase_encod_dir,flip,total_channel_no,ROW,COL,SLC)
%% 1. Decide whether input file is DICOM or DAT; Call then appropriate function


if(numel(strfind(image_file, '.dat')) > 0)
    
    % Assign standard values if variables are not inputted.
    if(~exist('fredir_desired','var'))
        fredir_desired = 128;
    end
    if(~exist('phadir_desired','var'))
        phadir_desired = 128;
    end    
    if(~exist('kspace_corr_flag','var'))
        kspace_corr_flag = 0;
    end
    if(~exist('fredir_shift','var'))
        fredir_shift = 0;
    end
    if(~exist('phase_encod_dir','var'))
        phase_encod_dir = 'AP';
    end    
    if(~exist('flip','var'))
        flip = 0;
    end    
    
    [image_complex,image_complex_kspace] = read_image_dat_1_8(image_file,fredir_desired, phadir_desired,kspace_corr_flag,fredir_shift,phase_encod_dir,flip);
    
else
    
    [image_complex, image_complex_kspace] = read_image_dicom_1_0(image_file,fredir_desired, phadir_desired,fredir_shift,total_channel_no,ROW,COL,SLC);
    
end



