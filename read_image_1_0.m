function [image_complex, image_complex_kspace] = read_image_1_0(image_file,fredir_desired, phadir_desired,kspace_corr_flag,fredir_shift,phase_encod_dir,flip,total_channel_no,ROW,COL,SLC)
%% 1. Decide whether input file is DICOM or DAT; Call then appropriate function


if(numel(strfind(image_file, '.dat')) > 0)
    
    [image_complex,image_complex_kspace] = read_image_dat_1_7(image_file,fredir_desired, phadir_desired,kspace_corr_flag,fredir_shift,phase_encod_dir,flip);
    
else
    
    [image_complex, image_complex_kspace] = read_image_dicom_1_0(image_file,fredir_desired, phadir_desired,fredir_shift,total_channel_no,ROW,COL,SLC);
    
end



