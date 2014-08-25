function [image_flag,image_flip_flag,image_weighting_flag,hamming_flag,phase_encod_dir_is,subbas_flag,dwelltime,water_frequency,pat_name,thickness,FOV_X,FOV_Y,freq_enc_dir_shift,freq_enc_dir_step, ...
         csi_file,basis_file,out_dir,image_file,image_flip_file,image_weighting_file,hamming_factor,sddegz,sddegp] = read_Parameter_txt_1_0(Parameter_txt)

     
% read Parameter.txt

%fid_Parameter = fopen('/data/mrsproc/Strasser_Bernhard/PROCESS/Multichannel_Combination/Multichannel_Combination_0.x11/tmp/Parameter.txt', 'r');
%fid_Parameter = fopen('./tmp/Parameter.txt','r');
fid_Parameter = fopen(Parameter_txt,'r');

%read in flags
dummy=fscanf(fid_Parameter,'%s',2);             image_flag=fscanf(fid_Parameter,'%d',1);
dummy=fscanf(fid_Parameter,'%s',1);             image_flip_flag=fscanf(fid_Parameter,'%d',1);
dummy=fscanf(fid_Parameter,'%s',1);             image_weighting_flag=fscanf(fid_Parameter,'%d',1);
dummy=fscanf(fid_Parameter,'%s',1);             hamming_flag=fscanf(fid_Parameter,'%d',1);
dummy=fscanf(fid_Parameter,'%s',3);             phase_encod_dir_is_RL_flag=fscanf(fid_Parameter,'%d',1);
dummy=fscanf(fid_Parameter,'%s',1);             subbas_flag=fscanf(fid_Parameter,'%d',1);

if(phase_encod_dir_is_RL_flag)
    phase_encod_dir_is = 'RL';
else
    phase_encod_dir_is = 'AP';
end



%read CSI Parameter
dummy=fscanf(fid_Parameter,'%s',3);              dwelltime = fscanf(fid_Parameter,'%d',1);            water_frequency = fscanf(fid_Parameter,'%d',1);
dummy=fscanf(fid_Parameter,'%s',1);              pat_name = fscanf(fid_Parameter,'%s',1);
dummy = fscanf(fid_Parameter,'%s',1);            thickness = fscanf(fid_Parameter,'%d',1);  % ZSTEP
dummy = fscanf(fid_Parameter,'%s',1);            FOV_X = fscanf(fid_Parameter,'%d',1);      % BFOV_X
dummy = fscanf(fid_Parameter,'%s',1);            FOV_Y = fscanf(fid_Parameter,'%d',1);      % BFOV_Y



%read IMAGING Parameter
dummy = fscanf(fid_Parameter,'%s',3);            freq_enc_dir_shift = fscanf(fid_Parameter,'%f',1);     
dummy = fscanf(fid_Parameter,'%s',1);            freq_enc_dir_step = fscanf(fid_Parameter,'%f',1);   


%Read FILE PATHS
dummy = fscanf(fid_Parameter,'%s',2);            csi_file = fscanf(fid_Parameter,'%s',1);   
dummy = fscanf(fid_Parameter,'%s',1);            basis_file = fscanf(fid_Parameter,'%s',1);  
dummy = fscanf(fid_Parameter,'%s',1);            out_dir = fscanf(fid_Parameter,'%s',1);
dummy = fscanf(fid_Parameter,'%s',1);            image_file = fscanf(fid_Parameter,'%s',1);
dummy = fscanf(fid_Parameter,'%s',1);            image_flip_file = fscanf(fid_Parameter,'%s',1);
dummy = fscanf(fid_Parameter,'%s',1);            image_weighting_file = fscanf(fid_Parameter,'%s',1);


%Read additional info
dummy = fscanf(fid_Parameter,'%s',4);            hamming_factor = fscanf(fid_Parameter,'%f',1);
dummy = fscanf(fid_Parameter,'%s',1);            sddegz = fscanf(fid_Parameter,'%f',1);
dummy = fscanf(fid_Parameter,'%s',1);            sddegp = fscanf(fid_Parameter,'%f',1);

fclose(fid_Parameter);