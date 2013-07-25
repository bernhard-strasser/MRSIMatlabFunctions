function [csi,NoiseCorrMat,csi_kspace] = read_csi_1_6(csi_file,zerofill_to_nextpow2_flag,zerofilling_fact, Hadamard_flag, x_shift,y_shift,NoFFT_flag,NoiseCorrMat)
%
% read_csi_x_x Read in csi-data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function only decides if input file is .dat file or DICOM according to its ending. It then calls the read_csi_dat_x_x or read_csi_dicom_x_x
% functions. Refer to these for more info
%
%
% [csi,csi_kspace] = read_csi_1_4(csi_path, zerofill_to_nextpow2_flag, zerofilling_fact, x_shift,y_shift)
%
% Input: 
% -         csi_file                    ...     Path of MRS(I) file.
% -         zerofill_to_nextpow2_flag   ...     Flag, if the MRSI data should be zerofilled to the next power of 2 in k-space (e.g. 42x42 sampled --> zf to 64x64?)
% -         zerofilling_fact            ...     Factor with which the MRSI data should be zerofilled in k-space for interpolation (e.g. zerofill from 64x64 to 128x128)
% -         Hadamard_flag               ...     If data is multislice hadamard encoded, perform hadamard-decoding function
% -         x_shift                     ...     Shift the MRSI data in the left-right direction ( = row direction of matrix) by x_shift voxels
% -         y_shift                     ...     Shift the MRSI data in anterior-posterior direction ( = column direction of matrix) by y_shift voxels
% -         NoFFT_flag                  ...     If true, no fft is performed when reading dat-files.
% -         NoiseCorrMat                ...     If size(NoiseCorrMat) = [cha cha]: the k-space Data gets decorrelated with this matrix. 
%                                               If NoiseCorrMat = 1: the end of the FIDs of the outermost k-space/image points are used for noise decorrelation.
%                                               If NoiseCorrMat = 0, or not existant: No Noise Decorrelation performed
%
% Output:
% -         csi                         ...     Output data in image domain. In case of Single Voxel Spectroscopy, this is the only output
% -         csi_kspace                  ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x vecSize
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: read_csi_dat_1_10, Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m, read_csi_dicom_1_1



%% 0. PREPARATIONS

% Assign standard values to variables if nothing is passed to function.
if(~exist('zerofill_to_nextpow2_flag','var'))
    zerofill_to_nextpow2_flag = 1;
end 
if(~exist('zerofilling_fact','var'))
    zerofilling_fact = 1;
end
if(~exist('Hadamard_flag','var'))
    Hadamard_flag = 0;
end
if(~exist('x_shift','var'))
    x_shift = 0;
end    
if(~exist('y_shift','var'))
    y_shift = 0;
end 
if(~exist('NoFFT_flag','var'))
    NoFFT_flag = false;
end
if(~exist('NoiseCorrMat','var'))
    NoiseCorrMat = false;
end

%% 1. Decide whether input file is DICOM or DAT; Call then appropriate function

if(numel(strfind(csi_file, '.dat')) > 0)
    % Avoid memory problems. If only image space data is needed, pass over only this.    
    if(nargout < 3)
        [csi,NoiseCorrMat] = read_csi_dat_2_4(csi_file, zerofill_to_nextpow2_flag,zerofilling_fact, Hadamard_flag, x_shift,y_shift,false,NoiseCorrMat);
    else
        [csi,NoiseCorrMat,csi_kspace] = read_csi_dat_2_4(csi_file, zerofill_to_nextpow2_flag,zerofilling_fact, Hadamard_flag, x_shift,y_shift,NoFFT_flag,NoiseCorrMat);        
    end
    
else
    
    if(nargout < 3)
        csi = read_csi_dicom_1_2(csi_file,zerofilling_fact, x_shift,y_shift);        
    else
        [csi, csi_kspace] = read_csi_dicom_1_2(csi_file,zerofilling_fact, x_shift,y_shift);
    end
    NoiseCorrMat = 0;
    
end



