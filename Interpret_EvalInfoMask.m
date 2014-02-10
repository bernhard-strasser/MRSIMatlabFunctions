function [Flags,FlagNamesSet] = Interpret_EvalInfoMask(EvalInfoMask_2xuint32)
%
% read_csi_dat Read in csi-data from Siemens raw file format
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [csi,NoiseCorrMat,Noise_mat,kSpace] = read_csi_dat(csi_path, zerofill_to_nextpow2_flag, zerofilling_fact, Hadamard_flag, x_shift,y_shift,NoFFT_flag, NoiseCorrMat)
%
% Input: 
% -         csi_path                    ...     Path of MRS(I) file.
% -         zerofill_to_nextpow2_flag   ...     Flag, if the MRSI data should be zerofilled to the next power of 2 in k-space (e.g. 42x42 sampled --> zf to 64x64?)
% -         zerofilling_fact            ...     Factor with which the MRSI data should be zerofilled in k-space for interpolation (e.g. zerofill from 64x64 to 128x128)
% -         Hadamard_flag               ...     If data is multislice hadamard encoded, perform hadamard-decoding function
% -         x_shift                     ...     Shift the MRSI data in the left-right direction ( = row direction of matrix) by x_shift voxels
% -         y_shift                     ...     Shift the MRSI data in anterior-posterior direction ( = column direction of matrix) by y_shift voxels
% -         NoFFT_flag                  ...     If this is true, don't perform any fft.
% -         NoiseCorrMat                ...     If size(NoiseCorrMat) = [cha cha]: the k-space Data gets decorrelated with this matrix. 
%                                               If NoiseCorrMat = 1: the end of the FIDs of the outermost k-space/image points are used for noise decorrelation.
%                                               If NoiseCorrMat = 0, or not existant: No Noise Decorrelation performed
%
% Output:
% -         csi                         ...     Output data in image domain. In case of Single Voxel Spectroscopy, this is the only output
% -         NoiseCorrMat                ...     The Noise Correlation Matrix in order to check if it was properly computed from the csi data. Is 0 if no decorrelation performed.
% -         Noise_mat                   ...     The Noise gathered from the CSI data. Noise_mat = 0, if not gathered.
% -         kSpace                  ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x vecSize x Averages
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux_1_0,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.






%% 1. Define Flags

FlagsString1 = dec2bin(EvalInfoMask_2xuint32(1),32);
FlagsString2 = dec2bin(EvalInfoMask_2xuint32(2),32);

Flags1 = zeros([1 numel(FlagsString1)]);
Flags2 = zeros([1 numel(FlagsString2)]);

    
for k = 1:numel(FlagsString1)
    Flags1(k) = logical(str2double(FlagsString1(k)));
    Flags2(k) = logical(str2double(FlagsString2(k)));        
end
    
Flags = cat(2,fliplr(Flags1),fliplr(Flags2));
always_zeros = [7,8,10,34,35,36,37,38,39,40,44,45,54,55,56,57,58,59,60,61,62,63,64];

Flags(always_zeros) = [];




%% 2. Associate Flags with Flagnames

numel_FlagNamesSet = 0;
if(Flags(4))
    FlagNamesSet{numel_FlagNamesSet+1} = 'ONLINE';
    numel_FlagNamesSet = numel_FlagNamesSet + 1;
end
if(Flags(7))
    FlagNamesSet{numel_FlagNamesSet+1} = 'LASTSCANINCONCAT';
    numel_FlagNamesSet = numel_FlagNamesSet + 1;

end
if(Flags(9))
    FlagNamesSet{numel_FlagNamesSet+1} = 'LASTSCANINMEAS';
    numel_FlagNamesSet = numel_FlagNamesSet + 1;
end
if(Flags(21))
    FlagNamesSet{numel_FlagNamesSet+1} = 'PATREFANDIMASCAN';
    numel_FlagNamesSet = numel_FlagNamesSet + 1;
end
if(Flags(23))
    FlagNamesSet{numel_FlagNamesSet+1} = 'NOISEADJSCAN';
    numel_FlagNamesSet = numel_FlagNamesSet + 1;
end





