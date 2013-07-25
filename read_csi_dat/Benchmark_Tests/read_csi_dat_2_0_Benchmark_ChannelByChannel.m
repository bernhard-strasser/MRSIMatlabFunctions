function [csi,csi_kspace] = read_csi_dat_2_0_Benchmark_ChannelByChannel(csi_path, zerofill_to_nextpow2_flag, zerofilling_fact, Hadamard_flag, x_shift,y_shift)
%
% read_csi_dat_x_x Read in csi-data from Siemens raw file format
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" 
% and perform some easy corrections like zerofilling, Hadamard decoding, etc.
%
%
% [csi,csi_kspace] = read_csi_dat_1_10(csi_path, zerofill_to_nextpow2_flag, zerofilling_fact, x_shift,y_shift)
%
% Input: 
% -         csi_path                    ...     Path of MRS(I) file.
% -         zerofill_to_nextpow2_flag   ...     Flag, if the MRSI data should be zerofilled to the next power of 2 in k-space (e.g. 42x42 sampled --> zf to 64x64?)
% -         zerofilling_fact            ...     Factor with which the MRSI data should be zerofilled in k-space for interpolation (e.g. zerofill from 64x64 to 128x128)
% -         Hadamard_flag               ...     If data is multislice hadamard encoded, perform hadamard-decoding function
% -         x_shift                     ...     Shift the MRSI data in the left-right direction ( = row direction of matrix) by x_shift voxels
% -         y_shift                     ...     Shift the MRSI data in anterior-posterior direction ( = column direction of matrix) by y_shift voxels
%
% Output:
% -         csi                         ...     Output data in image domain. In case of Single Voxel Spectroscopy, this is the only output
% -         csi_kspace                  ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x vecSize x Averages
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.



%% 0. Preparations

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


% Read info from mdh
ParList_mdh = Analyze_csi_mdh_1_4(csi_path,0);
SLC = ParList_mdh.SLC;
total_k_points = ParList_mdh.total_k_points;
vecSize = ParList_mdh.vecSize;
Averages = ParList_mdh.Averages;
clear ParList_mdh


% Read Info from ascconv Header
ParList = read_ascconv_1_2(csi_path);
ROW_raw = ParList.ROW_raw;
COL_raw = ParList.COL_raw;
if(not(isnan(ParList.ROW)))
    ROW_mdh = ParList.ROW;              % This gives the value FinalMatrixSize in the ascconv header. The k-space values in the mdh are according
else                                    % to that size. E.g. if 42x42 (ROW_raw x COL_raw) values were measured and 
    ROW_mdh = ROW_raw;                  % ROW_mdh x COL_mdh = 64x64, then the k-space center point in the mdh is not 42/2+1 = 22, but 64/2+1 = 33;
end                                     % In SVS these entries in the ascconv header do not exist.
if(not(isnan(ParList.COL)))
    COL_mdh = ParList.COL;
else
    COL_mdh = COL_raw;
end
total_channel_no = ParList.total_channel_no;
clear ParList


% zero filling to nextpow2 (e.g. from 42x42 to 64x64)
if(zerofill_to_nextpow2_flag)
    ROW = 2^nextpow2(ROW_raw);
    COL = 2^nextpow2(COL_raw);
else
    ROW = ROW_raw;
    COL = COL_raw;
end


% In case of additional zerofilling (e.g. from 64x64 to 128x128)
ROW = ROW * zerofilling_fact;
COL = COL * zerofilling_fact;




% Compute the shift between the mdh k-space-center (that numbers that are written in the mdh) and the k-space center of the desired matrix size
% This value is used to map the measured k-space to the wanted one.
kspaceCenter_ROW_desired = floor(ROW/2) + 1;
kspaceCenter_COL_desired = floor(COL/2) + 1;

kspaceCenter_ROW_mdh = floor(ROW_mdh/2) + 1;    % Compute the k-space center that is written in the mdh.
kspaceCenter_COL_mdh = floor(COL_mdh/2) + 1;

kspaceCenter_ROW_shift = kspaceCenter_ROW_mdh - kspaceCenter_ROW_desired;
kspaceCenter_COL_shift = kspaceCenter_COL_mdh - kspaceCenter_COL_desired;






%% 1. READ DATA

tic

csi_fid = fopen(sprintf('%s', csi_path),'r');
headersize = fread(csi_fid,1, 'uint32');
fseek(csi_fid, headersize,'bof');                                                                     
csi_kspace = zeros(total_channel_no,ROW,COL,SLC,vecSize,Averages); 

for ADC_MeasNo = 1:Averages*total_k_points*total_channel_no

    chak_header = fread(csi_fid, 64, 'int16');
    k_x = chak_header(17) + 1 - kspaceCenter_ROW_shift;                     
    k_y = chak_header(22) + 1 - kspaceCenter_COL_shift;        
    if(ge(chak_header(19),1))
        k_z = chak_header(19) +1;                                    % SAYS WHICH REPETITION FOR HADAMARD ENCODING OF THE SAME K-POINT IS MEASURED
    else
        k_z = 1;                                                     % some distinction between multislice/hadamard and real 3d-CSI necessary!
    end
    Avg = chak_header(24) + 1;                                       % Averages
    channel_no = chak_header(63) + 1;
    chak_data = fread(csi_fid, vecSize*2, 'float32');                % Read real & imaginary (--> vecSize*2) measured points
    csi_real = chak_data(1:2:end);
    csi_imag = chak_data(2:2:end);
    csi_kspace(channel_no,k_x,k_y,k_z,:,Avg) = complex(csi_real,csi_imag);        
            
end

fclose(csi_fid);

display([ 'Read in took ... ' num2str(toc) ' seconds' ])





%% 2. FFT FROM K_SPACE TO DIRECT SPACE


if(size(csi_kspace,2) > 1 || size(csi_kspace,3) > 1)                    % In case that MRSI data, not Single Voxel Spectroscopy data is inputted.
    
    
    toc_sum = 0;
    
    csi = csi_kspace;
    if(nargout == 1)
        clear csi_kspace;
    end    
    
    for channel_no = 1: size(csi,1)
        
        channel_no
        
        
        tic
        csi_channel = ifftshift(ifftshift(csi(channel_no,:,:,:,:,:),2),3);
        display([ 'ifftshift took     ... ' num2str(toc) ' seconds' ])
        toc_sum = toc_sum + toc;
        
        tic
        csi_channel = fft(fft(csi_channel,[],2),[],3);
        display([ 'fft took     ... ' num2str(toc) ' seconds' ])       
        toc_sum = toc_sum + toc;        
        
        tic
        csi_channel = conj(csi_channel);                                                    % the chem shift gets higher from right to left --> conj reverses that
        display([ 'conj took     ... ' num2str(toc) ' seconds' ])               
        toc_sum = toc_sum + toc;
         
        tic
        csi_channel = fftshift(fftshift(csi_channel,2),3);
        display([ 'fftshift took     ... ' num2str(toc) ' seconds' ])                       
        toc_sum = toc_sum + toc;
        
        tic
        csi(channel_no,:,:,:,:,:) = flipdim(csi_channel,2);  
        if(channel_no == size(csi,1))
            clear csi_channel
        end
        display([ 'flipdim took     ... ' num2str(toc) ' seconds' ])                              
        toc_sum = toc_sum + toc;
        
    end

    display([ 'FFT PROCESS took     ... ' num2str(toc_sum) ' seconds' ])
    
else
    
    csi = conj(csi_kspace);
    csi_kspace = 0;
    
end





%% 3. Spatial Shift

if(ne(x_shift,0))
    csi = circshift(csi, [0 x_shift 0 0 0]);
end

if(ne(y_shift,0))
    csi = circshift(csi, [0 0 y_shift 0 0]);
end





%% 4. Hadamard Decoding

if(Hadamard_flag == 1 && SLC>1)
    csi = hadamard_encoding_7(csi);
end



