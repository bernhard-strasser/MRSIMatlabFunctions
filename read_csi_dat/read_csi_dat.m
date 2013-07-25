function [csi,NoiseCorrMat,Noise_mat,csi_kspace] = read_csi_dat(csi_path, zerofill_to_nextpow2_flag, zerofilling_fact, Hadamard_flag, x_shift,y_shift,NoFFT_flag, NoiseCorrMat)
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
% -         NoFFT_flag                  ...     If this is true, don't perform any fft.
% -         NoiseCorrMat                ...     If size(NoiseCorrMat) = [cha cha]: the k-space Data gets decorrelated with this matrix. 
%                                               If NoiseCorrMat = 1: the end of the FIDs of the outermost k-space/image points are used for noise decorrelation.
%                                               If NoiseCorrMat = 0, or not existant: No Noise Decorrelation performed
%
% Output:
% -         csi                         ...     Output data in image domain. In case of Single Voxel Spectroscopy, this is the only output
% -         NoiseCorrMat                ...     The Noise Correlation Matrix in order to check if it was properly computed from the csi data. Is 0 if no decorrelation performed.
% -         Noise_mat                   ...     The Noise gathered from the CSI data. Noise_mat = 0, if not gathered.
% -         csi_kspace                  ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x vecSize x Averages
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux_1_0,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.





%% 0. Preparations


% Find out memory used by MATLAB
memused_before = memused_linux_1_0(1); 



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
if(~exist('Noise_mat','var'))
    Noise_mat = 0;
end


% Read info from mdh
ParList_mdh = Analyze_csi_mdh_1_4(csi_path,0);
SLC = ParList_mdh.SLC;
total_k_points = ParList_mdh.total_k_points;
vecSize = ParList_mdh.vecSize;
Averages = ParList_mdh.Averages;
clear ParList_mdh


% Read Info from ascconv Header
ParList = read_ascconv_1_8(csi_path);
ROW_raw = ParList.nPhasEnc;
COL_raw = ParList.nFreqEnc;
if(not(isnan(ParList.nPhasEnc_FinalMatrix)))
    ROW_mdh = ParList.nPhasEnc_FinalMatrix;              % This gives the value FinalMatrixSize in the ascconv header. The k-space values in the mdh are according
else                                                     % to that size. E.g. if 42x42 (ROW_raw x COL_raw) values were measured and 
    ROW_mdh = ROW_raw;                                   % ROW_mdh x COL_mdh = 64x64, then the k-space center point in the mdh is not 42/2+1 = 22, but 64/2+1 = 33;
end                                                      % In SVS these entries in the ascconv header do not exist.
if(not(isnan(ParList.nFreqEnc_FinalMatrix)))
    COL_mdh = ParList.nFreqEnc_FinalMatrix;
else
    COL_mdh = COL_raw;
end
total_channel_no = ParList.total_channel_no_reco;
ThreeD_flag = ParList.ThreeD_flag;


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





%% 1. Define k-space centers

% Compute the shift between the mdh k-space-center (that numbers that are written in the mdh) and the k-space center of the desired matrix size
% This value is used to map the measured k-space to the wanted one.
kspaceCenter_ROW_desired = floor(ROW/2) + 1;
kspaceCenter_COL_desired = floor(COL/2) + 1;

kspaceCenter_ROW_mdh = floor(ROW_mdh/2) + 1;    % Compute the k-space center that is written in the mdh.
kspaceCenter_COL_mdh = floor(COL_mdh/2) + 1;

kspaceCenter_ROW_shift = kspaceCenter_ROW_mdh - kspaceCenter_ROW_desired;
kspaceCenter_COL_shift = kspaceCenter_COL_mdh - kspaceCenter_COL_desired;





%% 2. READ DATA

tic
fprintf('\nReading data\t\t\t...')
        
csi_fid = fopen(sprintf('%s', csi_path),'r');
headersize = fread(csi_fid,1, 'uint32');
fseek(csi_fid, headersize,'bof');                                                                     
csi_kspace = zeros([total_channel_no,ROW,COL,SLC,vecSize,Averages]);

for ADC_MeasNo = 1:Averages*total_k_points
    for channel_no = 1:total_channel_no
        
        chak_header = fread(csi_fid, 64, 'int16');
        k_x = chak_header(17) + 1 - kspaceCenter_ROW_shift;                     
        k_y = chak_header(22) + 1 - kspaceCenter_COL_shift;        
        if(ge(chak_header(19),1))
            k_z = chak_header(19) +1;                                    % SAYS WHICH REPETITION FOR HADAMARD ENCODING OF THE SAME K-POINT IS MEASURED
        else
            k_z = 1;                                                     % some distinction between multislice/hadamard and real 3d-CSI necessary!
        end
        Avg = chak_header(24) + 1;                                       % Averages
        %channel_no = chak_header(63) + 1;
        chak_data = fread(csi_fid, vecSize*2, 'float32');                % Read real & imaginary (--> vecSize*2) measured points
        csi_real = chak_data(1:2:end);
        csi_imag = chak_data(2:2:end);
        csi_kspace(channel_no,k_x,k_y,k_z,:,Avg) = complex(csi_real,csi_imag);   
        
    end       
end

fclose(csi_fid);

fprintf('\ttook\t%10.6f seconds',toc)       



%% 3. Perform Noise Decorrelation

if(exist('NoiseCorrMat','var'))
    if(numel(NoiseCorrMat) > 1)
    
        fprintf('\nNoise Decorrelating Using\t...\tHanded Over NoiseCorr Matrix.')    
        %DecorrelationMatrix = inv(chol(Noise_CorrMat/2,'lower'));
        csi_kspace = reshape(chol(NoiseCorrMat/2,'lower') \ reshape(csi_kspace, [size(csi_kspace,1) numel(csi_kspace)/size(csi_kspace,1)]), size(csi_kspace));    % Matrix multiplication
        

    elseif(NoiseCorrMat)
        fprintf('\nNoise Decorrelating Using\t...\tCSI Data Itself.')    
        % Gather "Noise" of the CSI data (end of FIDs in the outermost acquired circle of k-space points). Only working for 2D and Multislice
        nFIDendpoints = 400;
        TakeNPointsOutOfEnd = 200;
        randy = randperm(nFIDendpoints); randy = randy(1:TakeNPointsOutOfEnd); % Take 'em randomly
        % Take random points at end of FID
        Noise_csi = csi_kspace(:,:,:,:,end - (nFIDendpoints - 1) : end); Noise_csi = Noise_csi(:,:,:,:,randy);

        % Take only csi-points which are farest away from k-space center (so a circle with radius 31 k-space points)
        if(ParList.Full_ElliptWeighted_Or_Weighted_Acq == 2)
            [Elliptical_dummy,OuterkSpace_mask1] = EllipticalFilter_1_1(squeeze(csi_kspace(1,:,:,1,1)),[1 2],[1 1 1 size(csi_kspace,3)/2-1],1);
            [Elliptical_dummy,OuterkSpace_mask2] = EllipticalFilter_1_1(squeeze(csi_kspace(1,:,:,1,1)),[1 2],[1 1 1 size(csi_kspace,3)/2-2],1);
            OuterkSpace_mask = OuterkSpace_mask1 - OuterkSpace_mask2;
        else
            OuterkSpace_mask = zeros([size(csi_kspace,2), size(csi_kspace,3)]);
            OuterkSpace_mask(1:end,1) = 1; OuterkSpace_mask(1,1:end) = 1; OuterkSpace_mask(end,1:end) = 1; OuterkSpace_mask(1:end,end) = 1;
        end
        PI_mask = abs(squeeze(csi_kspace(1,:,:,1,1))); PI_mask(PI_mask > 0) = 1;
        csi_mask = OuterkSpace_mask .* PI_mask;
        csi_mask = repmat(logical(csi_mask), [1 1 size(csi_kspace,1)*size(csi_kspace,4)*TakeNPointsOutOfEnd]);
        csi_mask = reshape(csi_mask, [size(csi_kspace,2) size(csi_kspace,3) size(csi_kspace,1) size(csi_kspace,4) TakeNPointsOutOfEnd]);
        csi_mask = permute(csi_mask, [3 1 2 4 5]);

        Noise_mat = Noise_csi(csi_mask);
        Noise_mat = reshape(Noise_mat, [total_channel_no numel(Noise_mat)/total_channel_no]);  

        % Compute NoiseDecorrelation of CSI data
        NoiseCorrMat = 1/(size(Noise_mat,2)) * (Noise_mat * Noise_mat');    

        % Apply Noise Decorrelation
        csi_kspace = reshape(chol(NoiseCorrMat/2,'lower') \ reshape(csi_kspace, [size(csi_kspace,1) numel(csi_kspace)/size(csi_kspace,1)]), size(csi_kspace));    % Matrix multiplication   

        clear nFIDendpoints TakeNPointsOutOfEnd randy Noise_csi Elliptical_dummy OuterkSpace_mask2 OuterkSpace_mask1 OuterkSpace_mask PI_mask csi_mask Noise_mat NoiseCorrMat
        
        
    end
    
end






%% 4. Hadamard Decoding

if(Hadamard_flag == 1 && SLC>1 && ~ThreeD_flag)
    for channel_no = 1:total_channel_no
        
        fprintf('\nHadamard decoding channel %02d\t...', channel_no)
        tic

        csi_kspace(channel_no,:,:,:,:,:) = hadamard_decoding_9(csi_kspace(channel_no,:,:,:,:,:));            
     
        fprintf('\ttook\t%10.6f seconds', toc)              
        
    end
end





%% 5. FFT FROM K_SPACE TO DIRECT SPACE

if(~NoFFT_flag)
    if(size(csi_kspace,2) > 1 || size(csi_kspace,3) > 1)                    % In case that MRSI data, not Single Voxel Spectroscopy data is inputted.

        csi = csi_kspace;


        if(nargout < 3)
            clear csi_kspace;
        end    

        
        % Sacher: Doing everythin at once, because enough memory is available.
        [bla, hostname] = unix('echo $(hostname)'); clear bla;
        if(~isempty(strfind(hostname,'sacher')))
 
                tic_overall = tic;
                csi = ifftshift(ifftshift(csi,2),3);
                csi = fft(fft(csi,[],2),[],3);
                if(~Hadamard_flag && SLC > 1)                       % This is a hack for now to perform fft also in z-direction, if data has several slices but not hadamard encoded.
                    csi = ifftshift(csi,4);                 % In the long run, this must be changed. It must be read out from the header if it is real 3D,2D or 2D+Hadamard.
                    csi = fft(csi,[],4); 
                end                
                csi = conj(csi);                            % the chem shift gets higher from right to left --> conj reverses that
                csi = fftshift(fftshift(csi,2),3); 
                if(~Hadamard_flag && SLC > 1)                       % Same hack as above
                    csi = fftshift(csi,4);
                end                
                csi = flipdim(csi,2);                       % THIS FLIPS LEFT AND RIGHT IN SPATIAL DOMAIN BECAUSE PHYSICIANS WANT TO SEE IMAGES FLIPPED 
                fprintf('\nOverall FFT Process\t\t...\ttook\t%10.6f seconds', toc(tic_overall))

        % All others: Everything channel by channel.
        else
          
            tic_overall = tic;
            for channel_no = 1: size(csi,1)

                tic_loop = tic;
                fprintf('\nFouriertransforming channel %02d\t...', channel_no)

          


                csi_channel = csi(channel_no,:,:,:,:,:);                    % This extra assignment proved to be faster than using always csi(channel_no,:,:,:,:,:). Seems that indexing is rather slow.

                csi_channel = ifftshift(ifftshift(csi_channel,2),3);
                csi_channel = fft(fft(csi_channel,[],2),[],3);

                if(~Hadamard_flag && SLC > 1)                       % This is a hack for now to perform fft also in z-direction, if data has several slices but not hadamard encoded.
                    csi_channel = ifftshift(csi_channel,4);                 % In the long run, this must be changed. It must be read out from the header if it is real 3D,2D or 2D+Hadamard.
                    csi_channel = fft(csi_channel,[],4); 
                end

                csi_channel = conj(csi_channel);                            % the chem shift gets higher from right to left --> conj reverses that
                csi_channel = fftshift(fftshift(csi_channel,2),3);
                if(~Hadamard_flag && SLC > 1)                       % Same hack as above
                    csi_channel = fftshift(csi_channel,4);
                end
                csi_channel = flipdim(csi_channel,2);                       % THIS FLIPS LEFT AND RIGHT IN SPATIAL DOMAIN BECAUSE PHYSICIANS WANT TO SEE IMAGES FLIPPED

                csi(channel_no,:,:,:,:,:) = csi_channel; 

                fprintf('\ttook\t%10.6f seconds', toc(tic_loop))       

            end

            clear csi_channel
            fprintf('\nOverall FFT Process\t\t...\ttook\t%10.6f seconds', toc(tic_overall))
        
        end

    else

        csi = conj(csi_kspace);
        csi_kspace = 0;

    end
    
else
    
    csi = 0;
    
end





%% 6. Spatial Shift

if(ne(x_shift,0))
    csi = circshift(csi, [0 x_shift 0 0 0]);
end

if(ne(y_shift,0))
    csi = circshift(csi, [0 0 y_shift 0 0]);
end







%% 7. Postparations

memused_after = memused_linux_1_0(1); 
display([char(10) 'The function used ' num2str(memused_after-memused_before) '% of the total memory.'])


