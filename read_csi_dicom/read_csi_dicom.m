function [csi, csi_kspace] = read_csi_dicom_1_2(csi_path, zerofilling_fact, x_shift,y_shift)
%
% read_csi_dicom_x_x Read in csi-data from DICOM file format.
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the DICOM file format.
%
%
% [csi,csi_kspace] = read_csi_dicom_1_1(csi_path, zerofilling_fact, x_shift,y_shift)
%
% Input: 
% -         csi_path                    ...     Path of spectroscopy (imaging) file.
% -         zerofilling_fact            ...     Factor with which the MRSI data should be zerofilled in k-space for interpolation (e.g. zerofill from 64x64 to 128x128)
% -         x_shift                     ...     Shift the MRSI data in the left-right direction ( = row direction of matrix) by x_shift voxels
% -         y_shift                     ...     Shift the MRSI data in anterior-posterior direction ( = column direction of matrix) by y_shift voxels
%
% Output:
% -         csi                         ...     Output data in image domain. In case of Single Voxel Spectroscopy, this is the only output
% -         csi_kspace                  ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x vecSize
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: read_ascconv_1_2

% Further remarks: This function uses IFFTs in ROW-dimension to get from image to k-space. This is mathematically wrong, but Siemens seems
% to do the same when creating DICOMS. The only difference is that the k-space is flipped up/down by that. In COL dimension it uses FFT. This was
% chosen in order to get the same result in k-space when the k-space of Siemens raw data (.dat-files) is read in.


%% 0. Preparations

% Assign standard values to variables if nothing is passed to function.
if(~exist('zerofilling_fact','var'))
    zerofilling_fact = 1;
end

if(~exist('x_shift','var'))
    x_shift = 0;
end

if(~exist('y_shift','var'))
    y_shift = 0;
end


% Read info from ascconv header
ParList = read_ascconv_1_8(csi_path);
total_channel_no = ParList.total_channel_no_reco;

if(not(isnan(ParList.nFreqEnc_FinalMatrix)))
    ROW = ParList.nFreqEnc_FinalMatrix;              % This gives the value FinalMatrixSize in the ascconv header, so the Matrix Size that is written to DICOM
else                                    % file. In SVS these entries in the ascconv header do not exist.
    ROW = ParList.nFreqEnc;                  
end                                     
if(not(isnan(ParList.nPhasEnc_FinalMatrix)))
    COL = ParList.nPhasEnc_FinalMatrix;
else
    COL = ParList.nPhasEnc;
end
if(not(isnan(ParList.nPhasEnc_FinalMatrix)))
    SLC = ParList.nSLC_FinalMatrix;
else
    SLC = ParList.nPartEnc;
end


SLC = ParList.nSLC * SLC;
vecSize = ParList.vecSize;
clear ParList



%% 1. READ & Reshape data

tic

% Read Data
csi_fid = fopen(csi_path,'r');
fseek(csi_fid, -((vecSize*ROW*COL*SLC*2*4)), 'eof');  % Find position where to start to read. Go from end of file back vecSize*ROW*COL*SLC * 2 (real & imaginary) * 4 (float32) bytes.
csi_data = fread(csi_fid,'float32');
fclose(csi_fid);


% Separate real & imaginary part
csi_real=csi_data(1:2:end);                           % Odd points are real data
csi_imag=csi_data(2:2:end);                           % Even imaginary
csi_complex=complex(csi_real,csi_imag);


% Reshape to 5D-Matrix
csi = zeros(total_channel_no,ROW,COL,SLC,vecSize);

k=-1;
for z=1:SLC
    for y=1:ROW
        for x=1:COL
            k=k+1;
            csi(1,x,y,z,:) = csi_complex(k*vecSize+1:(k+1)*vecSize);	% total channel number not yet implemented, since ICE always creates 1 fiile per channel
        end                                                             % there is anyway only one channel in the DICOM file.
    end
end


display(['Read in took ... ' num2str(toc) 'seconds.'])




%% 2. Spatial Shift

if(ne(x_shift,0))
    csi = circshift(csi, [0 x_shift 0 0 0]);
end

if(ne(y_shift,0))
    csi = circshift(csi, [0 0 y_shift 0 0]);
end



%% 3. Zero filling in kspace

if(ne(zerofilling_fact,1))
    
    tic
    
    % Go to k-space
    csi_kspace = ifftshift(ifftshift(csi,2),3);
    csi_kspace = ifft(ifft(csi_kspace,[],2),[],3);
    csi_kspace = fftshift(fftshift(csi_kspace,2),3);
    
    % Initialize zerofilled matrix
    csi_kspace_zf = zeros([total_channel_no ROW*zerofilling_fact COL*zerofilling_fact SLC vecSize]);
    
    % Write old k-space data in center of new k-space
    csi_kspace_zf(:,ROW*zerofilling_fact/4+1:ROW*zerofilling_fact*3/4,COL*zerofilling_fact/4+1:COL*zerofilling_fact*3/4,:,:) = csi_kspace;
    
    % overwrite variable name, clear other variable name
    %clear csi_kspace
    csi_kspace = csi_kspace_zf;
    clear csi_kspace_zf;
    
    
    % Go back to image space
    csi = ifftshift(ifftshift(csi_kspace,2),3);
    csi = ifft(ifft(csi,[],2),[],3);
    csi = fftshift(fftshift(csi,2),3);    
    
    
    display(['Zerofilling took ... ' num2str(toc) 'seconds'])
    
end



%% 4. Fourier transform to k-space if not done yet & if necessary

if(nargout > 1 && not(exist('csi_kspace','var')))
    
    tic
    
	csi_kspace = ifftshift(ifftshift(csi,2),3);
    csi_kspace = ifft(ifft(csi_kspace,[],2),[],3);
    csi_kspace = fftshift(fftshift(csi_kspace,2),3);
    csi_kspace = conj(csi_kspace);
    
    display([ 'IFFT took    ... ' num2str(toc) ' seconds' ])

end


