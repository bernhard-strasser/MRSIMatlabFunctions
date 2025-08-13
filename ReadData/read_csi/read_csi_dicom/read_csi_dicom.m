function MRStruct = read_csi_dicom(csi_path, Settings)
%
% read_csi_dicom Read in csi-data from DICOM file format.
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the DICOM file format.
%
%
% [MRStruct] = read_csi_dicom(csi_path, Settings)
%
% Input: 
% -         csi_path                    ...     Path of spectroscopy (imaging)DICOM  file or folder with DICOM files.
% -         Settings                    ...     Variable with fields controlling the read-in:
%       - .zerofilling_fact             ...     Factor with which the MRSI data should be zerofilled in k-space for interpolation (e.g. zerofill from 64x64 to 128x128)
%       - .x_shift                      ...     Shift the MRSI data in the left-right direction ( = row direction of matrix) by x_shift voxels
%       - .y_shift                      ...     Shift the MRSI data in anterior-posterior direction ( = column direction of matrix) by y_shift voxels
%       - .FTTokSpace_flag              ...     If true the data will be Fourier transformed to k-space, and the MRStruct.Data will be in k-Space
%       - .RemoveFIDZerofilling_flag    ...     Remove zeros at the end of FID which is done by ICE.
%
% Output:
% -         MRStruct                    ...     Output with subfields:
%                - .Data                ...     Data in image domain if Settings.FTTokSpace_flag is false, or in k-Space if it is true.
%                - .Par                 ...    Parameters of Data
%                - .RecoPar             ...    Parameters of Data after all reconstruction steps were done
%                - .RecoSteps           ...    Saves all the settings of all reconstruction steps
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: read_ascconv

% Further remarks: This function uses IFFTs in ROW-dimension to get from image to k-space. This is mathematically wrong, but Siemens seems
% to do the same when creating DICOMS. The only difference is that the k-space is flipped up/down by that. In COL dimension it uses FFT. This was
% chosen in order to get the same result in k-space when the k-space of Siemens raw data (.dat-files) is read in.


%% 0. Preparations

% Assign standard values to variables if nothing is passed to function.
if(~exist('Settings','var'))
    Settings = struct; 
end
if(~isfield(Settings,'zerofilling_fact'))
    Settings.zerofilling_fact = 1;
end

if(~isfield(Settings,'x_shift'))
    Settings.x_shift = 0;
end

if(~isfield(Settings,'y_shift'))
    Settings.y_shift = 0;
end

if(~isfield(Settings,'FTTokSpace_flag'))
    Settings.FTTokSpace_flag = 0;
end
if(~isfield(Settings,'RemoveFIDZerofilling_flag'))
    Settings.RemoveFIDZerofilling_flag = 1;
end

% Always make cell to be consistent
if(~iscell(csi_path))
    tmp{1} = csi_path; csi_path = tmp; clear tmp;
end

% If we get a list of all IMA files
if(numel(csi_path) > 1)
    csi_path_allfiles = csi_path;
    csi_first_file = csi_path{1};
    
% If we either get one IMA file or a folder
else
    
    % If we get one IMA file
    if(endsWith(csi_path{1},{'.IMA','.dcm'},'IgnoreCase',true))
        csi_first_file = csi_path{1};
        csi_path_allfiles = csi_path(1);
        
    % If we get a folder
    else

        csi_path_allfiles = dir( fullfile(csi_path{1},'*.*') );
        csi_path_allfiles = natsortfiles(csi_path_allfiles);
        csi_path_allfiles = {csi_path_allfiles.name}';
        DicomFoundVec = ~cellfun(@isempty,regexpi(csi_path_allfiles,'\.IMA|\.dcm'));
        csi_path_allfiles = csi_path_allfiles(DicomFoundVec);

        csi_path_allfiles = strcat(csi_path,'/',csi_path_allfiles);
        if(numel(csi_path_allfiles) == 0)
            fprintf('\nError in read_csi_dicom: Could not find any DICOM file. Abort.\n');
            MRStruct = struct;
            return;
        end
        csi_first_file = csi_path_allfiles{1};
    end
end


%% Parameters
% Read info from ascconv header
MRStruct.Par = read_ascconv(csi_first_file);
MRStruct.RecoPar = MRStruct.Par;
total_channel_no = MRStruct.Par.total_channel_no_reco;

% if(not(isnan(MRStruct.Par.nFreqEnc_FinalMatrix)))
%     ROW = MRStruct.Par.nFreqEnc_FinalMatrix;              % This gives the value FinalMatrixSize in the ascconv header, so the Matrix Size that is written to DICOM
% else                                    % file. In SVS these entries in the ascconv header do not exist.
    ROW = MRStruct.Par.nFreqEnc;                  
% end                                     
% if(not(isnan(MRStruct.Par.nPhasEnc_FinalMatrix)))
%     COL = MRStruct.Par.nPhasEnc_FinalMatrix;
% else
    COL = MRStruct.Par.nPhasEnc;
% end
% if(not(isnan(MRStruct.Par.nPhasEnc_FinalMatrix)))
%     SLC = MRStruct.Par.nSLC_FinalMatrix;
% else
%     SLC = MRStruct.Par.nPartEnc;
% end
SLC = numel(csi_path_allfiles);
if(ROW < 1)
    ROW = 1;
end
if(COL < 1)
    COL = 1;
end
if(MRStruct.Par.nPartEnc == 1)
    SLC = 1;
end
%SLC = MRStruct.Par.nSLC * SLC;  % Only correct for multislice data with nPartitions > 1 (never occurring!)
vecSize = MRStruct.Par.vecSize;



%% 1. READ & Reshape data

tic

MRStruct.Data = zeros(ROW,COL,SLC,vecSize);
MRStruct.Data = complex(MRStruct.Data,MRStruct.Data);


for CurSlc = 1:SLC

    % Read Data
    csi_fid = fopen(csi_path_allfiles{CurSlc},'r');
    % Find a certain bit pattern which indicates the end of the data
    % (For some reason, some of the VE spiral data sets have traling zeros at the end of the file
    % I dont understand the reason for that, but this code should work in both cases.)
    SearchBytesFromBack = 1024;
    FindBitPattern = [0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 0 0 1 0 0 0 0 1];  % 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 also needed?
    fseek(csi_fid, -SearchBytesFromBack, 'eof');
    SearchBitPattern = fread(csi_fid,'ubit1');
    FoundBitPattern_Ind = SearchBytesFromBack*8-strfind(SearchBitPattern',FindBitPattern);

    % Go to beginning of data: From the end of data (end - FoundBitPattern_Ind/8) back the amount of data we expect ( -((vecSize*ROW*COL*SLC*2*4)) )
    if(isempty(FoundBitPattern_Ind))
        pause(0.1)
        fseek(csi_fid, -((vecSize*ROW*COL*2*4)), 'eof');
    else
        fseek(csi_fid, -((vecSize*ROW*COL*2*4)) - FoundBitPattern_Ind/8, 'eof');   % /8 because I searched for bits, but fseek works with bytes
    end

    % Read data
    csi_data = fread(csi_fid,vecSize*ROW*COL*2,'float32');
    fclose(csi_fid);

    % Separate real & imaginary part
    csi_real=csi_data(1:2:end);                           % Odd points are real data
    csi_imag=csi_data(2:2:end);                           % Even imaginary
    csi_complex=complex(csi_real,csi_imag);

    % Replace NaNs by 0s
    csi_complex(isnan(csi_complex)) = 0;


    % total channel number not yet implemented, since ICE always creates 1 fiile per channel. there is anyway only one channel in the DICOM file.
    csi_complex2 = permute(reshape(csi_complex,[vecSize ROW COL]),[4 2 3 5 6 1]);
    MRStruct.Data(:,:,CurSlc,:) = csi_complex2;


end

fprintf('Read in took ... %f seconds.', toc)



%% 2. Spatial Shift

if(ne(Settings.x_shift,0))
    MRStruct.Data = circshift(MRStruct.Data, [Settings.x_shift 0 0 0]);
end

if(ne(Settings.y_shift,0))
    MRStruct.Data = circshift(MRStruct.Data, [0 Settings.y_shift 0 0]);
end



%% 3. Zero filling in kspace

if(ne(Settings.zerofilling_fact,1))
    
    tic
    
    % Go to k-space
    csi_kspace = ifftshift(ifftshift(MRStruct.Data,2),3);
    %csi_kspace = conj(csi_kspace);
    csi_kspace = ifft(ifft(csi_kspace,[],2),[],3);
    csi_kspace = fftshift(fftshift(csi_kspace,2),3);

    
    % Initialize zerofilled matrix
    csi_kspace_zf = zeros([total_channel_no ROW*Settings.zerofilling_fact COL*Settings.zerofilling_fact SLC vecSize]);
    
    % Write old k-space data in center of new k-space
    csi_kspace_zf(:,ROW*Settings.zerofilling_fact/4+1:ROW*Settings.zerofilling_fact*3/4,COL*Settings.zerofilling_fact/4+1:COL*Settings.zerofilling_fact*3/4,:,:) = csi_kspace;
    
    % overwrite variable name, clear other variable name
    %clear csi_kspace
    csi_kspace = csi_kspace_zf;
    clear csi_kspace_zf;
    
    
    % Go back to image space
    MRStruct.Data = ifftshift(ifftshift(csi_kspace,2),3);
    MRStruct.Data = fft(fft(MRStruct.Data,[],2),[],3);
    %MRStruct.Data = conj(MRStruct.Data);
    MRStruct.Data = fftshift(fftshift(MRStruct.Data,2),3);    
    
    
    fprintf('\nZerofilling took ... %f seconds', toc)
    
end


%%

if(Settings.RemoveFIDZerofilling_flag)
    RealVecSize = squeeze(sum(sum(sum(MRStruct.Data,1),2),3));
    RealVecSize = sum(abs(RealVecSize(:,1,1,1,1)) > 0);
    MRStruct.Data = MRStruct.Data(:,:,:,1:RealVecSize,:,:,:,:);
    MRStruct.RecoPar.vecSize = RealVecSize;
end



%% 4. Fourier transform to k-space if not done yet & if necessary

if(Settings.FTTokSpace_flag)
    
    if(~exist('csi_kspace','var'))
        tic

        MRStruct.Data = ifftshift(ifftshift(MRStruct.Data,2),3);
        MRStruct.Data = conj(MRStruct.Data);
        MRStruct.Data = ifft(ifft(MRStruct.Data,[],2),[],3);
        MRStruct.Data = fftshift(fftshift(MRStruct.Data,2),3);

        fprintf('\nIFFT took    ... %f seconds\n',toc)
    else
        MRStruct.Data = csi_kspace;
    end

end






%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);


