function ErrorOccurred = io_WriteMRSNiftiFile(MRStruct,NiftiFileToWrite,RealOrAbsFuncHandle)
% io_WriteNiftiFile Write an Array to minc file
%
% This function was written by Bernhard Strasser, September 2020.
%
%
% The function writes an array InArry to minc file NiftiFile using as header info either
% another minc-file or some parameters. It basically writes a raw file, and then converts this
% with rawtominc either with the -like option if a path to a minc file is given, or using the
% given parameters.
%
%
% ErrorOccurred = io_WriteNiftiFile(MRStruct,NiftiFile,LikeNiftiFileOrMincPar)
%
% Input: 
% -         MRStruct                       ...     Data to be written.
% -         NiftiFile                      ...     Path of minc file to which MRStruct should be written.
% -         LikeNiftiFileOrMincPar         ...     Path of other minc-file with same header, or Parameters containing
%                                                 infor for minc-header.
%
%
% Output:
% -        ErrorOccurred                          ...     1 if error occurred, otherwise 0.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None








%% 0. Preparations

Compression_flag = false;

if(endsWith(NiftiFileToWrite,'gz'))
    NiftiFileToWrite = regexprep(NiftiFileToWrite,'.nii.gz','.nii');
    Compression_flag = true;
end

if(isfield(MRStruct,'RecoPar'))
    CurPar = MRStruct.RecoPar;
else
    if(isfield(MRStruct,'Par'))
        CurPar = MRStruct.RecoPar;
    else
        error('No .Par or .RecoPar found in input struct.')
    end
end
if(~exist('RealOrAbsFuncHandle','var'))
    RealOrAbsFuncHandle = @real;
end


%%

    
   
MRStruct.Data=flip(flip(MRStruct.Data,1),2);
ppm_raw = compute_chemshift_vector(MRStruct);

SpecMap_raw = single(feval(RealOrAbsFuncHandle,fftshift(fft(MRStruct.Data,[],4),4)));

NiftiData = make_nii(SpecMap_raw);
NiftiData.hdr.dime.pixdim = [0 CurPar.FoV_Phase/size(SpecMap_raw,1) CurPar.FoV_Read/size(SpecMap_raw,2) CurPar.FoV_Partition/size(SpecMap_raw,3) ppm_raw(2)-ppm_raw(1) 1 1 1 1];

% Header_4d_raw = Header_3d;
% Header_4d_raw.ImageSize = size(SpecMap_raw);
% Header_4d_raw.PixelDimensions = [Header_4d_raw.PixelDimensions ppm_raw(2)-ppm_raw(1)];
NiftiData.hdr.dime.toffset = ppm_raw(1);
% Header_4d_raw.DisplayIntensityRange = [0 0.01];

% niftiwrite(SpecMap_raw,NiftiFileToWrite,Header_4d_raw,'compressed',Compression_flag)


save_nii(NiftiData,NiftiFileToWrite);


%%
if(Compression_flag)
    ErrorOccurred = unix(['gzip ' NiftiFileToWrite]);
end







