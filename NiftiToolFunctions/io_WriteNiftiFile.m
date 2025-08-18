function ErrorOccurred = io_WriteNiftiFile(MRStruct,NiftiFile,LikeMincFileOrMincPar)
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

if(~exist('LikeMincFileOrMincPar','var'))
    LikeMincFileOrMincPar = '';
end

if(~isstring(LikeMincFileOrMincPar) && (endsWith(LikeMincFileOrMincPar,'nii') || endsWith(LikeMincFileOrMincPar,'nii.gz')))
    error('Error: Cannot process nifti files as like')
end

Compression_flag = false;

if(endsWith(NiftiFile,'gz'))
    NiftiFile = regexprep(NiftiFile,'.nii.gz','.nii');
    Compression_flag = true;
end


%% Write Minc File


temp =  java.util.UUID.randomUUID;
myuuid = char(temp.toString);
MincFile = regexprep(NiftiFile,'.nii',['_' myuuid '.mnc']);
ErrorOccurred = io_WriteMincFile(MRStruct,MincFile,LikeMincFileOrMincPar);




%% Convert Nii to Minc

if(~ErrorOccurred)
    ErrorOccurred = unix(['mnc2nii ' MincFile ' ' NiftiFile]);
    if(Compression_flag && ~ErrorOccurred)
        ErrorOccurred = unix(['gzip ' NiftiFile]);
    end
end

%% Remove temporary Minc File

delete(MincFile);


