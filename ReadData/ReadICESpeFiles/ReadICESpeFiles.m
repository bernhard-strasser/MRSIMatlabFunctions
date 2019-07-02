function OutData = ReadICESpeFiles(InPath,MatrixSize_RO,MatrixSize_PE,MatrixSize_SL,MatrixSize_VecSize)
%
% ReadICESpeFiles Read in .spe-files as created by the spiral ICE.
%
% This function was written by Bernhard Strasser, Nov 2015.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% OutData = ReadICESpeFiles(InPath,MatrixSize_RO,MatrixSize_PE,MatrixSize_SL,MatrixSize_VecSize)
%
% Input: 
% -         InPath                      ...     Path ofspe file.
% -         MatrixSize_RO               ...     Size in ReadOut dimension
% -         MatrixSize_PE               ...     Size in Phase Encoding dimension
% -         MatrixSize_SL               ...     Size inslice dimension
% -         MatrixSize_VecSize          ...     Size in spectroscopic dimension
%
%
%
% Output:
% -         OutData                     ...     Output data
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: ,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.



%% 

if(~exist('InPath','var'))
    fprintf('\nError: No file was input.');
    return;
end
if(~exist(InPath,'file'))
    fprintf('\nError: Input file %s does not exist.', InPath);
    return;
end    
    
if(~exist('MatrixSize_RO','var')) 
    MatrixSize_RO = 16;
end
if(~exist('MatrixSize_PE','var')) 
    MatrixSize_PE = 16;
end
if(~exist('MatrixSize_SL','var')) 
    MatrixSize_SL = 1;
end
if(~exist('MatrixSize_VecSize','var')) 
    MatrixSize_VecSize = 512;
end




%%



% WriteToFile_0001.spe
fid = fopen(InPath);
OutData = fread(fid, [2,MatrixSize_RO*MatrixSize_PE*MatrixSize_SL*MatrixSize_VecSize], 'float32');
fclose(fid);

OutData = reshape(OutData,2,MatrixSize_VecSize,MatrixSize_RO,MatrixSize_PE,MatrixSize_SL);
OutData = squeeze( complex(OutData(1,:,:,:,:),OutData(2,:,:,:,:)));


