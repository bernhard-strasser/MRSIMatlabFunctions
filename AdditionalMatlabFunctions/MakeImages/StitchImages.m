function [OutImage,CellOfImages] = StitchImages(CellOfImages,SameExactSize,FillWithRGB)
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
% [OutImage,CellOfImages] = StitchImages(CellOfImages,SameExactSize)
%
% Input: 
% -         ?                    ...     ?
% Output:
% -         ?                         ...     ?

% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.





%% 0. Preparations


if(nargin < 1)
	fprintf('\nProblem: Cannot fft non-existant data. Provide some data.')
	FigHandle = 0;
	return;
end
if(~exist('FillWithRGB','var'))
    FillWithRGB = [102 102 102]; % some grey
end
if(exist('SameExactSize','var') && isempty(SameExactSize))
    clear SameExactSize;
end



%% 1.

CellSize = size(CellOfImages);

if(~exist('SameExactSize','var'))
    Sizes = cellfun(@size,CellOfImages, 'UniformOutput',false);
    Sizes = permute(reshape(cell2mat(Sizes),[CellSize(1) 3 CellSize(2)]),[1 3 2]);
    Sizes = Sizes(:,:,1:2);     % First entry in 3rd dim: height, Second entry: width of image

    TargetWidths = max(Sizes(:,:,2),[],1);
    TargetHeights = max(Sizes(:,:,1),[],2);
else
    TargetWidths = repmat(SameExactSize(2),[1 size(CellOfImages,2)]);
    TargetHeights = repmat(SameExactSize(1),[1 size(CellOfImages,1)]);
end

% Each Row has to have the same size in the col-dimension
% Each Column has to have the same size in the row-dimension
for CurRow = 1:size(CellOfImages,1)
    for CurCol = 1:size(CellOfImages,2)
        TargetSize(CurRow,CurCol,:) = [TargetHeights(CurRow); TargetWidths(CurCol); 3];
        
        % Fill figures with background
        % Example: Zerofill image of size [3 6] to [12 14].
        % Will place voxels [12/2+1 - 1 : 12/2+1 + 1, 14/2+1 - 2 : 14/2+1 + 3]
        % General: [floor(12/2)+1 - floor(3/2) : floor(12/2)+1 + floor(3/2), floor(14/2)+1 - floor(6/2) : floor(12/2)+1 + floor(3/2)]
        CellOfImages{CurRow,CurCol} = ...
        XFillOrCutData(CellOfImages{CurRow,CurCol},transpose(squeeze(TargetSize(CurRow,CurCol,:))),0,FillWithRGB);
        
    end
end
OutImage = cell2mat(CellOfImages);





