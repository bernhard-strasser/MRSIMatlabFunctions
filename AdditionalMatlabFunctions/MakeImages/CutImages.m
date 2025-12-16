function [CellOfFigures] = CutImages(FilePaths,ImageEdges,PixelsOrPercent)
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
% [CellOfFigures] = CutImages(FilePaths,ImageEdges)
%
% Input: 
% -         FilePaths                    ...     FilePaths which should be read in ---- OR ---- cell of figures.
% -         ImageEdges                   ...     The image edges [up down left right] so the images will be cut to Image(up:down,left:right).
%                                                Can be also specified in percent (see next input), which will be then given in 
%                                                [up_percent down_percent left_percent right_percent]
% -         PixelsOrPercent              ...     If 'Percent', then the ImageEdges are interpreted as percent. Values of [100 100 100 100] will then not cut the
%                                                images at all, [50 100 50 100] will cut the upper side of the image by half, while not cutting from down
%                                                (resulting in an image which is 75% in size in up-down direction, but cutting only occurs from top), and it will
%                                                cut 50 % from left, but nothing from right.
% Output:
% -         CellOfFigures                ...     The cell of all the figures, which can then combined to one image by StitchImages.

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

if(~exist('ImageEdges','var') && ~exist('PixelsOrPercent','var'))
    PixelsOrPercent = 'Percent';
    ImageEdges = 100;
end

if(~exist('PixelsOrPercent','var'))
    PixelsOrPercent = 'Pixels';
end

if(~exist('AllEdgesSame_flag','var'))
    AllEdgesSame_flag = false;
end

if(iscell(ImageEdges) && numel(ImageEdges) == 1)
    ImageEdges = repmat(ImageEdges,size(FilePaths));
elseif(~iscell(ImageEdges))
    ImageEdges = repmat({ImageEdges},size(FilePaths));
end


%%

if(ischar(FilePaths{1,1}))
    CellOfFigures = cell(size(FilePaths));
    for jj = 1:size(FilePaths,2)
        for ii = 1:size(FilePaths,1)
            [CellOfFigures{ii,jj}, bla2]  =  imread( FilePaths{ii,jj} );
        end
    end
else
    CellOfFigures = FilePaths;
end


%% Convert from Percent To Pixels, if necessary

if(strcmpi(PixelsOrPercent,'Percent'))
    for jj = 1:size(CellOfFigures,2)
        for ii = 1:size(CellOfFigures,1)
            CurSiz = size(CellOfFigures{ii,jj});
            CurSiz = [CurSiz(1),CurSiz(1),CurSiz(2),CurSiz(2)]/2;
            ImageEdges{ii,jj} = round([1 0 1 0] + CurSiz + [-1, 1, -1, 1].*CurSiz .* ImageEdges{ii,jj}/100);
        end
    end
end

%% 1.


for jj = 1:size(FilePaths,2)
    for ii = 1:size(FilePaths,1)

        if((~iscell(ImageEdges) == 1 && ischar(ImageEdges) && strcmpi(ImageEdges,'all') ) ...
        || (iscell(ImageEdges) && ischar(ImageEdges{ii,jj}) && strcmpi(ImageEdges{ii,jj},'all'))) 
%             CellOfFigures{ii,jj} = bla;    
        else
            CellOfFigures{ii,jj} = CellOfFigures{ii,jj}(ImageEdges{ii,jj}(1):ImageEdges{ii,jj}(2),ImageEdges{ii,jj}(3):ImageEdges{ii,jj}(4),:);
        end
    end
end


