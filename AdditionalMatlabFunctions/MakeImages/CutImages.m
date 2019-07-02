function [CellOfFigures] = CutImages(FilePaths,ImageEdges)
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


