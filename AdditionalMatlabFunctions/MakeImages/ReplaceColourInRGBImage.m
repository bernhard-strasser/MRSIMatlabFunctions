function [Figgy] = ReplaceColourInRGBImage(Figgy,ReplaceColour_RGB,ColourSimilarity,ReplaceColourWith_RGB)
%
% ReplaceTealColourInMincScreenshots Read in csi-data from Siemens raw file format
%
% This function was written by Bernhard Strasser, March 2019.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [OutImage,CellOfImages] = StitchImages(CellOfImages,SameExactSize)
%
% Input: 
% -         Figgy                    ...     The Figure where the colours should be replaced
% -         ReplaceColour_RGB        ...     The RGB colour which should be replaced
% -         ColourSimilarity         ...     The similarity. If 0, only exactly this colour will be replaced. In general, 
%                                            norm(ActualColour - ReplaceColour_RGB) < ColourSimilarity. Keep in mind that RGB values are uint8, so each
%                                            channel can have values between 0 and 255.
% -         ReplaceColourWith_RGB    ...     The colour which will replace the matching colours.
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
    st = dbstack; funcname = st.name;
	fprintf('\nProblem in %s: Not enough input.',funcname)
	Figgy = 0;
	return;
end

if(~exist('ReplaceColour_RGB','var'))
    ReplaceColour_RGB = [1 1 1];
end
if(~exist('ReplaceColourWith_RGB','var'))
    ReplaceColourWith_RGB = [0 0 0];
end
if(~exist('ColourSimilarity','var'))
    ColourSimilarity = 0;
end

%% 1.

CellSz = 1;
if(iscell(Figgy))
    CellSz = numel(Figgy);
end
for CellInd = 1:CellSz
    if(iscell(Figgy))
        CurFiggy = Figgy{CellInd};
    else
        CurFiggy = Figgy;
    end
    
    
    for ii = 1:size(CurFiggy,1)
        for jj = 1:size(CurFiggy,2)
            test = transpose(double(squeeze(CurFiggy(ii,jj,:))));

            if(norm(test - ReplaceColour_RGB) <= ColourSimilarity)
                CurFiggy(ii,jj,:) = uint8(ReplaceColourWith_RGB);
            end

    %         if(test(1) == test(2) && test(1) == test(3) && test(2) == test(3) && test(1) < 4)
    %             CurFiggy(ii,jj,:) = zeros([1 1 3]);
    %         end
        end
    end
    
    
    if(iscell(Figgy))
        Figgy{CellInd} = CurFiggy;
    else
        Figgy = CurFiggy;
    end
end



