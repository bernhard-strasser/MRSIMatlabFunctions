function [Figgy] = ReplaceTealColourInMincScreenshots(Figgy)
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
    st = dbstack; funcname = st.name;
	fprintf('\nProblem in %s: Not enough input.',funcname)
	Figgy = 0;
	return;
end




%% 1.

for ii = 1:size(Figgy,1)
    for jj = 1:size(Figgy,2)
        test = squeeze(Figgy(ii,jj,:));
        if(test(1) == test(2) && test(1) == test(3) && test(2) == test(3) && test(1) < 4)
            Figgy(ii,jj,:) = zeros([1 1 3]);
        end
    end
end




