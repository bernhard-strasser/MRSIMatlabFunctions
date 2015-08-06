function [SynthMap] = CalcSynthMap(MetMap,MetMapTitle,SynthMapWeights,SynthMapTitle,Clip_flag)
%
% CalcHagbergMaligMap Calculate Malignancy Maps According to Hagberg et al., 1995, MRM 34(2):242-52
%
% This function was written by Bernhard Strasser, [month] [year].
%
%
% The function can really do nothing, and more specifically, exactly nothing.
% 
%
%
% [A,B] = read_csi_dat_1_10(inputvar1,inputvar2)
%
% Input: 
% -         inputvar1                   ...    This is the first input
% -         inputvar2                   ...    And this the second
%
% Output:
% -         A                           ...     This is output A
% -         B                           ...     This is output B
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations

SynthMap = 0;

if(nargin < 2)
	fprintf('\nError: Too little input. Aborting.\n')
	return
end

if(~exist('Clip_flag','var'))
	Clip_flag = true;
end

% 0.3 Definitions
    





%% 1. Search for SynthMapTitles in MetMapTitles


MetaboInd = regexp_SearchCellInCell(MetMapTitle,SynthMapTitle);
if( sum(cellfun(@numel,MetaboInd) > 1) == 0 )
	MetaboInd = cell2mat(MetaboInd);
else
	fprintf('\nUnexpected Error: At least one entry of SynthMapTitle occurs several times in MetMapTitle.\nAbort.')
	return
end

% Check if all exists
if(numel(MetaboInd) < numel(SynthMapTitle))
	fprintf('\nError: Did not find all necessary metabolic maps. Aborting.\n')
	return
end



%% 2. Compute Map

MetMap_small = MetMap(:,:,:,MetaboInd);
SynthMapWeights_mat = myrepmat(SynthMapWeights,size(MetMap_small));

SynthMap = sum(SynthMapWeights_mat .* MetMap_small,4);




%% 3. Clip Maps

if(Clip_flag)   
    SynthMap_MedMad1 = transpose(nanmedian_own(SynthMap(:)) + 6*mad_own(SynthMap(:),1)); 
    SynthMap_MedMad2 = transpose(nanmedian_own(SynthMap(:)) - 6*mad_own(SynthMap(:),1)); 
    SynthMap(SynthMap > repmat(SynthMap_MedMad1,size(SynthMap))) = NaN;
	SynthMap(SynthMap < repmat(SynthMap_MedMad2,size(SynthMap))) = NaN;
end


%% 4. Postparations

% fclose(fid)






