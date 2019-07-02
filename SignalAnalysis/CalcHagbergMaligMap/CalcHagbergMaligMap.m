function [HagbergMaligMap1,HagbergMaligMap2] = CalcHagbergMaligMap(MetData_Ratio,MetData_Ratio_title)
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

HagbergMaligMap1 = 0;
HagbergMaligMap2 = 0;

if(nargin < 2)
	fprintf('\nError: Too little input. Aborting.\n')
	return
end

% Find out if it is Ratio to tCr or to Cr or to PCr
RatTo = regexpi(MetData_Ratio_title{1},'RatTo.*','match');
if(isempty(RatTo))
	fprintf('\nError: Are you sure that you gave me ratio maps? Aborting.\n')
	return
end

if(isempty(regexpi(RatTo{1},'Cr|Cr+PCr|PCr|tCr')))
	fprintf('\nError: Ratio maps seem to be not to Cr. Aborting.\n')
	return
end


% 0.3 Definitions
    





%% 1. Search for tNAA/Cr, Glx/Cr, MM/Cr, Cho/Cr, Ins/Cr in MetData_Ratio_title




SearchFor = {'(NAA\+NAAG|tNAA)_','(Glu\+Gln|Glx)_','mm_','(GPC\+PCh\+Cho|GPC\+PCh|GPC\+Cho|PCh\+Cho|tCh)_','(Ins|mI)_'};
SearchFor = strcat(SearchFor,regexprep(RatTo,'+','\\+'));

MetaboInd = zeros([1 numel(SearchFor)]);
for SearchInd = 1:numel(SearchFor)
	
	dummy = ~cellfun('isempty',regexp(MetData_Ratio_title,SearchFor{SearchInd}));
	if(sum(dummy) == 1)
		MetaboInd(SearchInd) = find(dummy);
	else
		MetaboInd(SearchInd) = 0;
	end		
	
end

% Check if all exists
if(sum(MetaboInd([1 2 4 5]) == 0) > 0)


	fprintf('\nError: Did not find all necessary metabolic ratio maps. Aborting.\n')
	return
end



%% 2. Compute Maps

Eigenvector1 = [0.4232, -0.3116, -0.3116, -0.6605, -0.5348];
Eigenvector2 = [0.4422, -0.3909, -0.3909,  0.7321, -0.3317];



for i = 1:numel(MetaboInd)
	
	if(MetaboInd(i))
		HagbergMaligMap1 = HagbergMaligMap1 + Eigenvector1(i)*MetData_Ratio(:,:,1,MetaboInd(i));
		HagbergMaligMap2 = HagbergMaligMap2 + Eigenvector2(i)*MetData_Ratio(:,:,1,MetaboInd(i));		
	end
	
end




%% 3. Postparations

% fclose(fid)






