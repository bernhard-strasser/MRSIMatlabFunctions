function Data_Clean = GaussSeidelDecontamination(Data_Polluted, PSF, NoOfIter,NonNegativityConstrain_Flag)
% ExponentialFilter Apply an Exponential filter to time-domain signals (e.g. FIDs)
%
%
%
%
% The function computes an exponential filter in Hertz
%
%
% [OutArray,exp_filter_funct] = ExponentialFilter(InArray,dwelltime,ApplyAlongDim,exp_filter_Hz)
%
% Input: 
% -         Data_Polluted                     ...    Polluted Image
% -         PSF                   ...    PointSpreadFunction
% -         B0
%
% Output:
% -         OutArray                    ...     bla
%
%
% File dependancy:

% Further remarks: 
% This was copied from: http://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method#Program_to_solve_arbitrary_no._of_equations_using_Matlab


%% 0. Declarations, Preparations, Definitions

if(~exist('NonNegativityConstrain_Flag','var'))
	NonNegativityConstrain_Flag = 0;
end


%% 1.


N = size(Data_Polluted,1);
M = size(Data_Polluted,2);
Data_Clean = Data_Polluted;


for k = 1:NoOfIter

	PSF_Cur = circshift(PSF,[-N/2-1 -M/2-1]);
	
	for i = 1:N
		PSF_Cur = circshift(PSF_Cur,[1 0]);
		for j = 1:M
			PSF_Cur = circshift(PSF_Cur,[0 1]);
			
			% For all smaller indices
			CorrLow = 0;
			for i_bar = 1:i-1
				for j_bar = 1:j-1
					CorrLow = CorrLow + PSF_Cur(i_bar,j_bar)*Data_Clean(i_bar,j_bar);
				end
			end
			% For all bigger indices
			CorrHigh = 0;
			for i_bar = i:N
				for j_bar = j:M
					CorrHigh = CorrHigh + PSF_Cur(i_bar,j_bar)*Data_Clean(i_bar,j_bar);					
				end
			end			
			
			Data_Clean(i,j) = Data_Polluted(i,j) - CorrLow - CorrHigh;
			
		end
	end
	

	if(NonNegativityConstrain_Flag)
		Data_Clean(Data_Clean < 0) = 0;
	end
	
	
	
end


