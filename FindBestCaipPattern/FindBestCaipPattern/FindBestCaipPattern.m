function [BestPatterns,BestXPercPatterns,no_Patterns,QualityMeasure, AllPatterns] = FindBestCaipPattern(cell_size, no_measured_points, no_WantedPatterns)
%
% kSpace_Distance Compute quality measure of kSpace Pattern.
%
% This function was written by Bernhard Strasser, July 2013.
%
%
% This function 1) calculates the distance of a non-measured kSpace point to all the measured ones.
%               2) Does this for all non-measured points.
%               3) Calculates a quality measure for those distances (like maximum, mean or combination).
%               4) Does this for all non-measured points.
%               5) Computes the same quality measure using all non-measured points (maximum, mean or combination).
%
%
% [csi,csi_kspace] = read_csi_1_4(csi_path, zerofill_to_nextpow2_flag, zerofilling_fact, x_shift,y_shift)
%
% Input: 
% -         PtsMeas                     ...     The measured points as linear index (e.g. if A=[1 0;0 1], 1...measured --> PtsMeas = [1 4])
% -         CellSize                    ...     The size of the elementary caipirinha cell. in the above example: CellSize = [2 2]
%                                               Use at least 2*[Rx,Ry], where Rx and Ry would be the GRAPPA acceleration factors in x and y dir.
%                                               Otherwise you get only 1 possible undersampling cell.
%
% Output:
% -         QualityMeasure              ...     The quality measure of the kSpace Pattern, e.g. the mean distance of non-measured to measured points,
%                                               or the maximum of the distances between measured and non-measured ones.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: none yet.




%% 0. DEFINITIONS, PREPARATIONS

% Initialize
BestPatterns.Indices = 0;
BestPatterns.Patterns = 0;
BestXPercPatterns.Indices = 0;
BestXPercPatterns.Patterns = 0;
QualityMeasure = 0;
no_Patterns = 0;


% Assign standard values to variables if nothing is passed to function.
if(nargin < 2)
    display([ char(10) 'Gimme more input, Ma''am!' char(10) ])
    return;
end 
if(~exist('no_WantedPatterns','var'))
    no_WantedPatterns = 0;
end

if(no_measured_points == 1)
	BestPatterns.Indices = 1;
	BestPatterns.Patterns = 1;
	BestPatterns.CellSizes = cell_size;
	BestXPercPatterns = BestPatterns;
	no_Patterns = 1;
    QualityMeasure = DistBtwMeasPts(1, cell_size,0);
    AllPatterns = 1;
    return
end




% Check if too many patterns are possible
no_Patterns = nchoosek(prod(cell_size)-1,no_measured_points-1);
MatrixSize_GBytes = 3*no_Patterns*(no_measured_points-1)*1/2^30;		% *1: uin8 use 8 bits = 1 byte; 3*: security measure.
fprintf('\n\n%d possible patterns.', no_Patterns)

try
	[memused,memfree] = memused_linux(1);
catch
	memfree = 8000;
end

if(MatrixSize_GBytes > memfree/2^10)
    fprintf('\nSorry, there is not enough memory free.\n')
    return
end


if(no_Patterns > 10^8)
    fprintf('\nThese are too many. I will quit here.')
    return
elseif(no_Patterns > 10^7)
    fprintf('\nThis may take really long. Overnight processing recommended!')    
elseif(no_Patterns > 10^6)
    fprintf('\nThis may take a while. Sit back, relax & drink a tea!')
else
    fprintf('\n')    
end

if(no_WantedPatterns <= 0)
    if(no_Patterns < 70)
        PercentageBestPatterns = 0.51;
    elseif(no_Patterns < 2000)
        PercentageBestPatterns = 0.1;
    else
        PercentageBestPatterns = 200/no_Patterns;
    end
else
    PercentageBestPatterns = no_WantedPatterns/no_Patterns;
end
PercentageBestPatterns(PercentageBestPatterns > 1) = 1;

fprintf('\n Use %10.8f %% for the BestXPercPatterns.\n', PercentageBestPatterns*100)    







%% 1. Create All Possible Patterns of the Measured Points.

if(prod(cell_size) > 255)
	fprintf('\nProblem in FindBestCaipPattern.m occurred: AllPatterns can encode numbers 0-255.\nHowever higher numbers are needed in AllPatterns. Abort')
	return
end

AllPatterns = uint8(nchoosek(uint8(2:prod(cell_size)),uint8(no_measured_points-1)));    % nchoosek: Binomial combinations. Distribute 1 point less, because one can always choose the point (1,1) w.l.o.g.

% if(no_Patterns > 10^9)
% 	AllPatterns2 = zeros([size(AllPatterns,1) size(AllPatterns,2)+1],'int16');
% 	for PattNo = 1:size(AllPatterns,1)
% 		AllPatterns2(PattNo,:) = cat(2,1,AllPatterns(PattNo,:));
% 	end
% 	clear AllPatterns; AllPatterns = AllPatterns2; clear AllPatterns2;
% else
	AllPatterns = cat(2, ones(size(AllPatterns,1),1,'uint8'), AllPatterns);          % Add the point (1,1) to all Patterns.
% end



%% 2. Loop over all Patterns, Compute the mean distance for each.

QltyMeas_dummy = kSpace_DistBtwMeasAndNonMeasPts2(squeeze(AllPatterns(1,:)), cell_size,1);
NoQltyMeas = size(QltyMeas_dummy,2);

QualityMeasure = zeros([size(AllPatterns,1) NoQltyMeas],'single');
clear QltyMeas_dummy NoQltyMeas

if(~logical(matlabpool('size')) && size(AllPatterns,1) > 10000)
	matlabpool 7;
end

ParallelEnabled = logical(matlabpool('size'));
if(ParallelEnabled)
	parfor Patt_no = 1:size(AllPatterns,1)    
		QualityMeasure(Patt_no,:) = single(kSpace_DistBtwMeasAndNonMeasPts2(squeeze(AllPatterns(Patt_no,:)), cell_size,1));
	end	
	matlabpool close;
else
	for Patt_no = 1:size(AllPatterns,1)    
		QualityMeasure(Patt_no,:) = single(kSpace_DistBtwMeasAndNonMeasPts2(squeeze(AllPatterns(Patt_no,:)), cell_size,1));
	end
end






%% 3. Find the Best & the Best 10% Patterns


% In Image Domain
% Index_BestXPerc = QualityMeasure(:,2) >= quantile(QualityMeasure(:,2),PercentageBestPatterns);
% Index_Best = QualityMeasure(:,2) == max(QualityMeasure(:,2));

% In kSpace domain
Index_BestXPerc = QualityMeasure(:,end) <= quantile(QualityMeasure(:,end),PercentageBestPatterns);
Index_Best = QualityMeasure(:,end) == min(QualityMeasure(:,end));

BestPatterns.Indices = find(Index_Best);
BestXPercPatterns.Indices = find(Index_BestXPerc);

Index_BestXPerc = repmat(Index_BestXPerc, [1 no_measured_points]);
Index_Best = repmat(Index_Best, [1 no_measured_points]);


BestPatterns.Patterns = AllPatterns(Index_Best);
BestXPercPatterns.Patterns = AllPatterns(Index_BestXPerc);
if(nargout < 5)
	clear AllPatterns;
end
	
BestPatterns.Patterns = reshape(BestPatterns.Patterns, [numel(BestPatterns.Patterns)/no_measured_points no_measured_points]);
BestXPercPatterns.Patterns = reshape(BestXPercPatterns.Patterns, [numel(BestXPercPatterns.Patterns)/no_measured_points no_measured_points]);

BestPatterns.Patterns = transpose(BestPatterns.Patterns);
BestXPercPatterns.Patterns = transpose(BestXPercPatterns.Patterns);


BestPatterns.CellSizes = transpose(squeeze(repmat(cell_size,[1 1 size(BestPatterns.Patterns,2)])));
BestXPercPatterns.CellSizes = transpose(squeeze(repmat(cell_size,[1 1 size(BestXPercPatterns.Patterns,2)])));






%% 4. Exclude similar patterns

% Not yet done -- Is it really necessary?
BestPatterns = ExcludeCaipiPatternsByCircshift(BestPatterns,false);

if(numel(BestXPercPatterns.Indices) < 500)
	BestXPercPatterns = ExcludeCaipiPatternsByCircshift(BestXPercPatterns,false);
end












end







%% L. Local Functions


function PatternOut = ExcludeCaipiPatternsByCircshift(PatternIn,DebugMode)




	PatternIndSimilar = zeros(numel(PatternIn.Indices));
	FoundInd = 0;
	for PatternNo1 = 1:numel(PatternIn.Indices)
		CompareMat1 = zeros(squeeze(PatternIn.CellSizes(PatternNo1,:)));
		CompareMat1(PatternIn.Patterns(:,PatternNo1)) = 1;
		if(ismember(PatternNo1,PatternIndSimilar))
			continue
		end
		
		for PatternNo2 = PatternNo1+1:numel(PatternIn.Indices)
			if(ismember(PatternNo2,PatternIndSimilar))
				continue
			end
			BreakOut = false;

			
			CompareMat2_dummy = zeros(squeeze(PatternIn.CellSizes(PatternNo2,:)));
			CompareMat2_dummy(PatternIn.Patterns(:,PatternNo2)) = 1;
			
			for xshift = 1:PatternIn.CellSizes(PatternNo1,1)
				for yshift = 1:PatternIn.CellSizes(PatternNo2,2)
					
					CompareMat2 = circshift(CompareMat2_dummy,[xshift yshift]);
					
					if(DebugMode)
						fprintf(  '\nSimilar %d: Patt %d and %d with [xshift,yshift] = [%d,%d]', ~ismember(0,CompareMat1 == CompareMat2), PatternNo1, PatternNo2,xshift,yshift  )
						PatternIndSimilar(~PatternIndSimilar == 0)
						figure; imagesc(CompareMat1), movegui(gcf,'northwest')
						figure; imagesc(CompareMat2), movegui(gcf,'northeast')
						waitforbuttonpress
						close all;
					end
					
					if(~ismember(0,CompareMat1 == CompareMat2))
						FoundInd = FoundInd + 1;
						PatternIndSimilar(FoundInd) = PatternNo2;
						BreakOut = true;
						break;
					end
					
					
					
					
					
					
				end
				if(BreakOut), break; end
			end
		end
	end
	PatternIndSimilar(PatternIndSimilar == 0) = [];
	
	PatternOut = PatternIn; PatternOut.Indices(PatternIndSimilar) = []; PatternOut.Patterns(:,PatternIndSimilar) = []; PatternOut.CellSizes(PatternIndSimilar,:) = [];
	
	
	
	
	
	

end
