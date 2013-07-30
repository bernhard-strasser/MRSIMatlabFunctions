%% -1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                 PROGRAM TO HELL AND BACK                           %%%%%%%%%%%%%%
%%%%%%%%%%%%%%     ENTSTANDEN DURCH DIE KUNST DES PROGRAMMIERENS DURCH KONSEQUENTES ANSTARREN     %%%%%%%%%%%%%%
%%%%%%%%%%%%%%     Find out the best CAIPIRINHA pattern by computing the average distance of      %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                        measured to not measured k-space points.                    %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. DEFINITIONS, PREPARATIONS

clear variables; clear functions; close all;
pause on


cell_size = [6 6];      % e.g. undersampling of factor 2 in both directions, and to get more patterns double that size in both directions.
no_measured_points = 4; % So R = 4*4/4 = 4;
R = prod(cell_size)/no_measured_points

no_Patterns = nchoosek(prod(cell_size)-1,no_measured_points-1);


%% 1. Create All Possible Patterns.


PtsMeas = nchoosek(2:prod(cell_size),no_measured_points-1);    % nchoosek: Binomial combinations. Distribute 1 point less, because one can always choose the point (1,1) w.l.o.g.
PtsMeas = cat(2, ones(size(PtsMeas,1),1), PtsMeas);          % Add the point (1,1) to all Patterns.



%% 2. Loop over all Patterns, Compute the mean distance for each.


QualityMeasure = zeros([size(PtsMeas,1) 3]);


for Patt_no = 1:size(PtsMeas,1)
    
    QualityMeasure(Patt_no,:) = kSpace_DistBtwMeasPts(squeeze(PtsMeas(Patt_no,:)), cell_size);
    
end




%% 3. Find the Best 10%

% Index_best10p_std = QualityMeasure(:,1) <= quantile(QualityMeasure(:,1),0.10);
% Index_best_std = QualityMeasure(:,1) == min(QualityMeasure(:,1));
% 
% Index_best10p_min = QualityMeasure(:,2) >= quantile(QualityMeasure(:,2),0.90);
% Index_best_min = QualityMeasure(:,2) == max(QualityMeasure(:,2));
% 
% Index_best10p = find(Index_best10p_std .* Index_best10p_min);
% Index_best = find(Index_best_std .* Index_best_min);


Index_best10p_comb = find(QualityMeasure(:,3) <= quantile(QualityMeasure(:,3),0.025));
Index_best_comb = find(QualityMeasure(:,3) == min(QualityMeasure(:,3)));



%% 4. Plot those



for Indexl = transpose(Index_best_comb)
    
    Indexl
    QualityMeasure(Indexl,:)
    
    ElemCell = zeros(cell_size);
    ElemCell(squeeze(PtsMeas(Indexl,:))) = 1;
    ElemCell_fig = figure; 
    imagesc(ElemCell)
    movegui(ElemCell_fig, 'northwest')
    
    
    ElemCell_recol = ElemCell; ElemCell_recol(ElemCell_recol == 0) = 2; ElemCell_recol(ElemCell_recol == 1) = 3;
    RepCell = repmat(ElemCell,[3 3]);
    RepCell(size(ElemCell,1)+1 : 2*size(ElemCell,1), size(ElemCell,2)+1 : 2*size(ElemCell,2)) = ElemCell_recol;
    RepCell_fig = figure; 
    imagesc(RepCell)
    movegui(RepCell_fig, 'north')    
    
    pause
    close all
    
end



for Indexl = transpose(Index_best10p_comb)
    
    Indexl
    QualityMeasure(Indexl,:)
    
    ElemCell = zeros(cell_size);
    ElemCell(squeeze(PtsMeas(Indexl,:))) = 1;
    ElemCell_fig = figure; 
    imagesc(ElemCell)
    movegui(ElemCell_fig, 'northwest')
    
    
    ElemCell_recol = ElemCell; ElemCell_recol(ElemCell_recol == 0) = 2; ElemCell_recol(ElemCell_recol == 1) = 3;
    RepCell = repmat(ElemCell,[3 3]);
    RepCell(size(ElemCell,1)+1 : 2*size(ElemCell,1), size(ElemCell,2)+1 : 2*size(ElemCell,2)) = ElemCell_recol;
    RepCell_fig = figure; 
    imagesc(RepCell)
    movegui(RepCell_fig, 'north')    
    
    pause
    close all    
        
end


pause off





