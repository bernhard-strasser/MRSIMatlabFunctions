%% -1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                 PROGRAM TO HELL AND BACK                           %%%%%%%%%%%%%%
%%%%%%%%%%%%%%     ENTSTANDEN DURCH DIE KUNST DES PROGRAMMIERENS DURCH KONSEQUENTES ANSTARREN     %%%%%%%%%%%%%%
%%%%%%%%%%%%%%     Find out the best CAIPIRINHA pattern by computing the average distance of      %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                        measured to not measured k-space points.                    %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. DEFINITIONS, PREPARATIONS

cell_size = [4 4];      % e.g. undersampling of factor 2 in both directions, and to get more patterns double that size in both directions.
no_measured_points = 4; % So R = 4*4/4 = 4;




%% 1. Create All Possible Patterns.


Patterns = nchoosek(2:prod(cell_size),no_measured_points-1);    % nchoosek: Binomial combinations. Distribute 1 point less, because one can always choose the point (1,1) w.l.o.g.
Patterns = cat(2, ones(size(Patterns,1),1), Patterns);          % Add the point (1,1) to all Patterns.



%% 2. Loop over all Patterns, Compute the mean distance for each.


for Patt_no = 1:size(Patterns,1)
    
    
    
    
    
    
    
    
    
end




