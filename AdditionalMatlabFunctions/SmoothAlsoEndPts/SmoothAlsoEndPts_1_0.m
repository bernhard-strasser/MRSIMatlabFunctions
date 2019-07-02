function SmoothedVec = SmoothAlsoEndPts_1_0(Vec,Span)
%
% 

%% 0. Declarations, Preparations, Definitions

% 0.1 Declarations


% 0.2 Definitions
    

% 0.3 Preparations
if(Span <= 0)
    Span = 1;
end
if(Span > numel(Vec))
    Span = numel(Vec);
end

if(mod(Span,2) == 0)
    Span = Span - 1;
end


pause on




%%


SmoothedVec = smooth(Vec,Span);
SmoothedVec(1:(Span+1)/2) = mean(Vec(1:Span));
SmoothedVec(end-(Span-1)/2 : end) = mean(Vec(end - Span + 1 : end));





%% 7. THE END

pause off







