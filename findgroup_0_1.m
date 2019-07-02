function Groups = findgroup_0_1(LogicalVector,MinGroupSize)
%
% 

%% 0. Declarations, Preparations, Definitions

% 0.1 Declarations


% 0.2 Definitions
    

% 0.3 Preparations

pause on




%%

%Indices = find(LogicalVector);

ExtendedLogicalVector = [0 LogicalVector 0];
DerivLogVec = diff(ExtendedLogicalVector);
Groups.StartPts = find(DerivLogVec == 1);
Groups.EndPts = find(DerivLogVec == -1) -1;

ExcludeGroups = Groups.EndPts - Groups.StartPts + 1 < MinGroupSize;
%ExcludeGroups = find(ExcludeGroups)

Groups.StartPts(ExcludeGroups) = [];
Groups.EndPts(ExcludeGroups) = [];



%% 7. THE END

pause off







