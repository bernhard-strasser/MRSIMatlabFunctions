function Groups = findgroup_1_1(LogicalVector,MinGroupSize)
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
NumberOfGroups = size(find(DerivLogVec == 1),2);
Groups = zeros(NumberOfGroups,2);

if(NumberOfGroups == 0)
    Groups = cell([1 0]);
    return
end
Groups(:,1) = find(DerivLogVec == 1);
Groups(:,2) = find(DerivLogVec == -1) -1;

ExcludeGroups = Groups(:,2) - Groups(:,1) + 1 < MinGroupSize;

Groups(ExcludeGroups,:) = [];
Groups = mat2cell(Groups,repmat(1,[1 size(find(~ExcludeGroups),1)]));
Groups = transpose(Groups);



%% 7. THE END

pause off







