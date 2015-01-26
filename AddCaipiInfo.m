function Info_pooled_updated = AddCaipiInfo(InfoPath,Info_pooled)
%
% AddCaipiInfo Add Info of FindBest21DCaipPattern_ByArtifactPower-scripts to existing ones.
%
% This function was written by Bernhard Strasser, Jannuary 2015.
%
%
% The function takes the Path to some Info-output of a FindBest21DCaipPattern_ByArtifactPower-script,
% checks if there is something that is not already existing, and adds the additional info if there is.
%
%
% Info_pooled_updated = AddCaipiInfo(InfoPath,Info_pooled)
%
% Input: 
% -         InfoPath                               ...    Path to the Info which should be added.
% -         Info_pooled                            ...    The already existing Info.
%
% Output:
% -         Info_pooled_updated                    ...    The updated pooled Info.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations

if(~exist('InfoPath','var'))
    fprintf('\nMissing InfoPath.')
    Info_pooled_updated = Info_pooled;
    return
end
if(~exist(InfoPath,'file'))
    fprintf('\nFile\n%s\ndoes not exist.',InfoPath)
    Info_pooled_updated = Info_pooled;
    return
end



% 0.3 Definitions
    






%% 1. Load Info

try
    load(InfoPath)
catch err
    fprintf('\nCould not load\n%s\ndue to the following error:\n%s',InfoPath,err.message)
    Info_pooled_updated = Info_pooled;
    return;
end



%% 2. Find out volunteer in InfoPath, search for that Volunteer in Info_pooled

Vol = regexp(regexpi(InfoPath,'Vol_\d\d|Vol\d\d|Proband_\d\d|Proband\d\d','match'),'\d\d','match');
Vol = Vol{1}{:};

VolsPooled = fields(Info_pooled);
VolInPooled = ~cellfun('isempty',regexp(VolsPooled,Vol));

if(sum(VolInPooled) > 0)
    VolInPooled = VolsPooled{VolInPooled};
else
    VolInPooled = ['Vol' Vol];
end



%% 3. Add Info

Rtotals = fields(Info);
Info_pooled_updated = Info_pooled;

for RTotal = transpose(Rtotals)
    for Slcflag = 1:2
        WholeMatrixAddedFlag = false;
        Slcstr = ['RSlc_' num2str(Slcflag)];
        
        % If the RSlc doesnt exist in the Info, we dont need to add it
        if(~isfield(Info.(RTotal{:}),Slcstr))
            continue
        end
        
        
        for RInPlaneRealiz = transpose(fields(Info.(RTotal{:}).(Slcstr)))
            RInPlaneRealizStr = RInPlaneRealiz{:};

            
            
            % Add Whole Slices If Possible
            if( (~isfield(Info_pooled_updated.(VolInPooled).(RTotal{:}),Slcstr) || ...
               (isfield(Info_pooled_updated.(VolInPooled).(RTotal{:}),Slcstr) && isempty(Info_pooled_updated.(VolInPooled).(RTotal{:}).(Slcstr)) ) ) && ...
               sum(Info_pooled_updated.(VolInPooled).(RTotal{:}).WholeMatrix(:,4) == Slcflag) == 0 )

                fprintf('\nAdding field %s.%s.%s to Info_pooled_updated.',VolInPooled,RTotal{:},Slcstr);
                Info_pooled_updated.(VolInPooled).(RTotal{:}).(Slcstr) = Info.(RTotal{:}).(Slcstr);
                fprintf('\nAdding entries with RSlc = %f in Info_pooled_updated.%s.%s.WholeMatrix to Info_pooled_updated.',Slcflag,VolInPooled,RTotal{:});
                Info_pooled_updated.(VolInPooled).(RTotal{:}).WholeMatrix = cat(1,Info_pooled_updated.(VolInPooled).(RTotal{:}).WholeMatrix, Info.(RTotal{:}).WholeMatrix);        
                break
            end
            
            
            
            % Add all RInPlaneRealiz
            if(isfield(Info_pooled_updated.(VolInPooled).(RTotal{:}).(Slcstr),RInPlaneRealizStr) && ~isempty(Info_pooled_updated.(VolInPooled).(RTotal{:}).(Slcstr).(RInPlaneRealizStr)) ) % If exists and has elements, continue
                fprintf('\n''Info_pooled_updated'' contains already field %s.%s.%s.%s. Skip adding.',VolInPooled,RTotal{:},Slcstr,RInPlaneRealizStr)
            else
                 fprintf('\nAdding field %s.%s.%s.%s to Info_pooled_updated.',VolInPooled,RTotal{:},Slcstr,RInPlaneRealizStr);
                 Info_pooled_updated.(VolInPooled).(RTotal{:}).(Slcstr).(RInPlaneRealizStr) = Info.(RTotal{:}).(Slcstr).(RInPlaneRealizStr);               
            end
            
            
            % Check if Pooled WholeMatrix contains already the Matrix which should be added
            if( ~WholeMatrixAddedFlag && (~isfield(Info_pooled_updated.(VolInPooled).(RTotal{:}), 'WholeMatrix') || isempty(Info_pooled_updated.(VolInPooled).(RTotal{:}).WholeMatrix) || ...
            sum(ismember(Info_pooled_updated.(VolInPooled).(RTotal{:}).WholeMatrix(:,7), Info.(RTotal{:}).WholeMatrix(:,7))) < size(Info.(RTotal{:}).WholeMatrix,1)) )
                fprintf('\nAdding entries with %s in Info_pooled_updated.%s.%s.WholeMatrix to Info_pooled_updated.',RInPlaneRealizStr,VolInPooled,RTotal{:});
                Info_pooled_updated.(VolInPooled).(RTotal{:}).WholeMatrix = cat(1,Info_pooled_updated.(VolInPooled).(RTotal{:}).WholeMatrix, Info.(RTotal{:}).WholeMatrix);
                WholeMatrixAddedFlag = true;
            else
                fprintf('\n''Info_pooled_updated.%s.%s.%s.WholeMatrix'' contains already entries with %s. Skip adding.',VolInPooled,RTotal{:},Slcstr,RInPlaneRealizStr)
            end            
            
        
        
        
        
        end
        
        
    end
    
    
end
fprintf('\n')




%% 4. Postparations

% fclose(fid)






