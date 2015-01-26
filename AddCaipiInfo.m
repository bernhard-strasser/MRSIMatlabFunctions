function Info_pooled_updated = AddCaipiInfo(InfoPath,Info_pooled)
%
% template_1_0 Do nothing specific
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
% -         InfoPath                               ...    This is the first input
% -         Info_pooled                            ...    And this the second
%
% Output:
% -         Info_pooled_updated                    ...     This is output A
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
    return
end
if(~exist(InfoPath,'file'))
    fprintf('\nFile %s does not exist.',InfoPath)
    return
end



% 0.3 Definitions
    






%% 1. Load Info

try
    load(InfoPath)
catch err
    fprintf('\nCould not load %s due to the following error:\n%s',InfoPath,err.message)
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
    WholeMatrixAddedFlag = false;
    for Slcflag = 1:2
        Slcstr = ['RSlc_' num2str(Slcflag)];
        
        % If the RSlc doesnt exist in the Info, we dont need to add it
        if(~isfield(Info.(RTotal{:}),Slcstr))
            continue
        end
        
        % Check for existence in Info_pooled of the Data which should be added
        if(isfield(Info_pooled_updated.(VolInPooled).(RTotal{:}),Slcstr) && numel(Info_pooled_updated.(VolInPooled).(RTotal{:}).(Slcstr)) > 0) % If exists and has elements, continue
            fprintf('\n''Info_pooled_updated'' contains already field %s.%s.%s. Skip adding.',VolInPooled,RTotal{:},Slcstr)
        else
            fprintf('\nAdding field %s.%s.%s to Info_pooled_updated.',VolInPooled,RTotal{:},Slcstr);
            Info_pooled_updated.(VolInPooled).(RTotal{:}).(Slcstr) = Info.(RTotal{:}).(Slcstr);
        end
        
        if(sum(Info_pooled_updated.(VolInPooled).(RTotal{:}).WholeMatrix(:,4) == Slcflag) == 0 && ~WholeMatrixAddedFlag)
            fprintf('\nAdding entries with RSlc = %f in Info_pooled_updated.%s.%s.WholeMatrix to Info_pooled_updated.',Slcflag,VolInPooled,RTotal{:});
            Info_pooled_updated.(VolInPooled).(RTotal{:}).WholeMatrix = cat(1,Info_pooled_updated.(VolInPooled).(RTotal{:}).WholeMatrix, Info.(RTotal{:}).WholeMatrix);
            WholeMatrixAddedFlag = true;            
        else
            fprintf('\n''Info_pooled_updated.%s.%s.WholeMatrix'' contains already entries with RSlc = %f. Skip adding.',VolInPooled,RTotal{:},Slcflag)
        end
        
        
    end
    
    
end
fprintf('\n')




%% 4. Postparations

% fclose(fid)






