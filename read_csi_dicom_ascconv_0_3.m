function ParList = read_csi_dicom_ascconv_0_3(csi_dicom_file)
%% 0. Preparations

% Search for these entries in the ascconv header part
ParList_Search = {'VectorSize', 'FinalMatrixSizePhase', 'FinalMatrixSizeRead', 'FinalMatrixSizeSlice', 'MaximumNofRxReceiverChannels'};
% Name the structure entries of ParList like this
ParList_Assign = {'VectorSize', 'ROW',                  'COL',                 'SLC',                  'total_channel_no'};

% Initialize ParList
for Par_no = 1:numel(ParList_Search)
    eval([ 'ParList.' ParList_Assign{Par_no} ' = 1;' ]);
end


csi_fid = fopen(csi_dicom_file,'r');

%% 1. Find & SAVE ASCCONV

begin_found = 0;
ascconv = [];
sLine = 0;

while(sLine > -1)
    
    sLine = fgets(csi_fid);
    
    if(not(begin_found))                                % If begin of ascconv not yet found
        
        
        if(strcmp(['### ASCCONV BEGIN ###' char(10)],sLine))
            begin_found = true;
        else
            continue                                % If begin not yet found and current line is also not begin --> read next line
        end
        
        
    else                                            % If begin of ascconv was found
        
        if(strcmp(['### ASCCONV END ###" ' char(10)] ,sLine))     % If the end was found --> stop while loop
            break
        else
            ascconv = [ascconv; {sLine}];           % If the begin was already found and the current line is not the end --> read in line and save it
        end
        
    end
   
        
end



%% 2. Display error & stop if no Ascconv found

if(not(begin_found))
    display('You gave me a file without any ascconv, I cannot digest that! I will stop here . . .')
    return
end


%% 3. Convert ascconv

% Convert cell array of strings containing n x 2 entries. The first entries containing the parts before the '=' (pre-=) 
% and the second parts containing the parts after the '=' (post-=)


% ascconv is here a cell array of strings (with lets say 348 entries)

% This regexp makes ascconv to a cell array with 348 entries, each of these on its own a cell array of 2 strings
ascconv = regexp(ascconv, '=','split');

% ascconv is then a 2*348 = 696x1 cell array of strings; All odd cells contain the parts before the '=', all even cells the part after the '='
ascconv = transpose([ascconv{:}]);

% Now seperate the pre-= and the post-= parts, remove all white spaces before and after the entries.
ascconv = strtrim([ascconv(1:2:end) ascconv(2:2:end)]);

% Now we are happy and have our 348x2 cell array of strings.



%% 4. Search certain entries & Save these

% Find Parameters within the pre-= entries. Result: cell array, containing [] if in the corresponding cell the Parametername was not found, and [n]
% If it was found in the corresponding cell on place n of the string; Convert to that index where the Parametername was found in the cell array


for Par_no = 1:numel(ParList_Search)
    Par_Logic = strfind(ascconv(:,1),ParList_Search{Par_no});    
    Par_Logic = not(cellfun('isempty',Par_Logic));
    eval([ 'ParList.' ParList_Assign{Par_no} ' = str2double(ascconv(Par_Logic,2));' ]);
end

% There is a bug with Slice number: It says 8, even if only 1 slice was measured. Correct that here.
ParList.SLC = 1;




%% 5. Postparations

fclose(csi_fid);
