function ShowSequenceArbGrads(PathToTxtFile,inputvar2)
%
% ShowSequenceArbGrads Do nothing specific
%
% This function was written by Bernhard Strasser, Jannuary 2015.
%
%
% The function can really do nothing, and more specifically, exactly nothing.
% 
%
%
% [A,B] = read_csi_dat_1_10(inputvar1,inputvar2)
%
% Input: 
% -         PathToTxtFile                   ...    Path to text file containing the arbitrary gradient values
%
% Output:
% -         A                           ...     This is output A
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations

if(~exist('PathToTxtFile','var'))
    fprintf('\nNo input file. Abort.')
    return
end



% 0.3 Definitions
    





%% 1. Read File

GradValues = importdata(PathToTxtFile,'\t');

GradValues.data(isnan(GradValues.data)) = [];






%% 2. Compute B





%% 3. Postparations

% fclose(fid)






