function Output = read_JmruiTxtFile(file)
%
% read_csi_dat Read in raw data from Siemens
%
% This function was written by Bernhard Strasser, Dec 2014.
%
%
% The function reads in the header info of a jmrui txt file.
%
%
% Output = read_JmruiTxtFile(file)
%
% Input: 
% -         file                    ...     jmrui txt file.
%
% Output:
% -         Output                  ...     consisting of fields .Data and .Header
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: 





%% 0. Preparations


% Find out memory used by MATLAB
memused_before = memused_linux(1); 








%% 1. Gather information from header



Output.Header = read_JmruiTxtHeader(file);





%% 2. READ DATA

tic
fprintf('\nReading data\t\t\t...')
        
Output.Data = read_JmruiTxtData(file);

fprintf('\n\t\t\t\t...took\t%10.6f seconds',toc)       




%% 7. Postparations

memused_after = memused_linux(1); 
display([char(10) 'The function used ' num2str(memused_after-memused_before) '% of the total memory.'])


