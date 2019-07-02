function Data = read_JmruiTxtData(file)
%
% read_csi_dat Read in raw data from Siemens
%
% This function was written by Bernhard Strasser, Dec 2014.
%
%
% The function reads in the header info of a jmrui txt file.
%
%
% Data = read_JmruiTxtData(file)
%
% Input: 
% -         file                    ...     jmrui txt file.
%
% Output:
% -         Data                    ...     The read in header info.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: 





%% 0. Preparation

if(~exist('file','var'))
	fprintf('\nError: No file passed over.\n')
	Data = 0;
	return
end



%% 1. Gather information from header







data_raw = importdata(sprintf('%s', file), '\t', 21);       %imports data from the in_file and gives an cell array with one struct for each metabolite
Data = data_raw.data;




