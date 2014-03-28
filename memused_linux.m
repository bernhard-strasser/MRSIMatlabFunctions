function memused = memused_linux_1_0(quiet_flag)
%
% memused_linux_1_0 Show memory usage of MATLAB
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function shows the memory that is used by MATLAB in percent of all memory. Should work on all linux systems, probably even on all UNIX
% 
%
%
% memused = memused_linux_1_0(quiet_flag)
%
% Input: 
% -         quiet_flag                  ...    If 1, nothing is printed to display.
%
% Output:
% -         memused                     ...     used memory by MATLAB, in percent of all memory
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations

if(~exist('quiet_flag','var'))
    quiet_flag = 0;
end


% 0.3 Definitions
    






%% 1. Find out used memory

[stat,memused] = unix('ps aux | awk ''{print $4"\t"$11}'' | grep -i "matlab" | cut -f 1');          % unix just performs the unix-command.
clear stat
memused = sum(str2num(memused));                                                                    % str2double does not work




%% 2. Display used memory

if(~quiet_flag)
    display([ char(10) num2str(memused) '% of total memory is used.' char(10)])
end




%% 3. Postparations

% fclose(fid)






