function [memused_inGB,memfree_inGB,memtot_inGB] = memused_inGB_linux2(quiet_flag)
%
% memused_inGB_linux Show memory usage of MATLAB and free memory in GiB.
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function shows the memory that is used by MATLAB in percent of all memory and the free system memory. Should work on all linux systems, probably even on all UNIX
% 
%
%
% memused_inGB = memused_inGB_linux(quiet_flag)
%
% Input: 
% -         quiet_flag                  ...    If 1, nothing is printed to display.
%
% Output:
% -         memused_inGB                     ...    used memory by MATLAB, in GB
% -         memfree_inGB                     ...    free memory on system, in GB
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

% Check if we are on a linux system
[UnixSystem,memfree_inGB] = unix('free -g'); UnixSystem = ~UnixSystem;
if(~UnixSystem)
	if(~quiet_flag)
		fprintf('\nWarning: memused_inGB_linux cannot run on non-linux systems')
	end
	memused_inGB=0;memfree_inGB=0;
	return;
end

% 0.3 Definitions
    




%% 1. Find out used memory


% Find out pid
pid = feature('getpid');


% Find out username.
[stat,uname] = unix('whoami');
[stat,uid] = unix('id');
uid = regexp(regexp(uid,'uid=\d+','match'),'\d+','match');
uid = uid{:}; uid = str2num(uid{:});

[stat,memused_inGB] = unix(['ps aux | grep -i "' num2str(pid) '" | grep "MATLAB" | awk ''{print $4"\t"$11}'' | cut -f 1']);     % unix just performs the unix-command.
memused_inGB = sum(str2num(memused_inGB));                                                               % str2double does not work

% free changed from version 3.3.9 to 3.3.10 outputting different fields
[stat,freeVersion] = unix('free -V');
freeVersion = freeVersion(regexp(freeVersion,'[\d+\.]+\d+'):end);
freeVersion = strsplit(freeVersion,'.');
freeVersion = str2double(regexprep(freeVersion,'\D',''));

[stat,memtot_inGB] = unix('free -m | awk ''NR==2'' | awk ''{print $2}'''); 
memtot_inGB = str2num(memtot_inGB)/2^10;
if(any(isnan(freeVersion)) || freeVersion(1) > 3 || freeVersion(2) > 3 || freeVersion(3) > 9)
    [stat,memfree_inGB] = unix('free -m | awk ''NR==2'' | awk ''{print $7}''');     
else
    [stat,memfree_inGB] = unix('free -m | awk ''NR==3'' | awk ''{print $4}''');  % unix just performs the unix-command. perform free -g command,
end                                                                         % take 3rd line of that, and 4th column, i.e. "-/+ buffers/cache" of free column
clear stat
memfree_inGB = sum(str2num(memfree_inGB))/2^10;                                                               % str2double does not work
memused_inGB = memused_inGB * memtot_inGB / 100;

%% 2. Display used memory

if(~quiet_flag)
    Datie = datetime('now','TimeZone','local','Format','yyyy-MM-dd HH:mm:ss.ms');
    fprintf('\nTimeStamp %s: Currently %f GB of %f GB (%f %%) are used by MATLAB. ',string(Datie),memused_inGB,memtot_inGB,memused_inGB/memtot_inGB*100)
    fprintf('%f GB are still free.\n',memfree_inGB)
end




%% 3. Postparations

% fclose(fid)






