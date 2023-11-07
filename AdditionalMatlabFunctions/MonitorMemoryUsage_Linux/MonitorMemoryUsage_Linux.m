function [UniqueMemID,memused, memtot] = MonitorMemoryUsage_Linux(UniqueMemID,StartOrStop,Logpath)




if(nargout < 1)
	fprintf('\n\nMonitorMemoryUsage_Linux: You need to get the UniqueMemID in order to kill the correct process again! Abort!')
	return;
end

memused = 0;
memtot = 0;

if(UniqueMemID <= 0)
	UniqueMemID = randi(10^6,1);
end
if(~exist('Logpath','var'))
	[err,Logpath] = unix('echo $HOME');
end
Logfile = strcat(Logpath, '/TempMemUsage', num2str(UniqueMemID), '.m');
killfile = strcat(Logpath,'/KillMemProcess', num2str(UniqueMemID), '.sh');



if(strcmpi(StartOrStop,'Start'))
	Start_flag = true;
else
	Start_flag = false;
end

if(Start_flag)
	
	if(exist(Logfile,'file'))
		delete(Logfile)
		fprintf('DOUBLE')
	end
	
	BashPath = which('MonitorMemoryUsage_Linux');
	BashPath = regexp(BashPath,'.*/','match'); BashPath = BashPath{1}; BashPath = BashPath(1:end-1);
	
	unix(['bash ' BashPath '/LogMemoryUsage -l ' Logfile ' &']);
	pause(0.5)
	run(Logfile)
	pause(0.5)
	killfid = fopen(killfile,'w+');
	fprintf(killfid,'kill -10 %s',num2str(PID));
	fclose(killfid);
	pause(1)
	
else
	
	CurMem = 0;
	[stat,memtot] = unix('free -g | awk ''NR==2'' | awk ''{print $2}'''); memtot = str2double(memtot);

	% Matlab has a funny behaviour with file caching: If the file is changed, MATLAB doesnt see these changes if it had read in this file already.
	% Therefore, create a dummy file with unique name, and run this file.
	Logfile2 = regexprep(Logfile,'\.m',''); Logfile2 = strcat(Logfile2,num2str(randi(10^6,1)), '.m');
	copyfile(Logfile,Logfile2)
	delete(Logfile)
	
	% Run logfile, kill the bash process, calculated used memory
	run(Logfile2)
	KillingFailed = unix(['kill -10 ' num2str(PID)]);
	memused = max(CurMem) * memtot / 100;
	
	% Delete dummy file
	delete(Logfile2)
	if(KillingFailed == 0)
		delete(killfile)
	end

end




