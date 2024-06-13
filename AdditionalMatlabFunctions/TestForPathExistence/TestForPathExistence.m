function [EverythingExists] = TestForPathExistence(InputStruct,quiet_flag,TestThoseLocally,ServerName,TestThoseOnServer)
%
% TestForPathExistence Tests if the entries of a struct exist as dirs or files.
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
% -         inputvar1                   ...    This is the first input
% -         inputvar2                   ...    And this the second
%
% Output:
% -         A                           ...     This is output A
% -         B                           ...     This is output B
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations
EverythingExists = true;
if(~exist('InputStruct','var'))
	fprintf('\nError: No InputStruct given to test for.\n')
	return;
end

if(~isstruct(InputStruct))
	fprintf('\nError: Variable InputStruct must be structure.\n')
	return
end
if(~exist('quiet_flag','var'))
	quiet_flag = false;
end

Fields = fieldnames(InputStruct);
if(~exist('TestThoseLocally','var') || isempty(TestThoseLocally))
	TestThoseLocally = Fields;
end
if(~exist('TestThoseOnServer','var') || isempty(TestThoseLocally))
	TestThoseOnServer = Fields;
end


% 0.3 Definitions
    





%% 1. Loop over all fields and subfields of those etc.

for CurField = transpose(Fields)
	CurFieldStr = CurField{:};
	
	
	IsCell = iscell(InputStruct.(CurFieldStr));
	if(IsCell)
		LoopMax = numel(InputStruct.(CurFieldStr));
	else
		LoopMax = 1;
	end

	
	for CurCellNo = 1:LoopMax

		if(IsCell)
			PathToTest = InputStruct.(CurFieldStr){CurCellNo};
		else
			PathToTest = InputStruct.(CurFieldStr);
		end


		TestLocal = sum(~cellfun('isempty',regexp(CurFieldStr,TestThoseLocally))) > 0;
		if(exist('ServerName','var'))
			TestOnServer = sum(~cellfun('isempty',regexp(CurFieldStr,TestThoseOnServer))) > 0;
		else
			TestOnServer = false;
        end
        
        [bla,hostname] = unix('hostname');
        if(strcmp(hostname(1:end-1),ServerName) || ~isempty(regexp(hostname(1:end-1),[ServerName '\.+'],'ONCE')) || ~isempty(regexp(ServerName,[hostname(1:end-1) '\.+'],'ONCE')))
            TestOnServer = false; TestLocal = true;
        end
        

		% Test On Server
		if(TestOnServer)
			[bla,BashReturn] = unix(['ssh ' ServerName ' "if ! [[ -e ' PathToTest ' ]]; then echo wuffwi; fi; exit"']);
			if(~isempty(strfind(BashReturn,'wuffwi')))
				EverythingExists = false;
				if(~quiet_flag)
					fprintf('\nError: The path of %s,\n%s\ndoes not exist on server "%s". Abort...\n',CurFieldStr,PathToTest,ServerName)
				end
				return
			end
		end
		
		% Test Local
		if(TestLocal)
			if(~exist(PathToTest,'file') && ~exist(PathToTest,'dir') && ~strcmp(PathToTest,'lcmodel'))
				EverythingExists = false;
				if(~quiet_flag)
					fprintf('\nError: The path of %s,\n%s\ndoes not exist. Abort...\n',CurFieldStr,PathToTest)
				end
				return
			end
		end


	end
	
	
end
	





%% 2. Compute B





%% 3. Postparations

% fclose(fid)






