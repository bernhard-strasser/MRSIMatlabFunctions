function GetMatlabToolbox(toolboxname,pauselength)
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
% -         toolboxname               ...          Name of the toolbox you want to check out.
% -                                                USE THE PROPER NAMES, IF A WRONG NAME IS GIVEN, IT WILL TRY TO CHECKOUT SOMETHING FOR INFINITY!
%                                                  Examples: 'statistics_toolbox', 'image_toolbox'
% -         pauselength               ...          Length of pause in seconds between trying to check out the toolbox (default: 5s).
%
% Output:
% -         None
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations

if(~exist('toolboxname','var'))
    fprintf('\nWhich toolbox do you want to check out?\nExamples: ''statistics_toolbox'', ''image_toolbox''\n')
	return;
end
if(~exist('pauselength','var'))
	pauselength = 5;
end


% 0.3 Definitions
    






%% 1. Compute A

checked_out = 0;
pause on;

while(checked_out == 0)
    checked_out = license('checkout',toolboxname) %#ok
    pause(pauselength);
end

pause off;





%% 3. Postparations

% fclose(fid)






