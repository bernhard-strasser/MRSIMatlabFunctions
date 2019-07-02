function ErrorOcurred = ConvertOldToNewMetMaps()
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
% -         inputvar1                   ...    This is the first input
% -         inputvar2                   ...    And this the second
%
% Output:
% -         A                           ...     AllMaps from the new
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations




%% 1. Define string what to do in base workspace


% We need to perform this action in the base workspace. Therefore only write strings here, and perform evalin at the end.
evalin_str = [];

% Metabos
evalin_str = [evalin_str 'AllMaps.Metabos.Normal = MetData_amp;'];
evalin_str = [evalin_str 'AllMaps.Metabos.Clip = MetData_amp_clipped;'];
evalin_str = [evalin_str 'AllMaps.Metabos.Res = MetData_amp_res;'];
evalin_str = [evalin_str 'AllMaps.Metabos.Res_Clip = MetData_amp_res_clipped;'];
evalin_str = [evalin_str 'AllMaps.Metabos.Title = MetData_amp_title;'];


% Extras
evalin_str = [evalin_str 'AllMaps.Extras.Normal = MetData_extra;'];
evalin_str = [evalin_str 'AllMaps.Extras.Res = MetData_extra_res;'];
evalin_str = [evalin_str 'AllMaps.Extras.Title = MetData_extra_title;'];


% CRLBs
evalin_str = [evalin_str 'AllMaps.CRLBs.Normal = MetData_sd;'];
evalin_str = [evalin_str 'AllMaps.CRLBs.Clip = MetData_sd_clipped;'];
evalin_str = [evalin_str 'AllMaps.CRLBs.Res = zeros(size(MetData_sd_res_clipped));'];			% Does not exist in old
evalin_str = [evalin_str 'AllMaps.CRLBs.Res_Clip = MetData_sd_res_clipped;'];
evalin_str = [evalin_str 'AllMaps.CRLBs.Title = MetData_sd_title;'];


% Ratios
evalin_str = [evalin_str 'AllMaps.Ratios.Normal.RatioToNAANAAG = MetData_Ratio;'];
evalin_str = [evalin_str 'AllMaps.Ratios.Res.RatioToNAANAAG = MetData_res_Ratio;'];
evalin_str = [evalin_str 'AllMaps.Ratios.Normal.RatioToNAANAAG_Title = MetData_Ratio_title;'];

evalin_str = [evalin_str 'AllMaps.Ratios.Normal.RatioToCrPCr = zeros(size(MetData_Ratio));'];	% Does not exist in old
evalin_str = [evalin_str 'AllMaps.Ratios.Normal.RatioToCrPCr_Title = MetData_Ratio_title;'];	% Does not exist in old
evalin_str = [evalin_str 'AllMaps.Ratios.Res.RatioToCrPCr = zeros(size(MetData_res_Ratio));'];	% Does not exist in old


% Masks
evalin_str = [evalin_str 'AllMaps.Masks.Normal = mask;'];
evalin_str = [evalin_str 'AllMaps.Masks.Res = mask_res;'];



% clear variables
evalin_str = [evalin_str 'clear MetData_amp MetData_amp_clipped MetData_amp_res MetData_amp_res_clipped MetData_amp_title MetData_extra MetData_extra_res MetData_extra_title'];
evalin_str = [evalin_str ' MetData_sd MetData_sd_clipped MetData_sd_res_clipped MetData_sd_title MetData_Ratio MetData_res_Ratio MetData_Ratio_title mask mask_res'];


% 0.3 Definitions
    






%% 2. Eval evalin_str in base workspace

ErrorOcurred = 0;
try
	evalin('base',evalin_str);
catch
	ErrorOcurred = 1;
end



%% 3. Postparations

% fclose(fid)






