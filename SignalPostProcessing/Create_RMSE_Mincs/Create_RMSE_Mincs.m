function Create_RMSE_Mincs(folder,folder_Reference)
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

	if(~exist('folder','var'))
		return
	end



	% 0.3 Definitions









	%% 1. Get Ref and other maps
	
	[tNAA_Ref, tCho_Ref, tCr_Ref] = GetTotalMetMaps(folder_Reference);
	[tNAA, tCho, tCr] = GetTotalMetMaps(folder);



	%% 2. Calculate RMSE-maps

	RMSE = sqrt(( (tNAA-tNAA_Ref).^2 + (tCho-tCho_Ref).^2 + (tCr-tCr_Ref).^2 )/3);




	%% 3. Write .mnc file
	
	OutFile = [folder '/maps/Clip/RMSE'];
	write_RawFiles(RMSE,[OutFile '.raw'])
	Raw2MncLike = [folder '/maps/Clip/NAA+NAAG_amp_clip_map.mnc'];
	unix(['rawtominc -float -clobber -like ' Raw2MncLike ' -input ' OutFile '.raw' ' ' OutFile '.mnc']);
	delete([OutFile '.raw']);
	
	
	%% 3. Postparations

	% fclose(fid)
end




function [tNAA, tCho, tCr] = GetTotalMetMaps(folder)

	%% 1. Load Maps
	load([folder '/maps/AllMaps.mat'])

	% Find out new or old version of AllMaps.mat
	NewVersion_flag = false;
	if(exist('AllMaps','var') && isstruct(AllMaps))
		NewVersion_flag = true;
	end

	
	%% 2. Define tNAA, tCho, tCr

	% Currently not implemented
	if(NewVersion_flag)
		fprintf('\nNew versions of AllMaps.mat not supported yet. Abort...\n')
		return
	else

		Pos_tNAA = find(~cellfun('isempty',(regexpi(MetData_amp_title,'NAA\+NAAG'))));
		tNAA = MetData_amp_clipped(:,:,:,Pos_tNAA);

		Pos_tCho = find(~cellfun('isempty',(regexpi(MetData_amp_title,'GPC\+PCh'))));
		tCho = MetData_amp_clipped(:,:,:,Pos_tCho);

		Pos_tCr = find(~cellfun('isempty',(regexpi(MetData_amp_title,'Cr\+PC'))));
		tCr = MetData_amp_clipped(:,:,:,Pos_tCr);

	end

end



