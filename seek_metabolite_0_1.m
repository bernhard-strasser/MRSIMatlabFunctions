function [metabolite_found,metpeak_specpoint] = seek_metabolite_0_1(spectrum_seekregion,seek_criteria,seek_criteria_values)
% Function that seeks for a metabolite within a given region with some criteria. 
% The criteria are: - 'metpeak_height': checks if the assumed metpeak (max of region) is metpeak_height_criterium_value times higher than the border of the region; this value must be given to the function
% 					- 'derivative': checks if the derivative of the spectrum is derivative_criterium_value times higher than at the border; This value and a averaging value must be given (compute the derivative between neighbouring points, or farther points) NOT IMPLEMENTED YET


metpeak_specpoint = find(max(spectrum_seekregion) == spectrum_seekregion);
check_metpeak_height_crit = strcmpi(seek_criteria,'metpeak_height');
check_metpeak_height_crit_flag = logical(sum(check_metpeak_height_crit));
check_metpeak_height_crit_value = seek_criteria_values(check_metpeak_height_crit);

if(check_metpeak_height_crit_flag)
	metpeak_higher_than_leftborder = logical(max(spectrum_seekregion) > check_metpeak_height_crit_value*spectrum_seekregion(1));			% max for left because high chemshift --> low point			
	metpeak_higher_than_rightborder = logical(max(spectrum_seekregion) > check_metpeak_height_crit_value*spectrum_seekregion(end));
	metpeak_height_crit = metpeak_higher_than_leftborder && metpeak_higher_than_rightborder;
else
	metpeak_height_crit = true; 			% If it gets not checked it should not have any influence for metabolite_found
end


% More criteria needed ???



metabolite_found = metpeak_height_crit;
