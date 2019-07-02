function [metbasis_upfield_found,metbasis_lowfield_found,metbasis_specpoint_upfield,metbasis_specpoint_lowfield] = ...
seek_basis_of_metabolite_0_1(spectrum_seekregion,chemshift_vector,seek_criteria,seek_criteria_values,metpeak_specpoint)
% Function that seeks for a metabolite within a given region with some criteria. 
% The criteria are: - 'metpeak_height': checks if the assumed metpeak (max of region) is metpeak_height_criterium_value times higher than the border of the region; this value must be given to the function
% 					- 'derivative': checks if the derivative of the spectrum is derivative_criterium_value times higher than at the border; This value and a averaging value must be given (compute the derivative between neighbouring points, or farther points) NOT IMPLEMENTED YET

    



check_fixed_distance_crit = strcmpi(seek_criteria,'fixed_distance_to_peak');
check_fixed_distance_crit_flag = logical(sum(check_metpeak_height_crit));


if(check_fixed_distance_crit_flag)
    
    fixed_distance_crit = true;     % There is no possibility that this criterium doesn't fail.
    metbasis_chemshift_offset_upfield = seek_criteria_values(check_fixed_distance_crit);
    metbasis_chemshift_offset_lowfield = seek_criteria_values(check_fixed_distance_crit + 1);  
    
    % translate chemshift to specpoint
    
    metbasis_specpoint_offset_upfield = find(min(abs(chemshift_vector - metbasis_chemshift_offset_upfield)) == abs(chemshift_vector - metbasis_chemshift_offset_upfield));
    metbasis_specpoint_offset_lowfield = find(min(abs(chemshift_vector - metbasis_chemshift_offset_lowfield)) == abs(chemshift_vector - metbasis_chemshift_offset_lowfield));
    metbasis_specpoint_upfield = metpeak_specpoint - metbasis_specpoint_offset_upfield;
    metbasis_specpoint_lowfield = metpeak_specpoint + metbasis_specpoint_offset_lowfield;
    
    if(metbasis_specpoint_upfield < 1)
        metbasis_specpoint_upfield = 1;
    elseif(metbasis_specpoint_upfield > size(chemshift_vector,2))
        metbasis_specpoint_upfield = size(chemshift_vector,2);
    end
    if(metbasis_specpoint_lowfield < 1)
        metbasis_specpoint_lowfield = 1;
    elseif(metbasis_specpoint_lowfield > size(chemshift_vector,2))
        metbasis_specpoint_lowfield = size(chemshift_vector,2);
    end   
    
    
end






% More criteria needed ???



metbasis_upfield_found = fixed_distance_crit;
metbasis_lowfield_found = fixed_distance_crit;
