function x = fixnan(x)
%will work no matter how many levels, including if cells are not consistent types
    if isnumeric(x)
        x(isnan(x)) = 0;
    elseif iscell(x)
        x = cellfun(@fixnan, x, 'uniform', 0);
    %else other types do nothing, return unchanged
    end
end