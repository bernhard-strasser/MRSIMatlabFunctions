function Structy = supp_SetAllFields(Structy,CellOfFieldsToSet,ValuesToSet)
%
% supp_SetAllFields Set(add) several fields of an existing structure at once
%
% This function was written by Bernhard Strasser, August 2020.
%
%
% The function sets several fields at once of a structure Structy with certain fieldnames saved as a cell in CellofFieldsToSet with values ValuesToSet. E.g. you can set many fieldnames at once to NaN or 0 etc.
%
%
% [MRStruct.Data, MRStruct.mdhInfo] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         Structy                    ...     Structure to which the fields should be added.
% -         CellOfFieldsToSet          ...     Cell of the field names that should be added. 
% -         ValuesToSet                ...     Vector of values that should be set.
%
% Output:
% -         Structy                    ...     Output structure, same as input but with added fields
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None






%% 0. Preparations

if(numel(ValuesToSet) == 1)
    ValuesToSet = repmat(ValuesToSet,[1 numel(CellOfFieldsToSet)]);
end

%%

for ii = 1:numel(CellOfFieldsToSet)
    Structy.(CellOfFieldsToSet{ii}) = ValuesToSet(ii);
end



      

