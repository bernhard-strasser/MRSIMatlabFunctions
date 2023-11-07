function FieldMBytes = CheckSizesOfFields(Struct)
%
% CheckSizesOfFields Calculate sizes of all fields of a structure
%
% This function was written by Bernhard Strasser, August 2022.
%
%
% The function takes a structure as input and calculates the size in MBytes of all fields.
%
%
% Input: 
% -         Struct                    ...               Structure for which the sizes of its fields should be calculated.
%
% Output:
% -         FieldMBytes                         ...     Structure with same fields as input structure, but with sizes in MBytes
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None



for fieldy = transpose(fieldnames(Struct))
    Tmp = Struct.(fieldy{1}); %#ok 
    Tmp2 = whos('Tmp'); clear Tmp;
    FieldMBytes.(fieldy{1}) = Tmp2.bytes / 2^20;
    
end


