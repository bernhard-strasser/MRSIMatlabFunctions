function InArray = ReplaceValuesInNumArray(InArray, ReplaceValues, ReplaceValuesWith,maxError)
    %
    % ReplaceValuesInNumArray Replace all values in an numerical array with another value
    %
    % This function was written by Bernhard Strasser, April 2023.
    %
    %
    % The function replaces all occurrences of 'ReplaceValues' with 'ReplaceValuesWith' if the occurrences in InArray are close enough to ReplaceValues 
    %
    %
    % InArray = ReplaceValuesInNumArray(InArray, ReplaceValues, ReplaceValuesWith,maxError)
    %
    % Input: 
    % -         InArray                 Input numerical array
    % -         ReplaceValues           Replace those values in InArray
    % -         ReplaceValuesWith       Replace found occurrences with those values
    % -         maxError:               Treat values in InArray and ReplaceValues as equivalent, if they are below maxError (   abs(InArray - ReplaceValues) < maxError   )
    %
    % Output:
    % -         InArray:                 Updated output array
    %
    %
    % Feel free to change/reuse/copy the function. 
    % If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
    % Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
    % File dependancy: None

    % Further remarks:



    %% 0. Preparations

    if(numel(ReplaceValues) ~= numel(ReplaceValuesWith))
        fprintf('\nError in ReplaceValuesInNumArray: Vector ReplaceValues and ReplaceValuesWith must have same length. Do nothing on InArray.\n')
        return;
    end
    
    if(~exist('maxError','var'))
        if(strcmp(class(InArray), 'double') && strcmp(class(ReplaceValues), 'double'))
            maxError = 1E-14;
        else
            maxError = 1E-7;
        end
    end
    
    %% 1.
    
    for ii = 1:numel(ReplaceValues)
        InArray(abs(InArray - ReplaceValues(ii)) < maxError) = ReplaceValuesWith(ii);
    end
    
end
