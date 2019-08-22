function [RGBImageCell] = MakeRGBImageCell(ImageCell,CLims,SetFunction,AxisSquare,colormaptype)
%
% MakeRGBImageCell This function takes data and imagesc's them to create cell of RGB-images.
%
% This function was written by Bernhard Strasser, August 2019.
%
%
% The function [Description].
% 
%
%
% [RGBImageCell] = MakeRGBImageCell(ImageCell,CLims,SetFunction,AxisSquare)
%
% Input: 
% -         ImageCell               ...    This is the data that should be imagesc'd. Should be cell of 2D-matrices
% -         CLims                   ...    This is the clim (imagesc(data,clim)), which either should be 2-element vector, or cell of 2-element vectors
% -         SetFunction             ...    This function will be called after imagesc, and can set all properties like "imagesc(...), set(gca,...).
%                                          Example: SetFunction = @(gca) set(gca,'xtick',[],'ytick',[]);
% -         AxisSquare              ...    If 'AxisSquare', the command 'axis square' will be called after imagesc, and in that case the resulting figure needs
%                                          to have same rows as columns.
%
% Output:
% -         RGBImageCell            ...     This is the cell of RGB images
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions

% 0.1 Preparations
if(~exist('ImageCell','var'))
    fprintf('Error in MakeImageMontage: Need at least one input.');
    RGBImageCell = 0;
    return
end
if(~exist('CLims','var') || isempty(CLims))
    CLims = [-Inf Inf];
end
if(~exist('SetFunction','var') || isempty(SetFunction))
    SetFunction = @(x) FunctionTemplateBstrasser(x);   % Just to have a dummy. This function will do nothing, and therefore not change the images.
end
if(~exist('AxisSquare','var') || isempty(AxisSquare))
    AxisSquare = 'AxisNormal';
end
if(~exist('colormaptype','var') || isempty(colormaptype))
    colormaptype = 'hot';  
end


if(~iscell(CLims))
    CLims = repmat({CLims},size(ImageCell));
elseif(numel(CLims) ~= numel(ImageCell))
    CLims = myrepmat(CLims,size(ImageCell));    
end
if(~iscell(SetFunction))
    SetFunction = repmat({SetFunction},size(ImageCell));
elseif(numel(SetFunction) ~= numel(ImageCell))
    SetFunction = myrepmat(SetFunction,size(ImageCell));    
end

% 0.3 Definitions


%% 1. Open Figure

RGBImageCell = cell(size(ImageCell));
fighaendel = figure;
for xInd = 1:size(ImageCell,1)
    for yInd = 1:size(ImageCell,2)
        for ContLoop = 1:500
            imagesc(squeeze(ImageCell{xInd,yInd}),CLims{xInd,yInd}), colormap(colormaptype),  SetFunction{xInd,yInd}(gca);
            if(strcmpi(AxisSquare, 'AxisSquare'))
                axis square
            end
            CurFiggy = getframe();

            % Check if we are satisfied with Figure
            if( strcmpi(AxisSquare, 'AxisSquare') )
                if(size(CurFiggy.cdata,1) == size(CurFiggy.cdata,2) )
                    break
                else
                    close(fighaendel)
                    fighaendel = figure;
                end
            end
                        
        end
        RGBImageCell{xInd,yInd} = CurFiggy.cdata;
        
    end
end
close(fighaendel)


%% 3. Postparations





