function [RGBImageCell] = MakeRGBImageCell(ImageCell,CLims,SetFunction,AxisAspect,colormaptype,DemandImageSize)
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
% [RGBImageCell] = MakeRGBImageCell(ImageCell,CLims,SetFunction,AxisAspect)
%
% Input: 
% -         ImageCell               ...    This is the data that should be imagesc'd. Should be cell of 2D-matrices
% -         CLims                   ...    This is the clim (imagesc(data,clim)), which either should be 2-element vector, or cell of 2-element vectors
% -         SetFunction             ...    This function will be called after imagesc, and can set all properties like "imagesc(...), set(gca,...).
%                                          Example: SetFunction = @(gca) set(gca,'xtick',[],'ytick',[]);
% -         AxisAspect              ...    If 'AxisSquare', the command 'axis square' will be called after imagesc, and in that case the resulting figure needs
%                                          to have same rows as columns. If set to [] or not input, nothing will be called. If a [1 3]-vector, this vector will
%                                          be called with daspect(AxisAspect).
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
if(~exist('AxisAspect','var') || isempty(AxisAspect))
    AxisAspect = 'AxisNormal';
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
DemImSiExists = exist('DemandImageSize','var') && ~isempty(DemandImageSize);
for xInd = 1:size(ImageCell,1)
    for yInd = 1:size(ImageCell,2)
        for ContLoop = 1:100
            imagesc(squeeze(ImageCell{xInd,yInd}),CLims{xInd,yInd}), colormap(colormaptype),  SetFunction{xInd,yInd}(gca);
            if(ischar(AxisAspect) && strcmpi(AxisAspect, 'AxisSquare'))
                axis square
            elseif(isnumeric(AxisAspect) && numel(AxisAspect) == 3)
                daspect(AxisAspect)
            end
            CurFiggy = getframe();

            % Check if we are satisfied with Figure sizes
            if( strcmpi(AxisAspect, 'AxisSquare') )
                ImageSizeRatio = 1;       % For AxisSquare, the size of the images must be square, i.e. same amount of x- and y-pixels
            else                          % Otherwise, we can caluclate the number of pixels by the sizes of the image in x and y, and the AxisAspect
                ImageSizeRatio = size(ImageCell{xInd,yInd},2)*AxisAspect(2)/(size(ImageCell{xInd,yInd},1)*AxisAspect(1));
            end
            BreakCondition = DemImSiExists && all(size_MultiDims(CurFiggy.cdata,[1 2]) == DemandImageSize);
            BreakCondition = BreakCondition || (~DemImSiExists && round(size(CurFiggy.cdata,1)*ImageSizeRatio) == size(CurFiggy.cdata,2));
            BreakCondition = BreakCondition || (xInd > 1 || yInd > 1);
            if(BreakCondition )
                break
            else
                close(fighaendel)
                fighaendel = figure;
            end
                                      
        end
        RGBImageCell{xInd,yInd} = CurFiggy.cdata;
        
    end
end
close(fighaendel)


%% 3. Postparations





