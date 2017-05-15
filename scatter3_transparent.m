function [PlotDidWork] = scatter3_transparent(Alpha,Position_x,Position_y,varargin)
%
% scatter3_transparent Plot similar as scatter3 but with possibility for transparent points.
%
% This function was written by Bernhard Strasser, [month] [year].
%
%
% This function
% 
%
%
% [A,B] = read_csi_dat_1_10(Position_x,Position_y)
%
% Input: 
% -         Alpha                        ...    The transparency (alpha-values) of the points.
% -         Position_x                   ...    The position of the points in x-direction
% -         Position_y                   ...    The position of the points in y-direction
% -         varargin                     ...    All the other inputs. Can be:
%                                               - Position_z, if scatter3, otherwise this function replaces scatter2
%                                               - Radius of the points. Same as in scatter3.
%
% Output:
% -         PlotDidWork                 ...     Bool if plotting worked or not.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations
PlotDidWork = false;

if(~exist('Position_x','var'))
    return;
end
if(~exist('Position_y','var'))
    return;
end
if(~exist('Alpha','var'))
    Alpha = 1;
end
if(numel(Alpha) == 1)
    Alpha = repmat(Alpha,[1 numel(Position_x)]);
end

% 0.3 Definitions
NoOfStringInput = 0;
if(exist('varargin','var') && ~ischar(varargin{1}))
    NoOfStringInput = 1;
    if(numel(varargin{1}) > 1)
        Position_z = varargin{1};
        if(numel(varargin) > 1 && ~ischar(varargin{2}) && numel(varargin{2}) == 1)
            NoOfStringInput = NoOfStringInput+1;
            Radius = varargin{2};
        end
    else
        Position_z = zeros(size(Position_x));
        Radius = varargin{1};
    end

end





%% 2. Plot

% Create a shpere
[x, y, z] = sphere(50);

for i=1:size(Position_x)
    if(numel(varargin) > NoOfStringInput)
        surface(x*sqrt(Radius/pi)/30+Position_x(i),y*sqrt(Radius/pi)/30+Position_y(i),z/2*sqrt(Radius/pi)/30+Position_z(i),varargin{NoOfStringInput+1:end},'FaceAlpha',Alpha(i));
        %pb=patch((sin(Phi)*sqrt(Radius/pi)/30+ Position_x(i)),(cos(Phi)*sqrt(Radius/pi)/30+Position_y(i)), (zeros(size(Phi))+Position_z(i)),varargin{NoOfStringInput+1:end},'FaceAlpha',Alpha(i));
    else
        surface(x*sqrt(Radius/pi)/30+Position_x(i),y*sqrt(Radius/pi)/30+Position_y(i),z/2*sqrt(Radius/pi)/30+Position_z(i),'FaceAlpha',Alpha(i));
        %pb=patch((sin(Phi)*sqrt(Radius/pi)/30+ Position_x(i)),(cos(Phi)*sqrt(Radius/pi)/30+Position_y(i)), (zeros(size(Phi))+Position_z(i)),'FaceAlpha',Alpha(i));
    end
    %alphamap(pb,Alpha)
    %alpha(pb,Alpha);
end



daspect([1 1 1])


%% 3. Postparations

% fclose(fid)






