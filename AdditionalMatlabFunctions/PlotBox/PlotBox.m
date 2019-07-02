function [PlotDidWork] = PlotBox(Positions_x,Positions_y,Positions_z,Alpha,varargin)
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

if(~exist('Positions_x','var'))
    return;
end
if(~exist('Positions_y','var'))
    return;
end
if(~exist('Positions_z','var'))
    return;
end
if(~exist('Alpha','var'))
    Alpha = 1;
end
% if(numel(Alpha) == 1)
%     Alpha = repmat(Alpha,[1 numel(Positions_x)]);
% end

% % 0.3 Definitions
% NoOfStringInput = 0;
% if(exist('varargin','var') && ~ischar(varargin{1}))
%     NoOfStringInput = 1;
%     if(numel(varargin{1}) > 1)
%         Positions_z = varargin{1};
%         if(numel(varargin) > 1 && ~ischar(varargin{2}) && numel(varargin{2}) == 1)
%             NoOfStringInput = NoOfStringInput+1;
%             Radius = varargin{2};
%         end
%     else
%         Positions_z = zeros(size(Positions_x));
%         Radius = varargin{1};
%     end
% 
% end





%% 2. Plot

% The lower surface
patchline([min(Positions_x) max(Positions_x) max(Positions_x) min(Positions_x) min(Positions_x)],[min(Positions_y) min(Positions_y) max(Positions_y) max(Positions_y) min(Positions_y)],[min(Positions_z) min(Positions_z) min(Positions_z) min(Positions_z) min(Positions_z)],'edgealpha',Alpha)
% The upper surface
patchline([min(Positions_x) max(Positions_x) max(Positions_x) min(Positions_x) min(Positions_x)],[min(Positions_y) min(Positions_y) max(Positions_y) max(Positions_y) min(Positions_y)],[max(Positions_z) max(Positions_z) max(Positions_z) max(Positions_z) max(Positions_z)],'edgealpha',Alpha)
% The 4 Connections between them
patchline([min(Positions_x) min(Positions_x)],[min(Positions_y) min(Positions_y)],[min(Positions_z) max(Positions_z)],'edgealpha',Alpha)
patchline([max(Positions_x) max(Positions_x)],[min(Positions_y) min(Positions_y)],[min(Positions_z) max(Positions_z)],'edgealpha',Alpha)
patchline([min(Positions_x) min(Positions_x)],[max(Positions_y) max(Positions_y)],[min(Positions_z) max(Positions_z)],'edgealpha',Alpha)
patchline([max(Positions_x) max(Positions_x)],[max(Positions_y) max(Positions_y)],[min(Positions_z) max(Positions_z)],'edgealpha',Alpha)




%% 3. Postparations

PlotDidWork = true;






