function handle = fig(Name,Value)


if(nargin > 1)
    handle = figure(Name,Value); 
else
    handle = figure;
end

if(nargout < 1)
    clear handle
end

addToolbarExplorationButtons(gcf);


