function figure_handle_out = Plot2DPatterns(Pattern,figure_handle_in,repmat_flag)
%
% deg2rad Convert from degree to radians.
%
% This function was written by Bernhard Strasser, July 2012.
%
%
%
% Input: 
% -         deg					            ...     An angle in degree.
%
% Output:
% -        rad                              ...     An angle in radians.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None





%% 0. Preparations

% Assign standard values to variables if nothing is passed to function.
if(nargin < 2 || figure_handle_in == 0)
	figure_handle_out = figure;
else
	figure_handle_out = figure_handle_in;
end

Pattern = rot90(Pattern,-1);

if(repmat_flag)
	Pattern = repmat(Pattern,[3 3]);
	Pattern_Center = zeros(size(Pattern));
	Pattern_Center( size(Pattern,1)/3+1 : 2*size(Pattern,1)/3, size(Pattern,2)/3+1 : 2*size(Pattern,2)/3 ) = ...
	Pattern(size(Pattern,1)/3+1 : 2*size(Pattern,1)/3, size(Pattern,2)/3+1 : 2*size(Pattern,2)/3);
	NotPattern_Center = zeros(size(Pattern));
	NotPattern_Center( size(Pattern,1)/3+1 : 2*size(Pattern,1)/3, size(Pattern,2)/3+1 : 2*size(Pattern,2)/3 ) = ...
	~Pattern(size(Pattern,1)/3+1 : 2*size(Pattern,1)/3, size(Pattern,2)/3+1 : 2*size(Pattern,2)/3);	
end


[Ax, Ay, Az] = ind2sub(size(Pattern), find(Pattern));
[notAx, notAy, notAz] = ind2sub(size(~Pattern), find(~Pattern));

if(repmat_flag)
	[Axctr, Ayctr, Azctr] = ind2sub(size(Pattern_Center), find(Pattern_Center));
	[notAxctr, notAyctr, notAzctr] = ind2sub(size(NotPattern_Center), find(NotPattern_Center));
end

% Pattern itself
scatter(Ax, Ay, 150,'r','fill','LineWidth',1.5);
hold on
%scatter(notAx, notAy, 150,'r','LineWidth',1.5);

if(repmat_flag)
	
	scatter(Axctr,Ayctr,150,'b','fill','LineWidth',1.5)
	%scatter(notAxctr,notAyctr,150,'b','LineWidth',1.5)

end

hold off;

% Semi-thick Grid Lines
grid on
axh = gca;
set(axh,'GridLineStyle','--')
set(axh,'LineWidth',1.5)

% ticks
%set(axh,'XTickMode','manual','FontSize',14,'ticklength',[0.03 0.007]);
set(axh,'YTickMode','manual');
set(axh,'XTick',1:size(Pattern,1),'FontSize',20)
set(axh,'YTick',1:size(Pattern,2),'FontSize',20)
set(axh,'YTickLabel',size(Pattern,2):-1:1)
set(axh,'xaxisLocation','top')


axis([0.5 size(Pattern,1)+0.5 0.5 size(Pattern,2)+0.5])
box on
axis square
xlabel('x','FontSize',20)
ylabel('y','FontSize',20)



