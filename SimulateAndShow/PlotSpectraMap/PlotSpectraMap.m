function FigureHandle1 = PlotSpectraMap(csi, ppm_axis,PlotFromTo_Points)
%
%




%%



ppm_split = max(max(max(max(ppm_axis(PlotFromTo_Points(1):PlotFromTo_Points(2)))))) - min(min(min(min(ppm_axis(PlotFromTo_Points(1):PlotFromTo_Points(2))))));
spec_split = max(max(max(max(csi(:,:,PlotFromTo_Points(1):PlotFromTo_Points(2)))))) -min(min(min(min(csi(:,:,PlotFromTo_Points(1):PlotFromTo_Points(2))))));
FigureHandle1 = figure(1);


%title('');          %write here which SLC is displayed
% SPECTRUM - black
for y=y_min:y_max
    for x=x_min:x_max
        figure(FigureHandle1)
        dat = csi(x,y,PlotFromTo_Points(1):PlotFromTo_Points(2));
        plot( -squeeze(ppm_axis(x,y,PlotFromTo_Points(1):PlotFromTo_Points(2)))+ppm_split*x*1.2,squeeze(dat)-(y-1)*spec_split,'k', 'LineWidth', 0.3 ); hold on; axis off

    end
end
hold off

close (FigureHandle1)


