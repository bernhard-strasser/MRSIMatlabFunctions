function[wdth]=fwhm_v2(x,y,SearchPeakFromTo_ppm)


y_max=max(y);
maxpos = find(y==y_max);
% PeakRegion_LowerPt = find(min(abs(x - 1.9)) == abs(x - 1.9));
PeakRegion_LowerPt = find(min(abs(x - max(SearchPeakFromTo_ppm))) == abs(x - max(SearchPeakFromTo_ppm)));
PeakRegion_UpperPt = find(min(abs(x - min(SearchPeakFromTo_ppm))) == abs(x - min(SearchPeakFromTo_ppm)));


    try
        data1(:,1) = x(PeakRegion_LowerPt:maxpos-1);
        data1(:,2) = y(PeakRegion_LowerPt:maxpos-1);
        data2(:,1) = x(maxpos+1:PeakRegion_UpperPt);
        data2(:,2) = y(maxpos+1:PeakRegion_UpperPt);   
    catch wdth = NaN;
    end    


    try
        x_halfmax1 = interp1(data1(:,2),data1(:,1),(y_max./2));
        x_halfmax2 = interp1(data2(:,2),data2(:,1),(y_max./2));

        wdth = x_halfmax1-x_halfmax2;
    catch 
        wdth = NaN;
    end
