function[wdth]=fwhm(x,y,pts)


y_max=max(y);
maxpos = find(y==y_max);
% NAARegion_SP = find(min(abs(x - 1.9)) == abs(x - 1.9));
NAARegion_SP = find(min(abs(x - 2.3)) == abs(x - 2.3));


    if pts<2000
        try
            data1(:,1) = x(NAARegion_SP:maxpos-1);
            data1(:,2) = y(NAARegion_SP:maxpos-1);
            data2(:,1) = x(maxpos+1:size(y,1)-1);
            data2(:,2) = y(maxpos+1:size(y,1)-1);   
        catch wdth = NaN;
        end    
    else
        try
            data1(:,1) = x(NAARegion_SP:maxpos-1);
            data1(:,2) = y(NAARegion_SP:maxpos-1);
            data2(:,1) = x(maxpos+1:size(y,1)-1);
            data2(:,2) = y(maxpos+1:size(y,1)-1);
        catch wdth = NaN;
    end
    end

    try
        x_halfmax1 = interp1(data1(:,2),data1(:,1),(y_max./2));
        x_halfmax2 = interp1(data2(:,2),data2(:,1),(y_max./2));

        wdth = x_halfmax1-x_halfmax2;
    catch 
        wdth = NaN;
    end
