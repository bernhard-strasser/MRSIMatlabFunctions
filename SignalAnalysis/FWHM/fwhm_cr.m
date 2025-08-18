function[wdth]=fwhm(x,y,pts)

%min for 1.3 ms cre
%display([ char(10) char(10) 'TEST2'])
y_min=min(y);
maxpos = (find(y==y_min));
%y_min=-y_min;
    if pts<2000
        try
            data1(:,1) = x(maxpos-45:maxpos-1);
            data1(:,2) = y(maxpos-45:maxpos-1);
            data2(:,1) = x(maxpos+1:maxpos+46);
            data2(:,2) = y(maxpos+1:maxpos+46);   
        catch wdth = NaN;
        end    
    else
        try
            data1(:,1) = x(maxpos-15:maxpos-1);
            data1(:,2) = y(maxpos-15:maxpos-1);
            data2(:,1) = x(maxpos+1:size(y,1)-2);
            data2(:,2) = y(maxpos+1:size(y,1)-2);
        catch wdth = NaN;
        end
    end

    try
        x_halfmax1 = interp1(data1(:,2),data1(:,1),(y_min./2));
        x_halfmax2 = interp1(data2(:,2),data2(:,1),(y_min./2));

        wdth = x_halfmax1-x_halfmax2;
    catch 
        wdth = NaN;
    end
