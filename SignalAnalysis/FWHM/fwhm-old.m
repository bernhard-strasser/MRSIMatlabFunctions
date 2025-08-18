function[wdth]=fwhm(data)
    
x = data(:,1);
y = data(:,2);
ymax=max(y);
maxpos = find(y==ymax);
data1(:,1) = x(160:maxpos);
data1(:,2) = y(160:maxpos);
data2(:,1) = x(maxpos:size(y,1));
data2(:,2) = y(maxpos:size(y,1));
x_halfmax1 = interp1(data1(:,2),data1(:,1),(ymax./2));
x_halfmax2 = interp1(data2(:,2),data2(:,1),(ymax./2));
wdth = x_halfmax1-x_halfmax2;
 end
