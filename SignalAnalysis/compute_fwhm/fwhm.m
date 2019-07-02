function[wdth]=fwhm(data)

xx = data(:,1);
yy = data(:,2);
ymax=max(yy);
maxpos = find(yy==ymax);

minpos1 = find(min(abs(data(:,1)-2.1)) == abs(data(:,1)-2.1));
minpos2 = find(min(abs(data(:,1)-1.9)) == abs(data(:,1)-1.9));

data1(:,1) = xx(minpos1:maxpos);
data1(:,2) = yy(minpos1:maxpos);
data2(:,1) = xx(maxpos:minpos2);
data2(:,2) = yy(maxpos:minpos2);
for i = 2:size(data1,1)
    if data1(i-1,2) == data1(i,2)
        data1(i,2) = data1(i-1,2)+0.0001;
    end
end

for i= 2:size(data2,1)
    if data2(i-1,2) == data2(i,2)
        data2(i,2) = data2(i-1,2)-0.0001;
    end
end
if (size(data1,1) > 1) && (size(data2,1) > 1)
    x_halfmax1 = interp1(data1(:,2),data1(:,1),(ymax./2));
    x_halfmax2 = interp1(data2(:,2),data2(:,1),(ymax./2));
    wdth = x_halfmax1-x_halfmax2;
else wdth = NaN;
end

 end
