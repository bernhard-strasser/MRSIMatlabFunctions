function IdF = compute_B0map_Hz(GradField_p)

Iph = double(dicomread(GradField_p));
Iph_info = read_ascconv(GradField_p);
dTE = abs(Iph_info.TEs(1)-Iph_info.TEs(2))/1e6; % 1e3 range in ms,1e6 range in s
Iph = (2*pi).*Iph/4096 - pi;
IdF = (Iph./(2*pi*dTE))';
