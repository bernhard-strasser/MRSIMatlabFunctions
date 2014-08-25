function specpoint_vector = compute_inverse_chemshift_vector_1_0(chemshift_vec,water_frequency,dwelltime)




%for assign_index = 1:3
%	eval([ sprintf('input%d_1


%if(ischar(input1_2))
%	eval([ sprintf(%s = %s, sprintf('input%d_1',assign_index),sprintf('input%d_2',assign_index) ])
%elseif(isfloat(input1_2))
%	eval([ sprintf(%s = %s, input1_1,input1_2) ])
%end





T_total = (vecSize - 1) * dwelltime;  % Total measurement time; The first point measuring at time 0, second point measuring at time dwelltime, ... So vecSize'th point measures at time (vecSize - 1) * dwelltime

TMS_frequency = water_frequency/(1+4.65*10^-6);


step_chemshift = 10^6/(TMS_frequency * T_total);	

% step_chemical_shift = chemshift_1 - chemshift_2 = (freq1 - TMS_freq)/TMS_freq*10^6 - (freq2 - TMS_freq)/TMS_freq*10^6 =
% = (freq1-freq2)/TMS_freq*10^6 = 1/T_total * 10^6/TMS_freq
% where chemshift_1 and 2 and freq1 and 2 are to chemshifts/freq of adjoining points

% so we know now the difference in chem shift between 2 points. But what chem. shift is eg the 512'th point?
% the most left point is the highes frequency, that is freq_max = water_frequency + bandwidth/2 (bandwidth symmetric about water_frequency)
% This corresponds to a chem shift of
% max_chemshift = (freq_max - TMS_freq)/TMS_freq * 10^6  = (water_freq - TMS_freq + bandwidth/2) / TMS_freq * 10^6
% min_chemshift = (freq_min - TMS_freq)/TMS_freq * 10^6  = (water_freq - TMS_freq - bandwidth/2) / TMS_freq * 10^6

bandwidth_frequency = 1/dwelltime;
max_chemshift = (water_frequency - TMS_frequency + bandwidth_frequency/2)/(TMS_frequency)*10^6; % From equation: chemshift = max_chemshift - (n-1) * step_chemshift

specpoint_vector = round(-chemshift_vec + max_chemshift)/step_chemshift + 1);





