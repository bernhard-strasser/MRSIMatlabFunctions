function chemshift_vector = compute_chemshift_vector_1_1(water_frequency,dwelltime,vecSize)
% In units of:
% water_frequency: Hz
% dwelltime: s
% vecSize: -

%T_total = (vecSize - 1) * dwelltime;  % Total measurement time; The first point measuring at time 0, second point measuring at time dwelltime, ... So vecSize'th point measures at time (vecSize - 1) * dwelltime

TMS_frequency = water_frequency/(1+4.65*10^-6);


%step_chemshift = 10^6/(TMS_frequency * T_total);	

% step_chemical_shift = chemshift_1 - chemshift_2 = (freq1 - TMS_freq)/TMS_freq*10^6 - (freq2 - TMS_freq)/TMS_freq*10^6 =
% = (freq1-freq2)/TMS_freq*10^6 = 1/T_total * 10^6/TMS_freq
% where chemshift_1 and 2 and freq1 and 2 are to chemshifts/freq of adjoining points


% so we know now the difference in chem shift between 2 points. But what chem. shift is eg the 512'th point?
% the most left point is the highes frequency, that is freq_max = water_frequency + bandwidth/2 (bandwidth symmetric about water_frequency)
% This corresponds to a chem shift of
% max_chemshift = (freq_max - TMS_freq)/TMS_freq * 10^6  = (water_freq - TMS_freq + bandwidth/2) / TMS_freq * 10^6
% min_chemshift = (freq_min - TMS_freq)/TMS_freq * 10^6  = (water_freq - TMS_freq - bandwidth/2) / TMS_freq * 10^6

bandwidth_frequency = 1/dwelltime;

max_chemshift = (water_frequency - TMS_frequency + bandwidth_frequency/2)/(TMS_frequency)*10^6;
min_chemshift = (water_frequency - TMS_frequency - bandwidth_frequency/2)/(TMS_frequency)*10^6;
step_chemshift = (max_chemshift - min_chemshift)/(vecSize-1);

% The chemichal shift is now from min_chemshift to max_chemshift with steps step_chemshift

chemshift_vector = max_chemshift:-step_chemshift:min_chemshift; % High chemshift means small point because the FID got conjugated so that the spectrum is flipped left/right
