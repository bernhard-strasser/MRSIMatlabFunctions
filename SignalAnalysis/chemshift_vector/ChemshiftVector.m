function chemshift_vector = ChemshiftVector(LarmorFreq,dwelltime,vecSize)
%
% chemshift_vector = compute_chemshift_vector_1_2(LarmorFreq,dwelltime,vecSize)
% In units of:
% water_frequency: Hz
% dwelltime: s
% vecSize: -

water_frequency = LarmorFreq * (1 + 4.65 * 10^-6);
bandwidth_frequency = 1/dwelltime;

% CAUTION: Only if spectrum is centered around water frequency!!!
max_chemshift = (water_frequency - LarmorFreq + bandwidth_frequency/2)/(LarmorFreq)*10^6;
min_chemshift = (water_frequency - LarmorFreq - bandwidth_frequency/2)/(LarmorFreq)*10^6;
step_chemshift = (max_chemshift - min_chemshift)/(vecSize-1);

% The chemichal shift is now from min_chemshift to max_chemshift with steps step_chemshift
chemshift_vector = max_chemshift:-step_chemshift:min_chemshift; % High chemshift means small point because the FID got conjugated so that the spectrum is flipped left/right
