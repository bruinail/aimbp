function [conc] = calc_drug_conc(doses, onset, halflife, dose_to_max_conc)
% Calculate filter based on onset and halflife
decay = 0:(halflife*16);
filter = [(0:onset-1)/onset, exp(-(log(2)/halflife)*decay)];
conc = conv(doses*dose_to_max_conc,filter,'full');
conc = conc(1:length(doses));