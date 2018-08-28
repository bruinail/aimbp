function [meds] = med_info(med_names)
% Nicardipine: Cook E, Clifton GG, Vargas R, et al. Pharmacokinetics, pharmacodynamics, and minimum effective clinical dose of intravenous nicardipine. Clin Pharmacol Ther. 1990;47(6):706-718. http://www.ncbi.nlm.nih.gov/pubmed/2357865. Accessed July 16, 2018.
% Labetalol: Saotome T, Minoura S, Terashi K, Sato T, Echizen H, Ishizaki T. Labetalol in hypertension during the third trimester of pregnancy: its antihypertensive effect and pharmacokinetic-dynamic analysis. J Clin Pharmacol. 1993;33(10):979-988. http://www.ncbi.nlm.nih.gov/pubmed/8227470. Accessed July 5, 2018.
%            Lalonde RL, O’Rear TL, Wainer IW, Drda KD, Herring VL, Bottorff MB. Labetalol pharmacokinetics and pharmacodynamics: Evidence of stereoselective disposition. Clin Pharmacol Ther. 1990;48(5):509-519. doi:10.1038/clpt.1990.187.
% V_c, onset, halflife info: NIH PubChem
meds = struct('name',{},'affects',{},'titratable',{},'V_c',{},'onset',{},'halflife',{},'dose_to_max_conc',{},'titrate_init',{},'titrate_step',{},'titrate_max',{},'titrate_maintain',{},'prior',{});
if(~iscell(med_names))
    med_names = {med_names};
end
for i = 1:length(med_names)
    if(strcmp(med_names{i},'labetalol_iv'))
        med.name = med_names{i};
        med.affects = 'bp';
        med.titratable = false;
        med.V_c = 1.1;
        med.onset = 0;
        med.halflife = 24;
        med.dose_to_max_conc = 12.5;
        med.titrate_init = NaN;
        med.titrate_step = NaN;
        med.titrate_max = NaN;
        med.titrate_maintain = NaN;
        med.prior.mu = [-30 110];
        med.prior.Sigma = [900 0; 0 2500];
        med.prior.Emax_range = [-100 0];
        med.prior.EC50_range = [0 500];
        meds = [meds med];
    elseif(strcmp(med_names{i},'labetalol_oral'))
        med.name = med_names{i};
        med.affects = 'bp';
        med.titratable = false;
        med.V_c = 4.0;
        med.onset = 12;
        med.halflife = 24;
        med.dose_to_max_conc = 3.5;
        med.titrate_init = NaN;
        med.titrate_step = NaN;
        med.titrate_max = NaN;
        med.titrate_maintain = NaN;
        med.prior.mu = [-30 110];
        med.prior.Sigma = [900 0; 0 2500];
        med.prior.Emax_range = [-100 0];
        med.prior.EC50_range = [0 500];
        meds = [meds med];
    elseif(strcmp(med_names{i},'nicardipine_iv'))
        med.name = med_names{i};
        med.affects = 'bp';
        med.titratable = true;
        med.V_c = 0.3;
        med.onset = 6;
        med.halflife = 3;
        med.dose_to_max_conc = 5;
        med.titrate_init = 5;
        med.titrate_step = 2.5;
        med.titrate_max = 15;
        med.titrate_maintain = 3;
        med.prior.mu = [-30 70];
        med.prior.Sigma = [900 0; 0 900];
        med.prior.Emax_range = [-100 0];
        med.prior.EC50_range = [0 500];
        meds = [meds med];
    end
end