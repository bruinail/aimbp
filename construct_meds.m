function [meds] = construct_meds(med_names,doses,Emax,EC50)

if(~iscell(med_names))
    med_names = {med_names};
    doses = {doses};
end
for i = 1:length(med_names)
    med = med_info(med_names{i});
    med.doses = doses{i};
    med.Emax = Emax(i);
    med.EC50 = EC50(i);
    meds(i) = med;
end