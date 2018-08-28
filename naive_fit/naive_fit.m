function [naive_model] = naive_fit(Y, M, init_params, priors)
init_Emax_EC50 = [init_params.U.Emax'; init_params.U.EC50'];
init_Emax_EC50 = init_Emax_EC50(:)';
x0 = [init_params.mu_0(1)-init_params.mu_0(3),init_params.U.r_B,init_params.mu_0(3),init_Emax_EC50];
lb = [-400,0,0];
ub = [400,1,400];
for i = 1:length(priors)
    lb = [lb, priors(i).Emax_range(1), priors(i).EC50_range(1)];
    ub = [ub, priors(i).Emax_range(2), priors(i).EC50_range(2)];
end
xdata.M = M;
xdata.is_missing = isnan(Y(1,:));
x = lsqcurvefit(@F,x0,xdata,Y(1,~xdata.is_missing),lb,ub);
naive_model.Bmax = x(1);
naive_model.r_B = x(2);
naive_model.sbp_baseline_mu_0 = x(3);
naive_model.Emax = x(4:2:end-1)';
naive_model.EC50 = x(5:2:end)';
end

function [Y] = F(x, xdata)
Bmax = x(1);
r_B = x(2);
M = xdata.M;
sbp_baseline_mu_0 = x(3);
n_drugs = size(M,1);
Emax = zeros(n_drugs,1);
EC50 = zeros(n_drugs,1);
for i = 1:n_drugs
    Emax(i) = x(3+2*(i-1)+1);
    EC50(i) = x(3+2*(i-1)+2);
end
n = size(M,2);
Y = sbp_baseline_mu_0 + stroke_perturbation(Bmax,r_B,n) + sum(drug_emax_model(Emax,EC50,M),1);
Y = Y(~xdata.is_missing);
end