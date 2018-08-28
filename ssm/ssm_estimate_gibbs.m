function [samples] = ssm_estimate_gibbs(Y, M, init_params, priors, n_burnin, n_samples, n_skip)
% MCMC estimation of model parameters $\Theta = 
% \{\mu_0, \Sigma_0, \Phi, A, Q, R, \Upsilon\}$ from data $Y$ and
% input $U$.
% init_params: Parameters for initializing a model using ssm_initialize.
% priors: Struct containing priors for each variable as follows (see
% doc_ssm_priors for more information):
%   mu_0: mu_0.m, mu_0.S
%   Sigma_0: Sigma_0.T, Sigma_0.nu
%   Phi: Phi.known, Phi.sums, Phi.norm_m, Phi.norm_S

n_iter = n_burnin + n_samples*n_skip;
n = size(Y,2);

Emax_estimation_method = false(size(M,1),1);
for i = 1:size(M,1)
    min_non_zero = M(i,:);
    min_non_zero = min(min_non_zero(min_non_zero > 0));
    if(isempty(min_non_zero))
        Emax_estimation_method(i) = -1;
    else
        M_hist = histcounts(M(i,:),10,'BinLimits',[min_non_zero, max(M(i,:))]);
        Emax_estimation_method(i) = max(M_hist) > 1*sum(M_hist);
    end
end

Upsilon_affects = init_params.Upsilon > 0;

U_all = zeros(2,0);
ll_all = zeros(1,0);
r_B_all = zeros(1,0);

samples = struct([]);

% Find initial parameters for the model $\Theta$.
[model_init] = ssm_initialize(init_params);

j = 1;
model = model_init;
for i = 1:n_iter
    if(mod(i,100) == 0)
        fprintf('Iteration %d...\n',i);
    end
    
    U = make_U(model,M);
    
    % We get estimated X from the Forward Filtering Backward Smoothing
    % algorithm
    [X, X_0] = ssm_ffbs(model, Y, U);
    X_prev = [X_0 X(:,1:end-1)];
    
    % We then simulate a new model conditioned on current X using Gibbs and
    % MH or IA2RMS for various variables
    
    % mu_0 - Gibbs step
    if(isfield(priors.mu_0,'known'))
        model.mu_0 = priors.mu_0.known;
    else
        model.mu_0 = mh_step_mu_0(model, priors.mu_0, X_0, X, U);
    end
    % Update Bmax and Upsilon
    model.U.Bmax = model.mu_0(1) - model.mu_0(3);
    SBP_DBP_ratio = (model.mu_0(2) - model.mu_0(4))/model.U.Bmax;
    for Ups_c = 1:size(Upsilon_affects,2)
        if(Upsilon_affects(1,Ups_c))
            model.Upsilon(1,Ups_c) = model.Phi(1,3);
        end
        if(Upsilon_affects(2,Ups_c))
            model.Upsilon(2,Ups_c) = model.Phi(2,4)*SBP_DBP_ratio;
        end
        % Make sure that the only drugs affecting HR have effect values M
        % that are the same scale as HR, not BP
        if(Upsilon_affects(5,Ups_c))
            model.Upsilon(5,Ups_c) = model.Phi(5,6);
        end
    end
    
    % Sigma_0 - Gibbs step
    if(isfield(priors.Sigma_0,'known'))
        model.Sigma_0 = priors.Sigma_0.known;
    else
        model.Sigma_0 = iwishrnd(priors.Sigma_0.T + (X_0 - model.mu_0)*(X_0 - model.mu_0)', priors.Sigma_0.nu + 1);
    end
    
    % Phi - MH step
    model.Phi = mh_step_Phi(model, priors.Phi, X_0, X, Y, U);
    
    % A - MH step
    model.A = mh_step_A(model, priors.A, X, Y);
    
    % Q - Gibbs step
    if(isfield(priors.Q,'known'))
        model.Q = priors.Q.known;
    else
        diff = (X - model.Phi*X_prev - model.Upsilon*U);
        model.Q = iwishrnd(priors.Q.T + diff*diff', priors.Q.nu + n);
        for k = 1:length(priors.Q.hold_zero)
            model.Q(:,priors.Q.hold_zero(k)) = 0;
            model.Q(priors.Q.hold_zero(k),:) = 0;
        end
    end
    
    % R - Gibbs step
    if(isfield(priors.R,'known'))
        model.R = priors.R.known;
    else
        diff = (Y - model.A*X);
        diff(isnan(diff)) = 0; % For the missing valued Y, set diff to 0
        model.R = iwishrnd(priors.R.T + diff*diff', priors.R.nu + n);
        model.R
    end
    
    % Upsilon - MH step
%     [model.Upsilon] = mh_step_Upsilon(model, priors.Upsilon, X_0, X, U);
    
    % Estimate parameters that govern U
    [model.U,ll] = mh_step_U(model, priors.U, X_0, X, U, M, Emax_estimation_method);
    ll_all = [ll_all ll];
    r_B_all = [r_B_all model.U.r_B];
    U_all = [U_all [model.U.Emax(1); model.U.EC50(1)]];
    
    % Save if we're on a post-burnin post-skip sample
    if(i == n_burnin + n_skip*j)
        samples(j).X = X;
        samples(j).X_0 = X_0;
        samples(j).model = model;
        j = j + 1;
    end
end

% figure;
% Emax_range = priors.U.Emax_model_params(1).Emax_range;
% EC50_range = priors.U.Emax_model_params(1).EC50_range;
% plot_Emax_EC50_hist(U_all',Emax_range,EC50_range);
% 
% figure;
% histogram(U_all(1,:)./(U_all(2,:)+mean(M(M(1,:) > 0))),100);

fprintf('Done.\n');