function [model, num_iters] = ssm_estimate_em(Y, U, init_params, constraints)
% EM estimation of model parameters $\Theta = 
% \{\mu_0, \Sigma_0, \Phi, A, Q, R, \Gamma, \Upsilon\}$ from data $Y$ and
% input $U$, with initial paramters init_params.

% NOTE: This is incomplete as it does not estimate Gamma and Upsilon.
max_iter = 100;

% From Shumway and Stoffer TSA4 Chapter 6 State Space Models page 341

% Find initial parameters for the model $\Theta$.
[model] = ssm_initialize(init_params);

n = size(Y,2);
lls = zeros(1,max_iter);
for i = 1:max_iter
    fprintf('Iteration %d...\n',i);
    % We run a Kalman smoother to obtain the initial estimates for $X_t^n$,
    % $P_t^n$, and the log likelihood.
    [X_smooth, P_smooth, X_smooth0, P_smooth0, P_lagone, loglike] = ssm_kalman_smoother(model, Y, U);
    lls(i) = loglike;
    
    %% Test for convergence
    if(i > 1 && -(lls(i) - lls(i-1)) < 1e-5)
        break;
    end
    
    %% Expectation step
    X_smooth_lag = [X_smooth0 X_smooth(:,1:end-1)];
    P_smooth_lag = cat(3, P_smooth0, P_smooth(:,:,1:end-1));
    S_11 = X_smooth*X_smooth' + sum(P_smooth,3);
    S_10 = X_smooth*X_smooth_lag' + sum(P_lagone,3);
    S_00 = X_smooth_lag*X_smooth_lag' + sum(P_smooth_lag,3);
    
    %% Maximization step
    model.Phi = S_10/S_00;
    model.Q = (S_11 - (S_10/S_00)*S_10')/n;
    model.R = 0;
    for t = 1:n
        err = Y(:,i) - model.A*X_smooth(:,i);
        model.R = model.R + err*err' + model.A*P_smooth(:,:,i)*model.A';
    end
    model.R = model.R/n;
    model.mu_0 = X_smooth0;
    model.Sigma_0 = P_smooth0;
end
num_iters = i;
fprintf('Done.\n');