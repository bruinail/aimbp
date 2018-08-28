function [X, P, X_0, P_0, P_lagone, loglike] = ssm_kalman_smoother(model, Y, U)
% Calculate most likely $X_i$ conditioned on $Y_n$ from model $\Theta = 
% \{\mu_0, \Sigma_0, \Phi, A, Q, R, \Upsilon\}$, data $Y$, and input
% $U$. Returns X, variances P, $X_0^n$ X_0, $P_0^n$ P_0, lag one covariance
% P_lagone, and negative log likelihood loglike

% From Shumway and Stoffer TSA4 Chapter 6 State Space Models Property 6.2
% and 6.3

Tn = size(Y,2);
Xn = size(model.mu_0,1);

[X_filter, P_filter, X_exp_prev, P_exp_prev, K, ~, ~, loglike] = ssm_kalman_filter(model, Y, U);
X = zeros(size(X_filter));
X(:,end) = X_filter(:,end);
P = zeros(size(P_filter));
P(:,:,end) = correct_cov(P_filter(:,:,end));

% Lag one covariance smoother $P_{t,t-1}^n$
P_lagone = zeros(size(P_filter));
K_curr = K(:,:,end);
P_lagone(:,:,end) = correct_cov((eye(Xn) - K_curr*model.A)*model.Phi*P_filter(:,:,end-1));

for i = (Tn-1):-1:1
    if(i < Tn-1)
        J_i = J_iminus1;
    else
        J_i = P_filter(:,:,i)*model.Phi'*pinv(P_exp_prev(:,:,i+1));
    end
    
    X(:,i) = X_filter(:,i) + J_i*(X(:,i+1) - X_exp_prev(:,i+1));
    P(:,:,i) = correct_cov(P_filter(:,:,i) + J_i*(P(:,:,i+1) - P_exp_prev(:,:,i+1))*J_i');
    
    if(i > 1)
        J_iminus1 = P_filter(:,:,i-1)*model.Phi'*pinv(P_exp_prev(:,:,i));
    else
        J_iminus1 = model.Sigma_0*model.Phi'*pinv(P_exp_prev(:,:,1));
    end
    P_lagone(:,:,i) = correct_cov(P_filter(:,:,i)*J_iminus1' + J_i*(P_lagone(:,:,i+1) - model.Phi*P_filter(:,:,i))*J_iminus1');
end

J_0 = J_iminus1;
X_0 = model.mu_0 + J_0*(X(:,1) - X_exp_prev(:,1));
P_0 = correct_cov(model.Sigma_0 + J_0*(P(:,:,1) - P_exp_prev(:,:,1))*J_0');