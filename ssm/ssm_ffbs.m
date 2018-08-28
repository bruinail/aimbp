function [X, X_0] = ssm_ffbs(model, Y, U)
%Forward filtering backward sampling algorithm

% From Shumway and Stoffer TSA4 Chapter 6 State Space Models 6.195

Tn = size(Y,2);
Xn = size(model.mu_0,1);

[X_filter, P_filter, X_exp_prev, P_exp_prev, K, ~, ~, loglike] = ssm_kalman_filter(model, Y, U);

% figure;
% hold on;
% plot(1:100,X_filter(1,:));
% plot(1:100,X_filter(3,:));
% plot(1:100,Y(1,:));
% legend('X','baseline','Y');

X = zeros(Xn,Tn);

X(:,end) = mvnrnd(X_filter(:,end),correct_cov(P_filter(:,:,end)));

for i = (Tn-1):-1:1
    J = P_filter(:,:,i)*model.Phi'*pinv(P_exp_prev(:,:,i+1));
    m = X_filter(:,i) + J*(X(:,i+1) - X_exp_prev(:,i+1));
    V = P_filter(:,:,i) - J*P_exp_prev(:,:,i+1)*J';
    X(:,i) = mvnrnd(m,correct_cov(V))';
end

J_0 = model.Sigma_0*model.Phi'*pinv(P_exp_prev(:,:,1));
m_0 = model.mu_0 + J_0*(X(:,1) - X_exp_prev(:,1));
V_0 = model.Sigma_0 - J_0*P_exp_prev(:,:,1)*J_0';
X_0 = mvnrnd(m_0, correct_cov(V_0))';