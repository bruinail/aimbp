function [Y,X,X_0] = ssm_gen(model,U,Tn)
% Takes SSM model $\Theta = \{\mu_0, \Sigma_0, \Phi, A, Q, R, 
% \Gamma, \Upsilon\}$ and generates hidden and observed data with input U
% for Tn time points.
X_0 = mvnrnd(model.mu_0',model.Sigma_0)';
X = zeros(size(X_0,1),Tn);
X(:,1) = model.Phi*X_0 + model.Upsilon*U(:,1) + mvnrnd(zeros(size(X_0))',model.Q)';
for i = 2:Tn
    X(:,i) = model.Phi*X(:,i-1) + model.Upsilon*U(:,i) + mvnrnd(zeros(size(X_0))',model.Q)';
end

Y = model.A*X + mvnrnd(zeros(Tn,size(model.A,1)),model.R)';