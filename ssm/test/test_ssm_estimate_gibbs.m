%% Generate DLM data

Tn = 100;

model.mu_0 = [5; 1];
model.Sigma_0 = [0.01 0; 0 0];
model.Phi = [0.5 0.5; 0 1];
model.Upsilon = [1; 0];
model.Q = [0.1 0; 0 0];
model.A = [0.5 0.5];
model.R = 0;

U = repmat([0 0 0 0 1],1,Tn/5);

[Y,X,X_0] = ssm_gen(model,U,Tn);

%% Construct priors

priors.mu_0.m = [4; 1];
priors.mu_0.S = [2 0; 0 1e-5];

priors.Sigma_0.T = [0.01 0; 0 1e-5].*10;
priors.Sigma_0.nu = 10;

priors.Phi.constant = [NaN NaN; 0 1];
priors.Phi.sums = [1,1, 1,2];

priors.Upsilon.constant = [1; 0];

priors.Q.T = [0.1 0; 0 1e-5].*10;
priors.Q.nu = 10;

priors.A.constant = [0.5 0.5];

priors.R.T = 1e-5*10;
priors.R.nu = 10;

%% Initial parameters

init_params.mu_0 = [3; 1];
init_params.Sigma_0 = [1 0; 0 0];
init_params.Phi = [0.1 0.9; 0 1];
init_params.Upsilon = [1; 0];
init_params.Q = [1 0; 0 0];
init_params.A = [0.5 0.5];
init_params.R = 0;

%% Test

n_burnin = 1;
n_samples = 1000;
n_skip = 1;

samples = ssm_estimate_gibbs(Y, U, init_params, priors, n_burnin, n_samples, n_skip);

%% Plot

Phi = zeros([size(model.Phi) length(samples)]);
for i = 1:length(samples)
    Phi(:,:,i) = samples(i).model.Phi;
end
Phi_one = squeeze(Phi(1,1,:));
figure;
subplot(2,1,1);
hist(Phi_one,10);
subplot(2,1,2);
plot(1:length(samples),Phi_one);