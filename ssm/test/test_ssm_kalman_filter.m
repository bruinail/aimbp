%% Generate DLM data

Tn = 100;

model.mu_0 = 1;
model.Sigma_0 = 2;
model.Phi = 1;
model.Upsilon = 0;
model.Q = 0.1;
model.A = 1;
model.R = 10;

U = repmat([0 0 0 0 1],1,Tn/5);

[Y,X,X_0] = ssm_gen(model,U,Tn);

[X_est, P_est, X_exp_prev, P_exp_prev, K, eps, Sigma_eps, loglike] = ssm_kalman_filter(model,Y,U);

%% Generate DLM data

Tn = 100;

model.mu_0 = [1; 1];
model.Sigma_0 = [1 0; 0 1];
model.Phi = [1 0; 0 1];
model.Upsilon = [0; 0];
model.Q = [1 0; 0 0.1];
model.A = [1 0];
model.R = 10;

U = repmat([0 0 0 0 1],1,Tn/5);

[Y,X,X_0] = ssm_gen(model,U,Tn);

[X_est, P_est, X_exp_prev, P_exp_prev, K, eps, Sigma_eps, loglike] = ssm_kalman_filter(model,Y,U);

%% Generate DLM data

Tn = 100;

model.mu_0 = 1;
model.Sigma_0 = 2;
model.Phi = 0.5;
model.Upsilon = 3;
model.Q = 1;
model.A = [0.5; 0.3];
model.R = [4 0; 0 2];

U = repmat([0 0 0 0 1],1,Tn/5);

[Y,X,X_0] = ssm_gen(model,U,Tn);

[X_est, P_est, X_exp_prev, P_exp_prev, K, eps, Sigma_eps, loglike] = ssm_kalman_filter(model,Y,U);
