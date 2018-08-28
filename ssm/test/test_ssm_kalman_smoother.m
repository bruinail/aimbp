%% Generate DLM data

Tn = 100;

model.mu_0 = 1;
model.Sigma_0 = 2;
model.Phi = 0.5;
model.Upsilon = 3;
model.Q = 0;
model.A = 0.5;
model.Gamma = 1.5;
model.R = 0;

U = repmat([0 0 0 0 1],1,Tn/5);

[Y,X,X_0] = ssm_gen(model,U,Tn);

[X_est, P_est, X_0, P_0, P_lagone, loglike] = ssm_kalman_smoother(model,Y,U);

%% Generate DLM data

Tn = 100;

model.mu_0 = [1; 1];
model.Sigma_0 = [1 0; 0 0];
model.Phi = [0.5 0.5; 0 1];
model.Upsilon = [3; 1];
model.Q = [1 0; 0 0];
model.A = [0.5 0.5];
model.Gamma = 1.5;
model.R = 0;

U = repmat([0 0 0 0 1],1,Tn/5);

[Y,X,X_0] = ssm_gen(model,U,Tn);

[X_est, P_est, X_0, P_0, P_lagone, loglike] = ssm_kalman_smoother(model,Y,U);

%% Generate DLM data

Tn = 100;

model.mu_0 = 1;
model.Sigma_0 = 2;
model.Phi = 0.5;
model.Upsilon = 3;
model.Q = 1;
model.A = [0.5; 0.3];
model.Gamma = [1.5; 1];
model.R = [4 0; 0 2];

U = repmat([0 0 0 0 1],1,Tn/5);

[Y,X,X_0] = ssm_gen(model,U,Tn);

[X_est, P_est, X_0, P_0, P_lagone, loglike] = ssm_kalman_smoother(model,Y,U);
