%% Test constant sums
% Generate DLM data

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

% Construct prior

Phi_prior.constant = [NaN NaN; 0 1];
Phi_prior.sums = [1,1, 1,2]; 

% Test

n = 2000;
model_iter = model;
model_iter.Phi = [0.1 0.9; 0 1];
Phi_samples = zeros([size(model.Phi) n]);
for i = 1:n
    Phi = mh_step_Phi(model_iter, Phi_prior, X_0, X, U);
    Phi_samples(:,:,i) = Phi;
    model_iter.Phi = Phi;
end

% Plot

Phi_one = squeeze(Phi_samples(1,1,:));
figure;
subplot(2,1,1);
hist(Phi_one,50);
subplot(2,1,2);
plot(1:n,Phi_one);

%% Test normals
% Generate DLM data

Tn = 100;

model.mu_0 = [5; 1];
model.Sigma_0 = [0.01 0; 0 0];
model.Phi = [0.5 0.5; 0 0.8];
model.Upsilon = [1; 0];
model.Q = [0.1 0; 0 0.00001];
model.A = [0.5 0.5];
model.R = 0;

U = repmat([0 0 0 0 1],1,Tn/5);

[Y,X,X_0] = ssm_gen(model,U,Tn);

% Construct prior

Phi_prior.constant = [0.5 0.5; 0 NaN];
Phi_prior.normals = [2,2, 1,2]; 

% Test

n = 2000;
model_iter = model;
model_iter.Phi = [0.5 0.5; 0 0.1];
Phi_samples = zeros([size(model.Phi) n]);
for i = 1:n
    fprintf('Sample %d...\n',i);
    Phi = mh_step_Phi(model_iter, Phi_prior, X_0, X, U);
    Phi_samples(:,:,i) = Phi;
    model_iter.Phi = Phi;
end

% Plot

Phi_one = squeeze(Phi_samples(2,2,:));
figure;
subplot(2,1,1);
hist(Phi_one,50);
subplot(2,1,2);
plot(1:n,Phi_one);