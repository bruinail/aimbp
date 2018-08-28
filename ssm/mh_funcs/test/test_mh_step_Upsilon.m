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

%% Construct prior

Upsilon_prior.constant = [NaN; 0];
Upsilon_prior.normals = [1,1, 0,20]; 

%% Test

n = 2000;
model_iter = model;
model_iter.Upsilon = [0.1; 0];
Upsilon_samples = zeros([size(model.Upsilon) n]);
for i = 1:n
    Upsilon = mh_step_Upsilon(model_iter, Upsilon_prior, X_0, X, U);
    Upsilon_samples(:,:,i) = Upsilon;
    model_iter.Upsilon = Upsilon;
end

%% Plot

Upsilon_one = squeeze(Upsilon_samples(1,1,:));
figure;
subplot(2,1,1);
hist(Upsilon_one,50);
subplot(2,1,2);
plot(1:n,Upsilon_one);