function [model] = simulate_stroke_model(initial,r_B,interventions,add_art_line,noise_scale)
% Constructs AIM-BP model

% Rates at which HR and BP converge to homeostasis
rates.hr = 0.9;
rates.bp = 0.9;

%% Simulating heart rate
% To simulate heart rate, we do something similar to BP where the
% body temp also wants to trend towards a homeostatic number. The variables
% are then [hr; hr homeostasis]

% Initial values for heart rate and homeostasis baseline
hr_mu_0 = initial.hr;
hr_Sigma_0 = diag([0 0].^2);

% Transition matrix (To keep homeostasis variable constant, last row must
% be identity
hr_transition = [(1-rates.hr) rates.hr; 0 1];
% Covariance matrix
hr_covariance = (noise_scale.*[4 0; 0 0]).^2;

% Our observed heart rate is a linear function of the actual heart rate
% Emission matrix
hr_emission = [1 0];
% Observation noise covariance
hr_obs_covariance = (noise_scale*3)^2;

%% Simulating BP
% To simulate BP, we take a model where the natural BP wants to trend
% towards a homeostatic number. The variables are [systolic BP; diastolic
% BP; systolic BP homeostasis; diastolic BP homeostasis]

% Initial values for BP and homeostasis baseline
bp_mu_0 = initial.bp;
bp_Sigma_0 = diag([0 0 0 0].^2);

% Transition matrix (To keep homeostasis variables constant, last two rows
% must be identity matrix)
bp_transition = [(1-rates.bp) 0 rates.bp 0; 0 (1-rates.bp) 0 rates.bp; 0 0 1 0; 0 0 0 1];
% Covariance matrix
bp_covariance = (noise_scale.*[9 0 0 0; 0 7 0 0; 0 0 0 0; 0 0 0 0]).^2;

% Our observed blood pressure is a linear function of the
% systolic/diastolic BP

% If we're adding an arterial line add the following lines to emission and
% covariance
if(add_art_line)
    % Emission matrix
    bp_emission = [1 0 0 0; 0 1 0 0; 1 0 0 0; 0 1 0 0];
    % Observation noise covariance
    bp_obs_covariance = (noise_scale.*[5 0 0 0; 0 4 0 0; 0 0 1 0; 0 0 0 1]).^2;
else
    % Emission matrix
    bp_emission = [1 0 0 0; 0 1 0 0];
    % Observation noise covariance
    bp_obs_covariance = (noise_scale.*[5 0; 0 4]).^2;
end



%% Combine all processes into a single SSM

% We want to include a constant term at the end to model non-zero mean
% measurement noise
model.mu_0 = [bp_mu_0; hr_mu_0; 1];
model.U.Bmax = bp_mu_0(1) - bp_mu_0(3);
model.U.r_B = r_B;
model.Sigma_0 = blkdiag(bp_Sigma_0,hr_Sigma_0,0);

model.Phi = blkdiag(bp_transition,hr_transition,1);
model.Q = blkdiag(bp_covariance,hr_covariance,0);

if(add_art_line)
    A_bias = [5; 5; 0; 0; 0];
else
    A_bias = [0; 0; 0];
end
model.A = [blkdiag(bp_emission,hr_emission) A_bias];
model.R = blkdiag(bp_obs_covariance,hr_obs_covariance);

% Add dependency of blood pressure on heart rate
model.Phi(1,5) = 0;
model.Phi(2,5) = 0;

% Construct Upsilon and fill out rest of model.U
n_drugs = length(interventions.meds);
model.U.Emax = zeros(n_drugs,1);
model.U.EC50 = zeros(n_drugs,1);
model.Upsilon = zeros(length(model.mu_0),1+n_drugs);
SBP_DBP_ratio = (bp_mu_0(2) - bp_mu_0(4))/model.U.Bmax;
% Stroke perturbation of SBP and DBP
model.Upsilon(1,1) = 1*rates.bp;
model.Upsilon(2,1) = SBP_DBP_ratio*rates.bp;
% Input from meds
for i = 1:n_drugs
    if(strcmp(interventions.meds(i).affects,'bp'))
        model.Upsilon(1,1+i) = 1*rates.bp;
        model.Upsilon(2,1+i) = SBP_DBP_ratio*rates.bp;
    end
    model.U.Emax(i) = interventions.meds(i).Emax;
    model.U.EC50(i) = interventions.meds(i).EC50;
end