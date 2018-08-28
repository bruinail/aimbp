function [model, num_iters] = ssm_estimate_matlab(Y, U, init_params, constraints)
% Maximum likelihood estimation of parameters for state space model based
% on observed data Y and input U. The Matlab estimator uses the innovations
% form of the state space model, with Kalman gain K and innovations e(t),
% instead of the general model with noise w(t) and v(t).

%% Convert constraints into idss model structure
mdl = idss(init_params.Phi,init_params.Gamma,init_params.A,init_params.Upsilon);
mdl.Structure.A.Free = constraints.Phi;
mdl.Structure.B.Free = constraints.Gamma;
mdl.Structure.C.Free = constraints.A;
mdl.Structure.D.Free = constraints.Upsilon;
mdl.Structure.K.Free = true(size(init_params.Xn,init_params.Yn));

%% Set options
opts = ssestOptions('InitialState',init_params.mu_0);

%% Find optimal parameters for the model
[sys, x0] = ssest(iddata(Y',U',1),mdl,opts,'Ts',1);
model.Phi = sys.A;
model.Gamma = sys.B;
model.A = sys.C;
model.Upsilon = sys.D;
model.R = eye(size(sys.NoiseVariance));
model.Q = sys.K*sys.K';

num_iters = 0;