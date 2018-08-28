function [model] = ssm_initialize(init_params)
% Initializes an SSM model $\Theta = \{\mu_0, \Sigma_0, \Phi, A, Q, R, 
% \Upsilon\}$.

% $\mu_0$ and $\Sigma_0$ are the initial mean and covariance of $x_0$.
% These should be passed in as initial parameters for what we think $X$ is
% like.
model.mu_0 = init_params.mu_0;

model.Sigma_0 = init_params.Sigma_0;

% $\Phi$ is the transition matrix for X.
if(isfield(init_params,'Phi'))
    model.Phi = init_params.Phi;
else
    model.Phi = rand(init_params.Xn, init_params.Xn);
end

% $A$ is the emission matrix from X to Y.
if(isfield(init_params,'A'))
    model.A = init_params.A;
else
    model.A = rand(init_params.Yn, init_params.Xn);
end

% $Q$ is the covariance matrix for $w$, the noise added to $X$.
if(isfield(init_params,'Q'))
    model.Q = init_params.Q;
else
    model.Q = init_params.Sigma_0;
end

% $R$ is the covariance matrix for $v$, the noise added to $Y$.
if(isfield(init_params,'R'))
    model.R = init_params.R;
else
    model.R = eye(init_params.Yn, init_params.Yn);
end

% $\Upsilon$ is the input matrix for $X$.
if(isfield(init_params,'Upsilon'))
    model.Upsilon = init_params.Upsilon;
else
    model.Upsilon = zeros(init_params.Xn, init_params.Un);
end

model.U.Bmax = init_params.U.Bmax;
model.U.r_B = init_params.U.r_B;
model.U.Emax = init_params.U.Emax;
model.U.EC50 = init_params.U.EC50;