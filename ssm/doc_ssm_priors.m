%% Priors
% The priors for each parameter are described here, in terms of specific
% variable structures.

%% mu_0
% mu_0.m and mu_0.S are the mean and covariance matrix for the multivariate
% normal distribution

%% Sigma_0
% Sigma_0.T and Sigma_0.nu are the scale matrix and degrees of freedom of 
% the inverse Wishart

%% Phi
% Phi.constant: A matrix the size of Phi with the constant values in 
% entries that are known constants and NaN in entries that are not.
% Phi.sums: A n by 2 matrix with n pairs of indices that correspond to
% pairs of entries in Phi that share a constant sum constraint.
% Phi.norm_m and Phi.norm_S: For each entry of Phi that has a normal prior,
% these cover the mean and variance of that entry. NaN otherwise.
