function [cov] = correct_cov(cov)
% Corrects covariance by making it positive semidefinite

% Remove negatives on diagonals
cov(diag(diag(cov) < 0)) = 0;

% Make values close to zero actually zero
cov(abs(cov) < 1e-10) = 0;

% Make sure it's symmetric
cov = (cov+cov')/2;