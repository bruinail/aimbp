function [E] = drug_emax_model(Emax, EC50, M)
% Calculates expected effects using the E_max model:
% Emax: A Mn x 1 vector of maximum effects for each of Mn drugs in M.
% EC50: A Mn x 1 vector of EC_50 half maximum concentrations for each of Mn
% drugs in M.
% M: A Mn x Tn matrix of drug plasma concentrations for each of Un drugs

Tn = size(M,2);
E = repmat(Emax,1,Tn).*M./(repmat(EC50,1,Tn)+M);
E(isnan(E)) = 0;