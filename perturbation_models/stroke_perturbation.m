function [B] = stroke_perturbation(Bmax, r_B, Tn)
% Calculates the amount of perturbation B to the homeostasis baseline due
% to the stroke. Bmax is the initial max amount of perturbation, r_B is
% the rate of decrease to 0 (where 1 is immediate and 0 is never), and Tn 
% is the number of time points to calculate. Since we assume X_0 is the
% full B_max, we start with Bmax*(1-r_B) at T=1.

B = cumprod(repmat(1-r_B,1,Tn));
B = Bmax*B;