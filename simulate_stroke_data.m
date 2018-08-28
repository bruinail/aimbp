function [Y, M, X, X_0] = simulate_stroke_data(model,Tn,interventions)
% Simulates stroke data for Tn discrete time steps.
% Other inputs:
% model: SSM model describing the patient
% interventions: Struct with interventions, including IV and oral

%% Generate hidden and observed data

U_L = stroke_perturbation(model.U.Bmax,model.U.r_B,Tn);
n_drugs = length(interventions.meds);
M = zeros(n_drugs,Tn);
doses = zeros(n_drugs,Tn);

% If we have it, calculate BP oral meds
for j = 1:n_drugs
    if(~interventions.meds(j).titratable)
        doses(j,:) = interventions.meds(j).doses;
        M(j,:) = calc_drug_conc(doses(j,:),interventions.meds(j).onset,interventions.meds(j).halflife,interventions.meds(j).dose_to_max_conc);
    end
end

X_0 = mvnrnd(model.mu_0',model.Sigma_0)';
X = zeros(size(X_0,1),Tn);

maintain_adjust = 0;
maintaining = false;
X(:,1) = model.Phi*X_0 + model.Upsilon*[U_L(:,1); drug_emax_model(model.U.Emax,model.U.EC50,M(:,1))] + mvnrnd(zeros(size(X_0))',model.Q)';
for i = 2:Tn
    % If we have a BP IV med and the BP is still too high, increase dose
    % from previous unless we're at max.  If the BP is too low, reduce dose
    % unless we're at 0. Start after 10 time points.
    if(i > 10)
        for j = 1:n_drugs
            if(interventions.meds(j).titratable)
                if(maintaining)
                    doses(j,i) = max(min(interventions.meds(j).titrate_maintain + maintain_adjust*interventions.meds(j).titrate_step, interventions.meds(j).titrate_max),0);
                    if(all(X(1,max(1,i-4):i-1) > 140))
                        maintain_adjust = maintain_adjust + 1;
                    elseif(all(X(1,max(1,i-4):i-1) < 120))
                        maintain_adjust = maintain_adjust - 1;
                    end
                else
                    if(all(X(1,max(1,i-4):i-1) > 140))
                        if(doses(j,i-1) == 0)
                            doses(j,i) = interventions.meds(j).titrate_init;
                        else
                            doses(j,i) = min(doses(j,i-1) + interventions.meds(j).titrate_step, interventions.meds(j).titrate_max);
                        end
                    else
                        doses(j,i) = interventions.meds(j).titrate_maintain;
                        maintaining = true;
                    end
                end
                    
                M(j,:) = calc_drug_conc(doses(j,:),interventions.meds(j).onset,interventions.meds(j).halflife,interventions.meds(j).dose_to_max_conc);
            end
        end
    end
    X(:,i) = model.Phi*X(:,i-1) + model.Upsilon*[U_L(:,i); drug_emax_model(model.U.Emax,model.U.EC50,M(:,i))] + mvnrnd(zeros(size(X_0))',model.Q)';
end

Y = model.A*X + mvnrnd(zeros(Tn,size(model.A,1)),model.R)';