function [U_params,ll] = mh_step_U(model, U_prior, X_0, X, U, M, estimation_method)

% Handle knowns
model.U.Bmax = model.mu_0(1) - model.mu_0(3);
model.U.r_B(~isnan(U_prior.r_B.known)) = U_prior.r_B.known(~isnan(U_prior.r_B.known));
model.U.EC50(~isnan(U_prior.EC50.known)) = U_prior.EC50.known(~isnan(U_prior.EC50.known));
model.U.Emax(~isnan(U_prior.Emax.known)) = U_prior.Emax.known(~isnan(U_prior.Emax.known));

% Handle r_B
if(isnan(U_prior.r_B.known))
    func = @(r_B) r_B_func(r_B, model, X_0, X, M); % Function to sample

    % Do MH sampling with proposal distribution as the normal prior
    % distribution
    r_B_prev = model.U.r_B;
    r_B_proposal = normrnd(r_B_prev,0.1);
    while(r_B_proposal < 0.005 || r_B_proposal > 0.2)
        r_B_proposal = normrnd(r_B_prev,0.1);
    end
    prob_prev = func(r_B_prev);
    prob_prop = func(r_B_proposal);
    accept_prob = exp(min(0,prob_prop - prob_prev));
%     fprintf('r_B Prev: %.2f (%.2f) Prop: %.2f (%.2f) Accept: %.2f\n',r_B_prev,prob_prev,r_B_proposal,prob_prop,accept_prob);
    if rand() < accept_prob
        r_B_sample = r_B_proposal;
    else
        r_B_sample = r_B_prev;
    end
    model.U.r_B = r_B_sample;
end

% Handle Emax and EC50 simultaneously
if(isfield(U_prior,'Emax_model_params'))
    for n_drug = 1:length(U_prior.Emax_model_params)
        if(~isnan(U_prior.Emax.known(n_drug)) && ~isnan(U_prior.EC50.known(n_drug)))
            ll = 0;
            continue;
        end
        
        func = @(Emax_model_params) Emax_model_func(Emax_model_params, U_prior.Emax_model_params(n_drug), n_drug, model, X_0, X, M); % Function to sample

        % Do MH sampling with proposal distribution as the normal prior
        % distribution
        Emax_model_params_prev = [model.U.Emax(n_drug) model.U.EC50(n_drug)];
        reset_walk = unifrnd(0,1) < 0;
        if(reset_walk)
            Emax_model_params_mean = U_prior.Emax_model_params(n_drug).mu;
        else
            Emax_model_params_mean = Emax_model_params_prev;
        end
        
        % If no doses of this drug were recorded (M is all zero), doesn't
        % matter what we estimate so just keep previous.
        if(estimation_method(n_drug) == -1)
            Emax_model_params_proposal = Emax_model_params_prev;
        % If we have a relatively constant drug dosage, we can't
        % calculate both Emax and EC50 very well. Just hold EC50
        % constant and calculate Emax and report the ratio later.
        elseif(estimation_method(n_drug) == 1)
            Emax_model_params_proposal = [normrnd(Emax_model_params_mean(1),U_prior.Emax_model_params(n_drug).Sigma(2,2)^0.5) U_prior.Emax_model_params(n_drug).mu(2)];
            while(Emax_model_params_proposal(1) < U_prior.Emax_model_params(n_drug).Emax_range(1) || ...
                  Emax_model_params_proposal(1) > U_prior.Emax_model_params(n_drug).Emax_range(2))
                Emax_model_params_proposal = [normrnd(Emax_model_params_mean(1),U_prior.Emax_model_params(n_drug).Sigma(2,2)^0.5) U_prior.Emax_model_params(n_drug).mu(2)];
            end
        else
            Emax_model_params_proposal = mvnrnd(Emax_model_params_mean,U_prior.Emax_model_params(n_drug).Sigma./9);
            while(Emax_model_params_proposal(1) < U_prior.Emax_model_params(n_drug).Emax_range(1) || ...
                  Emax_model_params_proposal(1) > U_prior.Emax_model_params(n_drug).Emax_range(2) || ...
                  Emax_model_params_proposal(2) < U_prior.Emax_model_params(n_drug).EC50_range(1) || ...
                  Emax_model_params_proposal(2) > U_prior.Emax_model_params(n_drug).EC50_range(2))
                Emax_model_params_proposal = mvnrnd(Emax_model_params_mean,U_prior.Emax_model_params(n_drug).Sigma./9);
            end
        end
        prob_prev = func(Emax_model_params_prev);
        prob_prop = func(Emax_model_params_proposal);
        accept_prob = exp(min(0,prob_prop - prob_prev));
%         fprintf('Emax/EC50 Prev: %.2f/%.2f (%.2f) Prop: %.2f/%.2f (%.2f) Accept: %.2f\n',Emax_model_params_prev(1),Emax_model_params_prev(2),prob_prev,Emax_model_params_proposal(1),Emax_model_params_proposal(2),prob_prop,accept_prob);
        if rand() < accept_prob
            Emax_sample = Emax_model_params_proposal(1);
            EC50_sample = Emax_model_params_proposal(2);
            ll = prob_prop;
        else
            Emax_sample = Emax_model_params_prev(1);
            EC50_sample = Emax_model_params_prev(2);
            ll = prob_prev;
        end
        model.U.Emax(n_drug) = Emax_sample;
        model.U.EC50(n_drug) = EC50_sample;
    end
end

U_params = model.U;
end

function [P] = r_B_func(r_B, model, X_0, X, M)


U_L = stroke_perturbation(model.U.Bmax,r_B,size(M,2));
U_M = drug_emax_model(model.U.Emax,model.U.EC50,M);
U = [U_L; U_M];
Qinv = pinv(model.Q);
X_prev = [X_0 X(:,1:end-1)];
X_diff = X - model.Phi*X_prev - model.Upsilon*U;
% Only consider the X_diff that r_B actually affects
keep = model.Upsilon(:,1) > 0;
X_diff(~keep,:) = 0;
log_cond = diag(X_diff'*Qinv*X_diff);
log_cond = -0.5*sum(log_cond);

% Uniform prior so no change to log_cond
P = log_cond;

end

function [P] = Emax_model_func(Emax_model_params, prior, n_drug, model, X_0, X, M)

Emax = model.U.Emax;
EC50 = model.U.EC50;
Emax(n_drug) = Emax_model_params(1);
EC50(n_drug) = Emax_model_params(2);

U_L = stroke_perturbation(model.U.Bmax,model.U.r_B,size(M,2));
U_M = drug_emax_model(Emax,EC50,M);
U = [U_L; U_M];

Qinv = pinv(model.Q);
X_prev = [X_0 X(:,1:end-1)];
X_diff = X - model.Phi*X_prev - model.Upsilon*U;
% Only consider the X_diff that our Upsilon entry actually affects
keep = model.Upsilon(:,1+n_drug) > 0;
X_diff(~keep,:) = 0;
X_diff(:,M(n_drug,:) == 0) = 0;
% Weight the X_diff by how little r_B and thus U_L would have
% interfered with baseline estimation
% X_diff = X_diff.*repmat(1-U_L/max(U_L),size(X_diff,1),1);
log_cond = diag(X_diff'*Qinv*X_diff);
log_cond = -0.5*sum(log_cond);

% Multiply by the prior; i.e. add log prior
P = log_cond + log(mvnpdf([Emax(n_drug) EC50(n_drug)],prior.mu,prior.Sigma));

end