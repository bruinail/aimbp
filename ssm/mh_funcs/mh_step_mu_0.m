function [mu_0] = mh_step_mu_0(model, mu_0_prior, X_0, X, U)

% Use the Gibbs step first to get samples for entries of mu_0 that are not
% constant
var_new = correct_cov(pinv(pinv(model.Sigma_0) + pinv(mu_0_prior.S)));
m_new = (var_new*(pinv(model.Sigma_0)*X_0+pinv(mu_0_prior.S)*mu_0_prior.m));
mu_0 = mvnrnd(m_new',var_new)';

% Because we hold certain rows of X_0 constant, we need to estimate
% mu_0 for those rows separately otherwise we'll never move away
% from the initial value

% Handle constant values
if(isfield(mu_0_prior,'constants'))
    for cs = 1:size(mu_0_prior.constants,1)
        % Sample first entry
        func = @(mu_0_entry) mu_0_func(mu_0_entry, mu_0_prior.constants(cs,:), mu_0, model, X_0, X, U); % Function to sample

        % Do MH sampling with proposal distribution as the normal prior
        % distribution, with limits of +/- 3 stdev from priors
        prior_mean = mu_0_prior.constants(cs,3);
        prior_std = mu_0_prior.constants(cs,4);
        mu_0_entry_prev = model.mu_0(mu_0_prior.constants(cs,1));
        mu_0_entry_proposal = normrnd(mu_0_entry_prev,10);
        while(mu_0_entry_proposal > prior_mean+3*prior_std || mu_0_entry_proposal < prior_mean-3*prior_std)
            mu_0_entry_proposal = normrnd(mu_0_entry_prev,10);
        end
        prob_prev = func(mu_0_entry_prev);
        prob_prop = func(mu_0_entry_proposal);
        accept_prob = exp(min(0,prob_prop - prob_prev));
%         % DEBUG print
%         if(cs == 1)
%             fprintf('mu_0 Previous: %.2f / %.2f | Proposal: %.2f / %.2f | Prob: %.2f\n',mu_0_entry_prev,func(mu_0_entry_prev),mu_0_entry_proposal,func(mu_0_entry_proposal),accept_prob);
%         end
        if rand() < accept_prob
            mu_0_entry_sample = mu_0_entry_proposal;
        else
            mu_0_entry_sample = mu_0_entry_prev;
        end
        mu_0(mu_0_prior.constants(cs,1)) = mu_0_entry_sample;
    end
end
end

function [P] = mu_0_func(mu_0_entry, constants_prior, mu_0, model, X_0, X, U)

mu_0_entry_row = constants_prior(1);
X_row = constants_prior(2);
prior_mean = constants_prior(3);
prior_std = constants_prior(4);

P = zeros(size(mu_0_entry));

for i = 1:length(mu_0_entry)
    mu_0(mu_0_entry_row) = mu_0_entry(i);
    X_0(mu_0_entry_row) = mu_0(mu_0_entry_row);
    X(mu_0_entry_row,:) = mu_0(mu_0_entry_row);
    
    Qinv = pinv(model.Q);
    X_prev = [X_0 X(:,1:end-1)];
    X_diff = X - model.Phi*X_prev - model.Upsilon*U;
    % Only consider the X_diff that our mu_0 entry actually affects
    keep = false(size(X_diff,1),1);
    keep(X_row) = true;
    X_diff(~keep,:) = 0;
    % Weight the X_diff by how little r_B and thus U_L would have
    % interfered with baseline estimation
    X_diff = X_diff.*repmat(1-abs(U(1,:))/max(abs(U(1,:))),size(X_diff,1),1);
    log_cond = diag(X_diff'*Qinv*X_diff);
    log_cond = -0.5*sum(log_cond);
    
    debug = false;
    if(debug)
        model.mu_0 = mu_0;
        [~,X_gen,~] = ssm_gen(model,U,100);
        X_calc = model.Phi*X_prev + model.Upsilon*U;
        figure;
        hold on;
        plot(0:100,[X_0(1) X(1,:)]);
        plot(1:100,X_calc(1,:));
        plot(1:100,X_gen(1,:));
        legend('X','X_{calc}','X_{gen}');
    end

    % Multiply by the prior; i.e. add log prior
    P(i) = log_cond + log(normpdf(mu_0_entry(i),prior_mean,prior_std));
end

end