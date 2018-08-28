function [Upsilon] = mh_step_Upsilon(model, Upsilon_prior, X_0, X, U)
% Because IA2RMS only works on a univariate level, we will estimate each
% entry in Upsilon separately.
Upsilon_sampled = false(size(Upsilon_prior.known));
Upsilon = model.Upsilon;

% Handle knowns
Upsilon_sampled(~isnan(Upsilon_prior.known)) = true;
Upsilon(~isnan(Upsilon_prior.known)) = Upsilon_prior.known(~isnan(Upsilon_prior.known));

% Handle normals
if(isfield(Upsilon_prior,'normals'))
    for cs = 1:size(Upsilon_prior.normals,1)
        % Sample first entry
        func = @(Upsilon_entry) Upsilon_func(Upsilon_entry, Upsilon_prior, cs, Upsilon, model, X_0, X, U); % Function to sample

        % Do MH sampling with proposal distribution as the normal prior
        % distribution
        Upsilon_entry_prev = model.Upsilon(Upsilon_prior.normals(cs,1),Upsilon_prior.normals(cs,2));
%         Upsilon_entry_proposal = unifrnd(-10,0);
        Upsilon_entry_proposal = normrnd(Upsilon_entry_prev,1);
        while(Upsilon_entry_proposal < -10 || Upsilon_entry_proposal > 0)
            Upsilon_entry_proposal = normrnd(Upsilon_entry_prev,1);
        end
        accept_prob = exp(min(0,func(Upsilon_entry_proposal) - func(Upsilon_entry_prev)));
%         % DEBUG print
%         if(cs == 1 || cs == 3)
%             fprintf('Upsilon Previous: %.2f / %.2f | Proposal: %.2f / %.2f | Prob: %.2f\n',Upsilon_entry_prev,func(Upsilon_entry_prev),Upsilon_entry_proposal,func(Upsilon_entry_proposal),accept_prob);
%         end
        if rand() < accept_prob
            Upsilon_entry_sample = Upsilon_entry_proposal;
        else
            Upsilon_entry_sample = Upsilon_entry_prev;
        end
        Upsilon_sampled(Upsilon_prior.normals(cs,1),Upsilon_prior.normals(cs,2)) = true;
        Upsilon(Upsilon_prior.normals(cs,1),Upsilon_prior.normals(cs,2)) = Upsilon_entry_sample;
    end
end
end

function [P] = Upsilon_func(Upsilon_entry, prior, prior_col, Upsilon, model, X_0, X, U)

Upsilon_entry_row = prior.normals(prior_col,1);
Upsilon_entry_col = prior.normals(prior_col,2);
prior_mean = prior.normals(prior_col,3);
prior_std = prior.normals(prior_col,4);

P = zeros(size(Upsilon_entry));

for i = 1:length(Upsilon_entry)
    Upsilon(Upsilon_entry_row,Upsilon_entry_col) = Upsilon_entry(i);
    
    Qinv = pinv(model.Q);
    X_prev = [X_0 X(:,1:end-1)];
    X_diff = X - model.Phi*X_prev - Upsilon*U;
    % Only consider the X_diff that our Upsilon entry actually affects
    keep = false(size(X_diff,1),1);
    keep(Upsilon_entry_row) = true;
    X_diff(~keep,:) = 0;
    log_cond = diag(X_diff'*Qinv*X_diff);
    % We want to weight the log_cond by the proportion by which this
    % particular Upsilon entry affects log_cond compared to other Upsilon
    % entries.  Can't weight by Upsilon itself, because then we'll just
    % drive Upsilon to zero.  Instead we'll use just the prior for Upsilon
    % as an estimate of the effect.
    Upsilon_est = Upsilon;
    Upsilon_est(sub2ind(size(Upsilon),prior.normals(:,1),prior.normals(:,2))) = prior.normals(:,3);
    U_weight = (abs(Upsilon_est(Upsilon_entry_row,Upsilon_entry_col))*U(Upsilon_entry_col,:))./(abs(Upsilon_est(Upsilon_entry_row,:))*U);
    U_weight(isnan(U_weight)) = 0;
    log_cond = log_cond.*U_weight';
%     % Accumulate log_cond
%     log_cond = cumsum(log_cond);
    log_cond = -0.5*sum(log_cond);
    
    debug = false;
    if(debug)
        model.Upsilon = Upsilon;
        [~,X_gen,~] = ssm_gen(model,U,100);
        X_calc = model.Phi*X_prev + Upsilon*U;
        figure;
        hold on;
        plot(1:100,X(1,:));
        plot(1:100,X_calc(1,:));
        plot(1:100,X_gen(1,:));
        legend('X','X_{calc}','X_{gen}');
    end

    % Multiply by the prior; i.e. add log prior
    P(i) = log_cond + log(normpdf(Upsilon_entry(i),prior_mean,prior_std));
end

end