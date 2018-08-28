function [Phi] = mh_step_Phi(model, Phi_prior, X_0, X, Y, U)
% Because IA2RMS only works on a univariate level, we will estimate each
% entry in Phi separately.
Phi_sampled = false(size(Phi_prior.known));
Phi = model.Phi; % We'll make a copy of Phi that we'll work with and update

% Handle knowns
Phi_sampled(~isnan(Phi_prior.known)) = true;
Phi(~isnan(Phi_prior.known)) = Phi_prior.known(~isnan(Phi_prior.known));

% Handle r_h
if(isfield(Phi_prior,'r_h'))
    r_h_vals = sub2ind(size(Phi),[1 2 5],[3 4 6]);
    one_minus_r_h_vals = sub2ind(size(Phi),[1 2 5],[1 2 5]);
    % Sample first entry
    func = @(Phi_entry) Phi_func_sums(Phi_entry, r_h_vals, one_minus_r_h_vals, Phi, model, X_0, X, Y, U); % Function to sample
    % Do MH sampling with proposal distribution Uniform(Phi_prior.r_h)
    Phi_entry_prev = model.Phi(r_h_vals(1));
    Phi_entry_proposal = unifrnd(Phi_prior.r_h(1), Phi_prior.r_h(2));
    prob_prev = func(Phi_entry_prev);
    prob_prop = func(Phi_entry_proposal);
    accept_prob = exp(min(0,prob_prop - prob_prev));
    if rand() < accept_prob
        Phi_entry_sample = Phi_entry_proposal;
    else
        Phi_entry_sample = Phi_entry_prev;
    end
%         % DEBUG print
%         if(cs == 1)
%             fprintf('Phi Previous: %.2f / %.2f | Proposal: %.2f / %.2f | Prob: %.2f\n',Phi_entry_prev,func(Phi_entry_prev),Phi_entry_proposal,func(Phi_entry_proposal),accept_prob);
%         end
    for i = 1:length(r_h_vals)
        Phi_sampled(r_h_vals(i)) = true;
        Phi(r_h_vals(i)) = Phi_entry_sample;
    end
    % Take second entry as 1-first entry
    for i = 1:length(one_minus_r_h_vals)
        Phi_sampled(one_minus_r_h_vals(i)) = true;
        Phi(one_minus_r_h_vals(i)) = 1 - Phi_entry_sample;
    end
end

end

function [P] = Phi_func_sums(Phi_entry, inds, one_minus_inds, Phi, model, X_0, X, Y, U)

P = -inf(size(Phi_entry));

for i = 1:length(Phi_entry)
    % Since we have a limit on [0,1] for the constant sum entry, we'll
    % enforce it here to simplify things
    if(Phi_entry(i) >= 0 && Phi_entry(i) <= 1)
        Phi(inds) = Phi_entry(i);
        Phi(one_minus_inds) = 1 - Phi_entry(i);
        
        model.Phi = Phi;

        Qinv = pinv(model.Q);
        X_prev = [X_0 X(:,1:end-1)];
        X_diff = X - model.Phi*X_prev - model.Upsilon*U;
        log_cond = -0.5*trace(X_diff'*Qinv*X_diff);

        debug = false;
        if(debug)
            [~,X_gen,~] = ssm_gen(model,U,100);
            X_calc = Phi*X_prev + model.Upsilon*U;
            figure;
            hold on;
            plot(0:100,[X_0(5) X(5,:)]);
            plot(1:100,X_calc(5,:));
            plot(1:100,X_gen(5,:));
            plot(1:100,Y(3,:));
            legend('X','X_{calc}','X_{gen}','Y');
        end
        
        % Since we're not really accounting for differences between A*X and
        % Y elsewhere (e.g. we're holding A and R pretty close to constant)
        % we'll enforce that as a prior here
%         diff = (Y - model.A*X);
%         diff(isnan(diff)) = 0;
%         keep_Y = logical(model.A*keep);
%         diff(~keep_Y,:) = 0;
%         Rinv = pinv(model.R);
%         prior = -0.5*trace(diff'*Rinv*diff);
        
        P(i) = log_cond;% + prior;
    end
end

end