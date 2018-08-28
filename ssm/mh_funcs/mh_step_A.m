function [A] = mh_step_A(model, A_prior, X, Y)
% Because IA2RMS only works on a univariate level, we will estimate each
% entry in A separately.
A_sampled = false(size(A_prior.known));
A = model.A;

% Handle knowns
A_sampled(~isnan(A_prior.known)) = true;
A(~isnan(A_prior.known)) = A_prior.known(~isnan(A_prior.known));

% Handle normals
if(isfield(A_prior,'normals'))
    for cs = 1:size(A_prior.normals,1)
        % Sample first entry
        func = @(A_entry) A_func(A_entry, A_prior.normals(cs,:), A, model, X, Y); % Function to sample

        Um = A_prior.normals(cs,3);
        US = A_prior.normals(cs,4);
%         support = [Um-3*US, Um, Um+3*US]; % Support points
%         A_entry_sample = IA2RMS(func, support, 10, 0);
%         A_entry_sample = A_entry_sample(10);

        % Do MH sampling
        A_entry_prev = model.A(A_prior.normals(cs,1),A_prior.normals(cs,2));
        A_entry_proposal = normrnd(A_entry_prev,US);
        prob_prev = func(A_entry_prev);
        prob_prop = func(A_entry_proposal);
        accept_prob = exp(min(0,prob_prop - prob_prev));
        if rand() < accept_prob
            A_entry_sample = A_entry_proposal;
        else
            A_entry_sample = A_entry_prev;
        end
        A_sampled(A_prior.normals(cs,1),A_prior.normals(cs,2)) = true;
        A(A_prior.normals(cs,1),A_prior.normals(cs,2)) = A_entry_sample;
    end
end
end

function [P] = A_func(A_entry, A_normal_prior, A, model, X, Y)
% IA2RMS likes to call this function on multiple A_entry values at once,
% so...

A_entry_row = A_normal_prior(1);
A_entry_col = A_normal_prior(2);
A_entry_m = A_normal_prior(3);
A_entry_S = A_normal_prior(4);

P = zeros(size(A_entry));

for i = 1:length(A_entry)
    A(A_entry_row,A_entry_col) = A_entry(i);

    Rinv = pinv(model.R);
    
    diff = Y - A*X;
    diff(isnan(diff)) = 0; % For the missing valued Y, set diff to 0
    log_cond = -0.5*trace(diff'*Rinv*diff);

    % Multiply by the prior; i.e. add log prior
    P(i) = log_cond + log(normpdf(A_entry(i),A_entry_m,A_entry_S));
end

end