function [X, P, X_exp_prev, P_exp_prev, K, eps, Sigma_eps, loglike] = ssm_kalman_filter(model, Y, U)
% Calculate most likely $X_i$ conditioned on $Y_i$ from model $\Theta = 
% \{\mu_0, \Sigma_0, \Phi, A, Q, R, \Upsilon\}$, data $Y$, and input
% $U$. Returns $X_n^n$ X, variances $P_n^n$ P, $X_n^{n-1}$ X_exp_prev,
% variances $P_n^{n-1}$ P_exp_prev, Kalman gain K, innovations eps,
% innovations covariance Sigma_eps, and the negative log likelihood loglike

Tn = size(Y,2);
Yn = size(Y,1);
Xn = size(model.mu_0,1);
Un = size(U,1);

% Do some dimension checking
if(size(U,2) ~= Tn)
    throw(MException('SSM:ssm_kalman_filter:WrongDimensions','U does not have the same number of time points as Y.'));
end
if(~all(size(model.Sigma_0) == [Xn Xn]))
    throw(MException('SSM:ssm_kalman_filter:WrongDimensions','model.Sigma_0 is the wrong size.'));
end
if(~all(size(model.Phi) == [Xn Xn]))
    throw(MException('SSM:ssm_kalman_filter:WrongDimensions','model.Phi is the wrong size.'));
end
if(size(model.A,1) ~= Yn || size(model.A,2) ~= Xn)
    throw(MException('SSM:ssm_kalman_filter:WrongDimensions','model.A is the wrong size.'));
end
if(~all(size(model.Q) == [Xn Xn]))
    throw(MException('SSM:ssm_kalman_filter:WrongDimensions','model.Q is the wrong size.'));
end
if(~all(size(model.R) == [Yn Yn]))
    throw(MException('SSM:ssm_kalman_filter:WrongDimensions','model.R is the wrong size.'));
end
if(size(model.Upsilon,1) ~= Xn || size(model.Upsilon,2) ~= Un)
    throw(MException('SSM:ssm_kalman_filter:WrongDimensions','model.Upsilon is the wrong size.'));
end

% Initialize

X = zeros(Xn,Tn);
P = zeros(Xn,Xn,Tn);
X_exp_prev = zeros(size(X));
P_exp_prev = zeros(size(P));
K = zeros(Xn,Yn,Tn);
eps = zeros(size(Y));
Sigma_eps = zeros(Yn,Yn,Tn);
loglike = 0;

% From Shumway and Stoffer TSA4 Chapter 6 State Space Models Property 6.1
for i = 1:Tn
    % Make substitutions for missing data if necessary
    Y_i = Y(:,i);
    ind_missing = isnan(Y_i);
    Y_i(ind_missing) = 0;
    A_i = model.A;
    A_i(ind_missing,:) = 0;
    R_i = model.R;
    R_i(ind_missing,:) = 0;
    R_i(:,ind_missing) = 0;
    R_i(sub2ind(size(R_i),find(ind_missing),find(ind_missing))) = 1;
    
    % x_prev and P_prev are $x_{i-1}$ and $P_{i-1}$
    if(i == 1)
        x_prev = model.mu_0;
        P_prev = model.Sigma_0;
    else
        x_prev = X(:,i-1);
        P_prev = P(:,:,i-1);
    end
    
    % x_exp_prev_i and P_exp_prev_i are the expected values of $x_i$ and $P_i$
    % conditional on $Y_{i-1}$
    x_exp_prev_i = model.Phi*x_prev + model.Upsilon*U(:,i);
    P_exp_prev_i = model.Phi*P_prev*(model.Phi') + model.Q;
    
    % Calculate innovations eps
    eps_i = Y_i - A_i*x_exp_prev_i;
    
    % Calculate Kalman gain K
    % In the situation in which both P_exp_prev_i and R_i are 0, the Kalman
    % gain converges to 0. We check to see if the denominator is invertible
    K_denom = (A_i*P_exp_prev_i*(A_i')+R_i);
    if(det(K_denom) == 0)
        K_i = zeros(Xn,Yn);
    else
        K_i = P_exp_prev_i*(A_i')/K_denom;
    end
    % Calculate current x and P
    % x_curr and P_curr are $x_i$ and $P_i$
    x_i = x_exp_prev_i + K_i*(Y_i - A_i*x_exp_prev_i);
    P_i = (eye(Xn) - K_i*A_i)*P_exp_prev_i;
    
    Sigma_eps_i = A_i*P_exp_prev_i*A_i' + R_i;
    % Make sure Sigma_eps_curr is symmetric
    Sigma_eps_i = (Sigma_eps_i + Sigma_eps_i')/2;
    
    X(:,i) = x_i;
    P(:,:,i) = correct_cov(P_i);
    X_exp_prev(:,i) = x_exp_prev_i;
    P_exp_prev(:,:,i) = correct_cov(P_exp_prev_i);
    K(:,:,i) = K_i;
    eps(:,i) = eps_i;
    Sigma_eps(:,:,i) = Sigma_eps_i;
    
    if(all(Sigma_eps_i(:) == 0))
        loglike = loglike + -Inf;
    else
        loglike = loglike + log(det(Sigma_eps_i)) + (eps_i'/Sigma_eps_i)*eps_i;
    end
end

loglike = -loglike/2;