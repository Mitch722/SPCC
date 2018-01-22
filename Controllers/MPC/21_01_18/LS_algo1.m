%% Optimisation step
% least squares a go

function [epsilon, confidence] = LS_algo1(x, len_sub)

% 1. Take the innovation data stored as a struct in x.
% 2. Least Squares model it. Producing an upper bound.
% 3. Compare this model to x.sample
% 4. If epsilon fits outside 95% of the probability then use model from
% smaller subsample
% 4. i) Start sampling from k
% 5. Else, continue with sample

%% Take the last sub_sample and Least Squares it

% [~, len_samples] = size(x.sample);

t = 1: len_sub;
n = 1;

% coef = zeros(no_sub_samps, n+1);
% max_dist = zeros(1, no_sub_samps);
% phi_dist = max_dist;

% 
x_opt.sample = x.sample(:, end - len_sub + 1: end);
x_opt.dim = x.dim;

%% 
ep_lo = 0.05;
ep_hi = 0.10;

Ppor = 0.9;
Ppst = 0.95;

%%

[~, ~, Ntrail, q_min, q_max] = algo1(x_opt, ep_lo, ep_hi, Ppor, Ppst);

[x_hat, ~] = algo1_return(x_opt, Ntrail, q_min, q_max);

max_dist = x_hat(2);

%% Compare Model to x.sample and work out violation probability

% remove samples used to train model
x.sample(:, end - len_sub: end) = [];

violations = x.sample.^2 >= max_dist^2;

epsilon = sum(violations, 2) / length(x.sample);
q_eps = length(x.sample) - sum(violations, 2);

%% Work out if epsilon fits within the 95% interval 
% Posterior distribution

% eps_k = 1 - q_hat / (length(x.sample) + len_sub);
confidence = binocdf(q_eps - 2, length(x.sample), 1 - epsilon);

