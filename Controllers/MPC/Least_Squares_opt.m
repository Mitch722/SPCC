%% Optimisation step
% least squares a go

function [epsilon, probab] = Least_Squares_opt(x, len_sub, q_hat)

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

x_opt = x.sample(:, end - len_sub + 1: end);

coef = polyfit(t, x_opt, n); 

c = [coef(1, 1), 1]';

% Take last sub-samples worth
x_opt = x_opt - coef(1, 2);

x_opt = [t; x_opt];

phi = atan(c(1));
% rotation matrix
% R = [c -s; c  s];
% phi = 0;
R = [cos(phi), -sin(phi);    
    sin(phi),   cos(phi)];

x_R = R * x_opt;

phi_dist = phi;
max_dist = max(x_R(2, :));

%% Compare Model to x.sample and work out violation probability

x.sample(:, end - len_sub: end) = [];

violations = x.sample.^2 <= max_dist^2;

q_eps = sum(violations, 2);
epsilon = 1 -  sum(violations, 2) / length(x.sample);

%% Work out if epsilon fits within the 95% interval 
% Posterior distribution

eps_k = 1 - q_hat / (length(x.sample) + len_sub);
probab = binopdf(q_eps - 2, length(x.sample), 1 - eps_k);



