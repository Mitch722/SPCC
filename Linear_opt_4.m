%% Linear Optimisation with Chance Constraints
% Moving plane problem

tic
%% Gaussian distribution Samples and Subsamples

% generate a sample of N constraints fromt this m samples will be selected

N = 1000000;             % number of global samples

P = [0 , 0]';           % defines the mean Point for the Gaussian

D = randn(2,N) + P;         % generates 2 by N Gaussian distributed number with mean at P
D = abs(D);             % keeps only positive numbers


m = 1000;                         % how many samples in a subset    

if m >= N
    error('m is the size of the subset of N: it must be small than N')
end


num_subset = N / m;               % the number of subsets

if num_subset ~= abs(num_subset)
    error('m1 must be a non-negative integer')
end


%% Define plane

c = [1 1]';                         % the normal that defines the plane

d = c'*D;                           % distances from origin to plane

distances = vec2mat(d,m);            % splits the matrix up into subsets with m columns and num_subset rows

[max_value,max_Index] = max(distances,[],2);                % finds maximum value and relevent index   


Plane_grad1 = -c(1) / c(2);     % finds the gradient of line of the plane

%% Violation: find the probability that 1% of the global violate the subset sample

logic_big_than = max_value < d;

one_col = ones(length(logic_big_than),1);

num_violates = logic_big_than * one_col;


violation_factors = num_violates ./ N;

if max(violation_factors) > 1
    error('Too few samples: more points violate the subset than are in the subset. Increase m.')
end

%% Compare each violation factor in order to work out the probability that violation occurs

eta = linspace(0,0.1,length(violation_factors));

logic_big_eta = violation_factors > eta;


col_len =  size( logic_big_eta );
col_len = col_len(1);


one_col_2 = ones(col_len,1);

num_big_eta = logic_big_eta'*one_col_2;

probab_Violate = num_big_eta ./ m;

plot(eta,probab_Violate,'.')

grid on
toc        