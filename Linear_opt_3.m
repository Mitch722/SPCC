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

distances = zeros(num_subset,m);            % distances from the plane 


for i = 1 : num_subset
    
    distances(i, (1 : m)) = d( (i - 1)*m + 1 : i*m  );      % splits up the distances into 10 rows
    
end

[max_value,max_Index] = max(distances,[],2);                % finds maximum value and relevent index   


Plane_grad1 = -c(1) / c(2);     % finds the gradient of line of the plane

%% Violation: find the probability that 1% of the global violate the subset sample

num_violate = zeros( num_subset,1 );

% Loop to test each maximum value in order to see how many global points
% violate the subset 
for i = 1 : num_subset
    
    counter = 0;
    
    maxI = max_value(i,1);
    length_d = length(d);
    
    for j = 1 : length_d
    
        if d(1,j) > maxI
            counter = counter +1;
        end
    end
    num_violate(i,1) = counter;
    
end

% Need to divide all NUMBER OF TIMES that the condition is VIOLATED  in
% order to find the violation factor this can be compared to eta

violation_factors = num_violate ./ N;


if max(violation_factors) > 1
    error('Too few samples: more points violate the subset than are in the subset. Increase m.')
end

%% Compare each violation factor in order to work out the probability that violation occurs

eta = linspace(0,0.1,m);
num_bigger_eta = zeros(length(eta) ,1);

for i = 1 : length(eta)
    
    eta_comp = eta(1,i);
    counter = 0;
    
    for j = 1 : num_subset
        
        viol_comp = violation_factors(j,1);
        
        if viol_comp > eta_comp
            
            counter = counter + 1;
        end 
        num_bigger_eta(i,1) = counter;   
    end
    
end

one_minus_eta = 1 - eta;
probab_Violate = num_bigger_eta ./ m;

probab_V = 1 - probab_Violate;

plot(eta,probab_V,'.')   

grid on
toc        