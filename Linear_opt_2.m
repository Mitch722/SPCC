%% Linear Optimisation with Chance Constraints
% Moving plane problem


%% Gaussian distribution Samples and Subsamples

% generate a sample of N constraints fromt this m samples will be selected

N = 1000;  % number of global samples

P = [0 , 0]';           % defines the mean Point for the Gaussian

D = randn(2,N) + P;         % generates 2 by N Gaussian distributed number with mean at P
D = abs(D);             % keeps only positive numbers

m1 = 10;                % subset of samples

m = N / m1;         

if m >= N
    error('m is the size of the subset of N: it must be small than N')
end

%% Define plane

c = [1 1]';                         % the normal that defines the plane

d = c'*D;                           % distances from origin to plane

distances = zeros(m1,m);            % distances from the plane 

max = zeros(m1,1);
index = zeros(m1,1);

for i = 1 : m1
    
    distances(i, (1 : m)) = d( (i - 1)*m + 1 : i*m  );
    
end


Plane_grad1 = -c(1) / c(2);     % finds the gradient of line of the plane


