%% Linear Optimisation with Chance Constraints
% Moving plane problem


%% Gaussian distribution Samples and Subsamples
% generate a sample of N constraints fromt this m samples will be selected

N = 1000;  % number of global samples

P = [0 , 0]';           % defines the mean Point for the Gaussian

D = randn(2,N) + P;         % generates 2 by N Gaussian distributed number with mean at P
D = abs(D);             % keeps only positive numbers

m = 100;                % subset of samples
D1 = D(:,1:m);        % first subset 

%% Define plane

c = [1 1]';

d1 = c'*D1;

dmax1 = max(d1);

[~,I1] = max(d1);
x_d1 = D1(:,I1);

Plane_grad1 = -c(1) / c(2);

plane_Y_intercept = x_d1(2) - Plane_grad1*x_d1(1);
plane_X_intercept = -plane_Y_intercept / Plane_grad1;

plane = [0 plane_X_intercept  ; plane_Y_intercept  0];

%% Plots

figure
axis square
plot(D(1,:),D(2,:),'.')         % plots Largest set of random samples

figure
axis square

plot(D1(1,:),D1(2,:),'.')         % plots Largest set of random samples
hold on 
plot(plane(1,:),plane(2,:),'r')

hold off





