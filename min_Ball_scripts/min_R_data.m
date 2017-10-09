% This script runs min_R_SOCP_1 n times in order to gather data on how well
% the gurobi optimsation gets to the sample mean which should be (0,0)
%% Gather data

n = 1000;
% n: no. of samples
m = 1;
% m: number of repeats
dim_x = 2;
centres = zeros(m, dim_x);

for run = 1 : m
    
    x = randn(n, dim_x);
    
    [R, c, no_viol] = min_R_SOCP(x);
    centres(run, :) = c;
       
end

%% Manipulate data using MLE for a Gaussian

Mean_estimate = sum(centres) ./ n;
Variance_estimate = sum( (centres - Mean_estimate).^2 )./ n;

norm_data = normpdf(centres);

figure
plot(centres(:,1), norm_data(:,1),'.')
hold on
plot(centres(:,2), norm_data(:,2),'+')
grid on
