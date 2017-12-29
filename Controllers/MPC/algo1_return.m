function [x_hat, q_hat] = algo1_return(x, Ntrail, q_min, q_max)

x.M = length(x.sample);
no_sub_samps = x.M / Ntrail;

% Take a round number which is smaller than 
if no_sub_samps >= round(no_sub_samps)
    no_sub_samps = round(no_sub_samps);
else
   no_sub_samps = round(no_sub_samps) - 1; 
end

% find the remainder and remove these from x
r = rem(x.M, no_sub_samps);

if r ~= 0
   x.sample(:, end-r+1:end) = []; 
end

x_data = x.sample;
x_reshape = reshape(x_data, [], no_sub_samps);

%% Optimisation step
% least squares a go
[~, len_samples] = size(x_reshape);

t = 1: len_samples;
n = 1;

coef = zeros(no_sub_samps, n+1);
max_dist = zeros(no_sub_samps, 1);
phi_dist = max_dist;

for i = 1: no_sub_samps
       
   coef_inter = polyfit(t, x_reshape(i, :), n); 
    
   coef(i, :) = coef_inter;
   % y_model = coef_inter*[t; ones(1, length(t))];
   
   c = [coef_inter(1, 1), 1]';
   
   x_opt = x_reshape(i, :);
   x_opt = x_opt - coef_inter(1, 2);
   
   x_opt = [t; x_opt];
   x_opt = x_opt;
   
   phi = -atan(c(1));
   % rotation matrix
   % R = [c -s; c  s];
   R = [cos(phi), -sin(phi);    
        sin(phi),   cos(phi)];
    
   x_R = R * x_opt;
   
   phi_dist(i, 1) = phi;
   max_dist(i, 1) = max(norm(x_R(2, :)));
   
end



