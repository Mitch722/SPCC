% Make a Global Multi Sample of M points 
% Take n samples from it. n < M
% Run SOCP to get Radius and centre
% Store above analyse how many and which points violate model produced by
% the sub-sample
%% Generate the Multisample
% x.M, no. of global samples in the Multi-sample
x.M = 100000;
% x.no_sub_samp, the no of sets of sub-samples
x.no_sub_samp = 1000;
% x.n, no. of samples in the sub-sample
x.n = x.M / x.no_sub_samp;
% x.dim, the dimension of the samples
x.dim = 2;
% x.sample, the global mutlisample
x.sample = randn(x.M, x.dim);

%% Begin the analysis
% Preallocate space
Radii = zeros(x.no_sub_samp, 1);
centres = zeros(x.no_sub_samp, x.dim);
no_outside = zeros(x.no_sub_samp, 1);

% cell array containing the sub samples
Output_data{4, x.no_sub_samp} = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output data is stored in the rows of a cell array
sub_samples_entry = 1;
no_violations_entry = 2;
violation_points_entry = 3;
violation_factors_entry = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run through the sub sample data

tic
% Runs the for loop in parallel
for i = 1 : x.no_sub_samp
    % make a cell array for the subsamples
    sub_samp_mat = x.sample( (1 + (i - 1)*x.n ) : i*x.n , : );
    Output_data{sub_samples_entry, i} = sub_samp_mat;
    
    [ Rad, cen, no ] = min_R_SOCP(sub_samp_mat);
    
    Radii(i, :) = Rad;
    centres(i, :) = cen;
    no_outside(i, :) = no;
    
    global_lite = x.sample;
    global_lite( (1 + (i - 1)*x.n ) : i*x.n , : ) = [];
    
    [ no_violate, output_points, viol_fact ] = violation_function( Rad, cen, global_lite );
    
    Output_data{no_violations_entry, i} = no_violate;
    Output_data{violation_points_entry, i} = output_points;
    Output_data{violation_factors_entry, i} = viol_fact;
end
toc    

%% plot epsilon against 

epsilon = linspace(0,0.1,x.no_sub_samp);

e_data = Output_data(violation_factors_entry,:);
e_data = cell2mat(e_data);

% find how many points are within 0.1 of each value of epsilon
% ep_repeated = repmat(epsilon, x.no_sub_samp, 1);

diff = e_data' - epsilon;
diff = abs(diff);

bound = (epsilon(end) - epsilon(1)) / length(epsilon);

logic_number = diff < bound;

freq = sum(logic_number, 1);

figure
plot(epsilon, freq, '+')
grid on

pdf = cumsum(freq,2);

pdf = pdf./max(pdf);

figure 
plot(epsilon, pdf)
grid on