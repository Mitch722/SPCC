% Make a Global Multi Sample of M points 
% Take n samples from it. n < M
% Run SOCP to get Radius and centre
% Store above analyse how many and which points violate model produced by
% the sub-sample
%% Generate the Multisample
% x.M, no. of global samples in the Multi-sample
x.M = 10000;
% x.no_sub_samp, the no of sets of sub-samples
x.no_sub_samp = 100;
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

tic
% Runs the for loop in parallel
for i = 1 : x.no_sub_samp
    % make a cell array for the subsamples
    sub_samp_mat = x.sample( (1 + (i - 1)*x.n ) : i*x.n , : );
    Output_data{1, i} = sub_samp_mat;
    
    [ Rad, cen, no ] = min_R_SOCP(sub_samp_mat);
    
    Radii(i, :) = Rad;
    centres(i, :) = cen;
    no_outside(i, :) = no;
    
    global_lite = x.sample;
    global_lite( (1 + (i - 1)*x.n ) : i*x.n , : ) = [];
    
    [ no_violate, output_points, viol_fact ] = violation_function( R, c, global_lite );
    
    Output_data{2, i} = no_violate;
    Output_data{3, i} = output_points;
    Output_data{4, i} = viol_fact;
end
toc    
        
    
    
