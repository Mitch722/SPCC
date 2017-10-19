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
    
    [ no_violate, output_points, viol_fact ] = violation_function( Rad, cen, global_lite, x);
    
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
title('Frequency plot of Epsilon')
ylabel('Frequency of Violation Factor')
xlabel('Violation Factor')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdf = cumsum(freq,2);

pdf = pdf./max(pdf);

x_str = 'Violation Factor, Bin Size =';
bin_str = num2str(bound);
x_str = strcat(x_str,{' '}, bin_str);

titlestr = 'Cumulative Probability Distribution, n, M, runs,';
nstr = num2str(x.n);
Mstr = num2str(x.M);
setstr = num2str(x.no_sub_samp);

titlestr = strcat(titlestr,{' '},nstr,{','}, Mstr, {','}, setstr); 

figure 
plot(epsilon, pdf)
grid on
title(titlestr)
ylabel('Probability')
xlabel(x_str)

zeta = [2, x.dim+1];

% Theoretical Distributions 
% Theoretical distribution for when there are only limiting solutions for x.dim points 
y_xdim0 = 1 - binocdf(zeta(1) - 1, x.n, epsilon);

% Theoretical distribution for when there are only limiting solutions for x.dim points 

y_xdim1 = 1 - binocdf(zeta(2) - 1, x.n, epsilon);

hold on 
plot(epsilon, y_xdim0)

plot(epsilon, y_xdim1)

%% Unviolated Points Q 
% Find the frequency and PDF by binning the number of unviolated in the Multisample

% q: the number of unviolated points in x.M the global sample
Q = x.M .*(1 - cell2mat(Output_data(violation_factors_entry, :) ));

% Length of the array of bins
len_arrBound = 1000;

% bound: the bin size
bound = ( max(Q) - min(Q) ) / len_arrBound;

% gives the Frequency of Q
[freqQ, cumFreqQ, arrBQ, indices] = bin_var(Q, bound);

cumQ_norm = cumFreqQ ./ max(cumFreqQ);

arr_norm = arrBQ / x.M;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot frequency of q
figure
hold off
plot(arrBQ, freqQ, '.')
grid on

% poly_freqQ = freqQ;
% logic_poly = poly_freqQ == 0;
% poly_freqQ( logic_poly ) = [];
% 
% poly_arrBQ( logic_poly ) = [];
% 
% coeffsP = polyfit(poly_arrBQ, poly_freqQ, 2);
% 
% poly_arrBQ = [ poly_arrBQ .^2; poly_arrBQ ];
% 
% hold off
% figure
% plot( coeffs'*poly_arrBQ, poly_arrBQ(2,:) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold off
plot(arr_norm, cumQ_norm)
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strXlab = 'Number of Unviolated Points, bound size: ';
strBound = num2str(bound);
strXlab = strcat(strXlab,{' '},strBound);

title('Cumulative Distribution compared with Violation Probability')
ylabel('Cumulative Probability')
xlabel(strXlab)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zeta again
zeta = [2, x.dim+1];
% epsilon is 1 - q/x.M
% q/x.M = epsilon, epsilon = arr_norm
theo_3 = 1 - binocdf(Q - zeta(1), x.M, Q);
theo_4 = 1 - binocdf(3 - 1, x.n, arr_norm);

hold on
plot(arr_norm, theo_3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This bins all the data that is betweent the bound in Q
output_bin = cell(2, length(arrBQ));

for i = 1 : length(arrBQ)
    
    output_bin{1, i} = Q( indices(i, :) );
    output_bin{2, i} = arrBQ(i);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Save current workspace variables

save('out_put_data.mat','-v7.3')

