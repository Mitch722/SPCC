% Violation Probability 
%% Load Data or Generate it

try load('out_put_data_5.mat')

catch 
    disp('No output from previous run')
    
    [Radii, centres, x] = multisample_5;
    
    save('out_put_data_5.mat', '-v7.3')
    
end

fprintf('The multisample contains %d and each subset is %d \n', x.M, x.n )

% Output data entries:
% sub_samples_entry = 1;       sub samples are stored here 

%% Generate infinite distribution

y.M = 100000;
% x.dim, the dimension of the samples
y.dim = x.dim;
% x.sample, the global mutlisample
y.sample = randn(y.M, y.dim);


%% Generate Q
zeta = [2, x.dim+1];

no_violate = zeros(size(Radii'));
viol_fact = zeros(size(Radii'));

fprintf('Progress:')
parfor i  = 1: length(Radii)
    
    [ no_violate(i), ~, viol_fact(i) ] = violation_function( Radii(i), centres(i), y.sample, y);
    if (i / 100) == round(i / 100)
        fprintf('#')
    end
end

fprintf('\n')
%%
% Q the number of points which violate
Q = y.M * ( 1 - viol_fact);
% Length of the array of bins
len_arrBound = 100;

% bound: the bin size
bound = ( max(Q) - min(Q) ) / len_arrBound;
% make bound an integer number 
bound = round(bound);

% gives the Frequency of Q
[freqQ, cumFreqQ, arrBQ, indices] = bin_var(Q, bound);

% normalizes the cumlatiave frequency of Q
cumQ_norm = cumFreqQ ./ max(cumFreqQ);
% normalizes the array of Q points
arr_norm = arrBQ / y.M;

%% Theo34:
Theo34 = zeros(size(arrBQ));
fprintf('Progress:')
% Work using betaln
for i = 1 : length(arrBQ)
   
    q_fun = arrBQ( i );
    
    if q_fun > y.M
        q_fun = y.M;
    end
    
    a = q_fun - 2;
    
    b = y.M;
    % c is equivalent to 1 -  epsilon
    c = q_fun / y.M;
    
    Theo34(i) = bigbino(a, b, c);
    % run an assertion error in Theo34 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % assert(isnotnan(Theo34(i)), 'Assertion failed at run i=%d', i)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    fprintf('#')
        
end

fprintf('\n')

Theo34 = 1 - Theo34/max(Theo34);

% Theo35:
Theo35 = zeros(size(arrBQ));
% Work using betaln
fprintf('Progress:')
for i = 1 : length(arrBQ)
   
    q_fun = arrBQ( i );
    
    if q_fun > y.M
        q_fun = y.M;
    end
    
    a = q_fun - (x.dim + 2);
    
    b = y.M;
    % c is equivalent to 1 -  epsilon
    c = q_fun ./ y.M;
    
    Theo35(i) = bigbino( a,b,c );
   
    fprintf('#')
end

fprintf('\n')

Theo35 = Theo35/max(Theo35);

% Epsilon:

cumEps = cumFreqQ./max(cumFreqQ);

cumEps = 1 - cumEps;

epsilon = (1 - arr_norm);

% Plots
figure
% Plot theorem 3.4
plot(epsilon, Theo34, 'r')
grid on
hold on
% plot worked out data
plot( epsilon, cumEps, 'b')
% Plot theorem 3.5
plot(epsilon, Theo35, 'm')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% String Manipulation
strXlab = 'Epsilon, bound size: ';
strBound = num2str(bound/x.M);
strXlab = strcat(strXlab,{' '},strBound);

title('Probability of Violation is less than Epsilon')
ylabel('Cumulative Probability')
xlabel(strXlab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

