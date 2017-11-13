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
tic
y.M = 1000000;
% x.dim, the dimension of the samples
y.dim = x.dim;
% x.sample, the global mutlisample
y.sample = randn(y.M, y.dim);


%% Generate Q
zeta = [2, x.dim+1];

no_violate = zeros(size(Radii'));
viol_fact = zeros(size(Radii'));

samp = y.sample;

fprintf('Progress:')
parfor i  = 1: length(Radii)
    
    [ no_violate(i), ~, viol_fact(i) ] = violation_function( Radii(i), centres(i), samp, y);
    if (i / 100) == round(i / 100)
        fprintf('#')
    end
end

fprintf('\n')

%% number of points in the 

bound = 0.0001;
% this data is theoretical is compared to an 'infinte' sized distributionl
[freq, cumFreq, epsilon, indices_trans] = bin_var(viol_fact, bound);


%% Theo34:
ep_Theo = linspace(min(epsilon), max(epsilon), 100);

Theo34 = zeros(size(ep_Theo));
fprintf('Progress:')
% Work using betaln
for i = 1 : length(ep_Theo)
   
    eps_val = ep_Theo(i);
    
    a = x.n - zeta(1);
    
    b = x.n;
    % c is equivalent to 1 -  epsilon
    c = 1 - eps_val;
    
    Theo34(i) = binocdf(a, b, c);
    % run an assertion error in Theo34 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % assert(isnotnan(Theo34(i)), 'Assertion failed at run i=%d', i)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if i / 500 == round(i / 500)
       fprintf('#')
    end
        
end

fprintf('\n')

Theo35 = zeros(size(ep_Theo));
fprintf('Progress:')
% Work using betaln
for i = 1 : length(ep_Theo)
   
    eps_val = ep_Theo(i);
    
    a = x.n - zeta(2);
    
    b = x.n;
    % c is equivalent to 1 -  epsilon
    c = 1 - eps_val;
    
    Theo35(i) = binocdf(a, b, c);
    % run an assertion error in Theo34 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % assert(isnotnan(Theo34(i)), 'Assertion failed at run i=%d', i)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if i / 500 == round(i / 500)
       fprintf('#')
    end
        
end

fprintf('\n')

% Plots

figure

plot(epsilon, cumFreq./max(cumFreq), 'b')
hold on
plot(ep_Theo, Theo34/max(Theo34), 'm')
hold on
plot(ep_Theo, Theo35/max(Theo35), 'r')
grid on

% String Manipulation
strXlab = 'Epsilon, bound size: ';
strBound = num2str(bound/x.M);
strXlab = strcat(strXlab,{' '},strBound);

title('Probability of Violation is less than Epsilon')
ylabel('Cumulative Probability')
xlabel(strXlab)

hold on
plot(1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
