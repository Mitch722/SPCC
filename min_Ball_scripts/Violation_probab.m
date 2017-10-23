% Violation Probability 
%% Load Data or Generate it

try load('out_put_data_2.mat')

catch 
    disp('No output from previous run')
    
    multisample_2
    
end
% Output data entries:
% sub_samples_entry = 1;       sub samples are stored here 
% no_violations_entry = 2;     the number that violate model 
% violation_points_entry = 3;  points that violate
% violation_factors_entry = 4; epsilon values
%% Generate Q
zeta = [2, x.dim+1];
% Q the number of points which violate
Q = x.M .*(1 - cell2mat(Output_data(violation_factors_entry, :) ));
% Length of the array of bins
len_arrBound = 50;

% bound: the bin size
bound = ( max(Q) - min(Q) ) / len_arrBound;
% make bound an integer number 
bound = round(bound);

% gives the Frequency of Q
[freqQ, cumFreqQ, arrBQ, indices] = bin_var(Q, bound);

% normalizes the cumlatiave frequency of Q
cumQ_norm = cumFreqQ ./ max(cumFreqQ);
% normalizes the array of Q points
arr_norm = arrBQ / x.M;

%% Theo34:
Theo34 = zeros(size(arrBQ));
% Work using betaln
for i = 1 : length(arrBQ)
   
    q_fun = arrBQ( i );
    a = q_fun - 90;
    
    b = x.M;
    % c is equivalent to 1 -  epsilon
    c = q_fun ./ x.M;
    
    Theo34(i) = binopdf( a,b,c );
    
end
Theo34 = Theo34 ./ max(Theo34);

%% Theo35:
Theo35 = zeros(size(arrBQ));
% Work using betaln
for i = 1 : length(arrBQ)
   
    q_fun = arrBQ( i );
    a = q_fun - 98;
    
    b = x.M;
    % c is equivalent to 1 -  epsilon
    c = q_fun ./ x.M;
    
    Theo35(i) = binopdf( a,b,c );
end

Theo35 = Theo35 / max(Theo35);

%% Epsilon:

cumEps = cumFreqQ./max(cumFreqQ);

cumEps = 1 - cumEps;

epsilon = (1 - arr_norm);

%% Plots
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

