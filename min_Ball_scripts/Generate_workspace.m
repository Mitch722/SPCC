% Try to load previous workspace
% if no previously stored workspace is avaliable then make and save a new
% one

try load('out_put_data.mat')

catch 
    disp('No output from previous run')
    
    multisample
    save('out_put_data.mat')
end

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

theo_3 = zeros(size(Q));
for i = 1 : length(arrBQ)
    
    q_var = round(arrBQ(i)); 
    theo_3(1, i) = binocdf( q_var - zeta(1), x.M,  q_var./x.M );
    
end

theo_4 = 1 - binocdf(3 - 1, x.n, arr_norm);

figure
plot(arr_norm, theo_3 );

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

% save('out_put_data.mat','-v7.3')
