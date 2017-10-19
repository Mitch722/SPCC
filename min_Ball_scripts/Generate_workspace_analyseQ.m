% Try to load previous workspace
% if no previously stored workspace is avaliable then make and save a new
% one

try load('out_put_data_1.mat')

catch 
    disp('No output from previous run')
    
    multisample_1
    
end

%% Unviolated Points Q 
% Find the frequency and PDF by binning the number of unviolated in the Multisample

% q: the number of unviolated points in x.M the global sample
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot frequency of q
figure
hold off
plot(arrBQ, freqQ/max(freqQ), '.')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold off
plot(arrBQ, cumQ_norm)
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%String Manipulation for plots
strXlab = 'Number of Unviolated Points, bound size: ';
strBound = num2str(bound);
strXlab = strcat(strXlab,{' '},strBound);

title('Cumulative Distribution compared with Violation Probability')
ylabel('Cumulative Probability')
xlabel(strXlab)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This bins all the data that is betweent the bound in Q
output_bin = cell(2, length(arrBQ));

for i = 1 : length(arrBQ)
    
    output_bin{1, i} = Q( indices(i, :) );
    output_bin{2, i} = arrBQ(i);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Work out Theoretical plots || Try Theorem 4:
% zeta again
zeta = [2, x.dim+1];
% epsilon is 1 - q/x.M
% q/x.M = epsilon, epsilon = arr_norm

% r parameter
r = 14;

n = x.M - r;

Theo36 = zeros(size(arrBQ));

zed = zeta(1);
for i = 1: length(arrBQ)
    
    k = arrBQ(i) - r;
    beta_coefln =  betaln(n - k + 1, k+1);
    beta_coefln = log(1 / (n+1)) - beta_coefln  ;
    
    a = x.M - arrBQ(i) + zed;
    b = arrBQ(i) - zed + 1;
    
    c = zed;
    d = r - zed + 1;
    
    if d < 0
        d = 1;
    end
    
    bulk = betaln(a,b) - betaln(c,d);
    assert(bulk <= eps, 'bulk less than or equal to eps')
    
    Theo36_inter = beta_coefln + bulk;
    Theo36(i) = exp(Theo36_inter);
    
end

figure
plot( arrBQ, Theo36./max(Theo36) ,'r')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logic_freqQ = freqQ == 0;
% compress freqQ_comp making it sparse
freqQ_comp = freqQ;
freqQ_comp(logic_freqQ) = [];
% 
arrBQ_comp = arrBQ;
% compress arrBQ to make it sparse
arrBQ_comp(logic_freqQ) = [];

n = 3;
coeffs = polyfit(arrBQ_comp, freqQ_comp, n);

freqQ_poly = [arrBQ.^3; arrBQ.^2 ; arrBQ; ones( size (arrBQ ) )];
freqQ_poly = coeffs * freqQ_poly;

hold on
plot(arrBQ, freqQ_poly/max(freqQ_poly),'b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
r = 3;

n = x.M - r;

Theo37 = zeros(size(arrBQ));

zed = zeta(2);
for i = 1: length(arrBQ)
    
    k = arrBQ(i) - r;
    beta_coefln =  betaln(n - k + 1, k+1);
    beta_coefln = log(1 / (n+1)) - beta_coefln  ;
    
    a = x.M - arrBQ(i) + zed;
    b = arrBQ(i) - zed + 1;
    
    c = zed;
    d = r - zed + 1;
    
    if d < 0
        d = 1;
    end
    
    bulk = betaln(a,b) - betaln(c,d);
    assert(bulk <= eps, 'bulk less than or equal to eps')
    
    Theo37_inter = beta_coefln + bulk;
    Theo37(i) = exp(Theo36_inter);
    
end

hold on
plot( arrBQ, Theo37./max(Theo37) ,'r')

