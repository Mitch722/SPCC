% Try to load previous workspace
% if no previously stored workspace is avaliable then make and save a new
% one

try load('out_put_data_1.mat')

catch 
    disp('No output from previous run')
    
    multisample_1
    
end

fprintf('The multisample contains %d and each subset is %d \n', x.M, x.n )

% Output data entries:
% sub_samples_entry = 1;       sub samples are stored here 
% no_violations_entry = 2;     the number that violate model 
% violation_points_entry = 3;  points that violate
% violation_factors_entry = 4; epsilon values

%% Unviolated Points Q 
% Find the frequency and PDF by binning the number of unviolated in the Multisample

% q: the number of unviolated points in x.M the global sample
Q = x.M .*(1 - cell2mat(Output_data(violation_factors_entry, :) ));
% find the size of each essential set r
R = x.M - cell2mat( Output_data(no_violations_entry, :) );
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
%% String Manipulation for plots
strXlab = 'Number of Unviolated Points, bound size: ';
strBound = num2str(bound);
strXlab = strcat(strXlab,{' '},strBound);

title('Cumulative Distribution compared with Violation Probability')
ylabel('Cumulative Probability')
xlabel(strXlab)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Laplacian to find PDF from CDF
% first smooth with Gaussian
% pad cdf losing info
pad = 100;
cumFreqQ2 = [cumQ_norm, ones(1, pad)];
% remove added values less 1 for derivative
sdtdev = 2;

smooth_cdf = smoothts(cumFreqQ2, 'g', 5, sdtdev);
smooth_cdf(end - pad +1: end) = [];

figure 
plot(arrBQ, smooth_cdf)
grid on

hold off
% applies derivative
pdf = diff( smooth_cdf ) ./ diff(arrBQ);

arrBQ2 = arrBQ;
arrBQ2(end) = [];
arrBQ2 = arrBQ2 + diff(arrBQ)/2;

areapdf = sum( smooth_cdf, 2 ) * bound;

plot(arrBQ2/x.M, pdf./max(pdf),'.')
grid on
hold on
%% Fit a polynomial to the diff data
n = 2;

norm_arrBQ2 = arrBQ2/x.M;
norm_pdf = pdf./areapdf;
coef_pdf = polyfit(norm_arrBQ2, norm_pdf, n);

poly_pdf = coef_pdf * [ norm_arrBQ2.^2; norm_arrBQ2; ones(size(norm_arrBQ2)) ];

% plot(norm_arrBQ2, poly_pdf, 'b')

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
r = x.n;

n = x.M - r;

Theo36 = zeros(size(arrBQ2));

zed = zeta(1);
for i = 1: length(arrBQ2)
    
    % r =  arrBQ(i);
    
    k = arrBQ2(i) - r;
    beta_coefln =  betaln(n - k + 1, k+1);
    beta_coefln = log(1 / (n+1)) - beta_coefln  ;
    
    a = x.M - arrBQ2(i) + zed;
    b = arrBQ2(i) - zed + 1;
    
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
plot( arrBQ2/x.M, Theo36,'r')
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% r = 25;



Theo37 = zeros(size(arrBQ2));

zed = zeta(2);
for i = 1: length(arrBQ2)
    
    % r = arrBQ(i);
    
    n = x.M - r;
    
    k = arrBQ2(i) - r;
    beta_coefln =  betaln(n - k + 1, k+1);
    beta_coefln = log(1 / (n+1)) - beta_coefln  ;
    
    a = x.M - arrBQ2(i) + zed;
    b = arrBQ2(i) - zed + 1;
    
    c = zed;
    d = r - zed + 1;
    
    if d < 0
        d = 1;
    end
    
    bulk = betaln(a,b) - betaln(c,d);
    assert(bulk <= eps, 'bulk less than or equal to eps')
    
    Theo37_inter = beta_coefln + bulk;
    Theo37(i) = exp(Theo37_inter);
    
end

hold on
plot( arrBQ2/x.M, Theo37, 'm')
grid on

hold on
plot(arrBQ2/x.M, pdf./max(pdf) * max(Theo37))
grid on

hold off

tit_str = 'Pdf of Q with r =';
str_r = num2str(r);
tit_str = strcat(tit_str, str_r);

title(tit_str)
xlabel('Q/M')
ylabel('Probability')
%% Plots






%% Work out the a polynomial fit for the Cumulative plot
% n = 3;
% concA = x.M:bound:2000;
% arrBQ1 = [arrBQ, concA];
% cumQ_poly = [cumQ_norm, ones(size(concA))];
% 
% coef_CF = polyfit(arrBQ1, cumQ_poly, n);
% 
% arrQ = [arrBQ.^3 ; arrBQ.^2 ; arrBQ ; ones( size(arrBQ) )];
% 
% cumFreq_poly = coef_CF * arrQ;
% 
% hold on
% plot(arrBQ, cumFreq_poly)
% 
% %% Derivatives of Cumulative Plot: pdf
% 
% coef_PF = [coef_CF(1)*3, coef_CF(2)*2, coef_CF(3)];
% 
% arrQF = [arrBQ.^2 ; arrBQ ; ones( size(arrBQ) )];
% 
% cdfQ = coef_PF * arrQF;
% areaCDF = sum(cdfQ, 2);
% 
% figure 
% plot(arrBQ, cdfQ/max(cdfQ))
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

