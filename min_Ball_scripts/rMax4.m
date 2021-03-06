%% Generate a plot to see what r looks like

% x.M, no. of global samples in the Multi-sample
x.M = 1000;
% x.no_sub_samp, the no of sets of sub-samples
x.no_sub_samp = 1;
% x.n, no. of samples in the sub-sample
x.n = x.M / x.no_sub_samp;

assert(x.n == norm(x.n), 'The number of Samples is not divisable by the number of sub-samples')

% x.dim, the dimension of the samples
x.dim = 2;
% x.sample, the global mutlisample
x.sample = randn(x.M, x.dim);

save('findR2.mat', '-v7.3')
    

%%
% produce r
r = 10;

r_1 = r;

q = r : x.M;

zeta = [2, x.dim+1];
n2 = 1;
zed = zeta(n2);

n = x.M - r;
k = q - r;

beta_coefln =  betaln(n - k + 1, k+1);
beta_coefln = log(1 ./ (n+1)) - beta_coefln  ;

nom = betaln(x.M - q + zed, q - zed + 1);
denom = betaln(zed, r - zed + 1);

total = beta_coefln + nom - denom;

func = exp(total);

figure
plot(q, func)
