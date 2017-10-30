%% Generate the Multisample
% x.M, no. of global samples in the Multi-sample
x.M = 600;
% x.no_sub_samp, the no of sets of sub-samples
x.no_sub_samp = 1;
% x.n, no. of samples in the sub-sample
x.n = x.M / x.no_sub_samp;

assert(x.n == norm(x.n), 'The number of Samples is not divisable by the number of sub-samples')

% x.dim, the dimension of the samples
x.dim = 2;
% x.sample, the global mutlisample
x.sample = randn(x.M, x.dim);

save('findR.mat', '-v7.3')

%%
zeta = [2, x.dim + 1];
n = 1;

zeta = zeta(n);

fin = x.M;

r = zeta : x.M;

rhs = 1 ./ (x.M - r + 1);

lhs = psi( fin - r + 1 ) - psi( fin - r + 2) - psi(r - zeta + 1 ) + psi( r + 1 );

hs = rhs + lhs;

plot(r, hs)
grid on
