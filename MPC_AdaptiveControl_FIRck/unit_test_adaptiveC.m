load('output_input.mat')

params.m = 6;
params.n = 200;
params.Ts = 0.01;

% p must be divisible by three
p = 6;

Q_bar = eye(p); 

y(:, end) = [];
Y = y(:, 1: 300);
Ck1 = Ck(:, end - 300: end); 

bnds = [0.8, 0.15, 0.4]';

c = adaptiveControl2(Q_bar, Y, Ck1, p, params, bnds);
