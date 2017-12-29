%% test 

M = 100;

x = struct;
x.sample = randn(1, M);

[x_hat, q_hat] = algo1_return(x, 6, 0, 0);

