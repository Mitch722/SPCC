%% test 

M = 1000;

x = struct;
x.sample = randn(1, M);
x.dim = 2;

ep_lo = 0.05;
ep_hi = 0.10;

Ppor = 0.9;
Ppst = 0.95;

%% Run Algorithm 1

[r_star, Ptrail, Ntrail, q_min, q_max] = algo1(x, ep_lo, ep_hi, Ppor, Ppst);

[x_hat, q_hat] = algo1_return(x, Ntrail, q_min, q_max);

%%
figure
plot(x.sample,'.')
hold on
plot([1, M], [x_hat(2, 1), x_hat(2, 1)])
hold on
plot([1, M], [-x_hat(2, 1), -x_hat(2, 1)])
