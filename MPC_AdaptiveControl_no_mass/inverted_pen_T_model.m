function sysd = inverted_pen_T_model(Ts, M, m)

%% Make the Inverted Pendulum model
% Physical data stored in a struct
% m = 0.2
data.m = m;
% M = 1.5
data.M = M;
data.I = 0.005;
data.l = 0.5;
data.b = 0.1;
data.g = 9.81;

[~, sysd] = pendulum_ss(data, Ts);

end