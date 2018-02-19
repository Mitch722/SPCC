try 
    try_statement = A;
    
catch
    MPC_adaptive_1
    
end

params.m = 50;
params.n = 100;
params.Ts = 0.01;

% p must be divisible by three
p = 6;

Q_bar = eye(p); 

Y = y(:, end - 300: end);
Ck1 = Ck(:, end - 300: end); 

bnds = [0.8, 0.15, 0.4]';

c = adaptiveControl(Q_bar, Y, Ck1, p, params, bnds);
