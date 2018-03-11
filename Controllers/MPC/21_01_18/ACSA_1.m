function [ck, Mmodels] = ACSA_1(M, y, u, p, params, Q_bar, bnds)


%% check that there is enough data
assert(length(y) >= M*(params.n+params.m), 'There is not enough data to make M models')

% amount of data reqired to create the model
n = params.n;
% the order fo the genereated system
m = params.m;

seglen = n + m + 1;
%%  Run algorithm 1 to gather relative data
ep_lo = 0.05;
ep_hi = 0.10;

Ppor = 0.80;
Ppst = 0.85;

x.dim = 2;
x.sample = zeros(2, M);
x.M = M;

[rstar, ~, Ntrial, q_min, q_max] = algo1(x, ep_lo, ep_hi, Ppor, Ppst);

%% Define the M models

Mmodels = cell(8, M);

entry.PhiP = 1;
entry.Bp = 2;
entry.Cp = 3;
entry.P = 4;
entry.H = 5;
entry.Ac = 6;
entry.Ax = 7;
entry.b1 = 8;
% [PhiP, Bp, Cp, P, H, Ac, Ax, b1] = makeModelandConstraints(y, u, p, params, bnds);
% use a parfor loop to create M models by calling the above command

%% flip the data, simplifies which part of the data to use
y_flip = fliplr(y);
u_flip = fliplr(u);

for i = 1: M
    
    % initialise the last n + m data, seglen = n+m+1
    
    y0 = (i-1)*seglen + 1;
    y1 = i*seglen;
    
    y_data = fliplr(y_flip(:, y0: y1));
    u_data = fliplr(u_flip(:, y0: y1));
    
    [PhiP, Bp, Cp, P, H, Ac, Ax, b1] = makeModelandConstraints(y_data, u_data, p, params, bnds);
    % Store the data in the cell array
    Mmodels{entry.PhiP, i} = PhiP;  Mmodels{entry.Bp, i} = Bp; Mmodels{entry.Cp, i} = Cp;
    Mmodels{entry.P, i} = P;  Mmodels{entry.H, i} = H; Mmodels{entry.Ac, i} = Ac;
    Mmodels{entry.Ax, i} = Ax;  Mmodels{entry.b1, i} = b1;
    
   
end

ck = NaN;
