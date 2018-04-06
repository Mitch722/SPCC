function [uk, c] = RmodelOutput(Q_bar, RstarModel, entry, y, u, p)

Ac = RstarModel(entry.Ac, :)';
Ac = cell2mat(Ac);

Ax = blkdiag(RstarModel{entry.Ax, :});

b1 = RstarModel(entry.b1, :)';
b1 = cell2mat(b1);


[~, rowsAx] = size(Ax);
[~, Rstar] = size(RstarModel);

noStates = rowsAx/Rstar; 

Xp = zeros(noStates*Rstar, 1);

% entry.PhiP = 1;
% entry.Bp = 2;
% entry.Cp = 3;
% entry.P = 4;
% entry.H = 5;
% entry.Ac = 6;
% entry.Ax = 7;
% entry.b1 = 8;
% entry.K_opt = 9;

for i = 1: Rstar
    
    Phi = RstarModel{entry.PhiP, i};
    B = RstarModel{entry.Bp, i};
    C = RstarModel{entry.Cp, i};
    
    Xp0 = (i-1)*noStates +1;
    Xp1 = i*noStates;
    
    Xp( Xp0: Xp1, 1) = getState_n_4(y, u, Phi, B, C);
    
end

b = b1 + Ax*Xp;


[L2, ~] = chol(Q_bar, 'lower');
Linv = inv(L2);
options = mpcqpsolverOptions;

c = mpcqpsolver(Linv, zeros(p, 1), Ac, b, [], zeros(0,1), false(size(b)), options);

K_opt = cell2mat(RstarModel(entry.K_opt, :));

uk = K_opt*Xp + c(1);

