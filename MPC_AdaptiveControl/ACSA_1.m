function [ck, Mmodels, RstarModel, entry] = ACSA_1(M, y, u, p, params, Q_bar, bnds, Xp)


% check that there is enough data
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

x.dim = 1;
x.sample = zeros(2, M);
x.M = M;

[rstar, ~, Ntrial, q_min, q_max] = algo1(x, ep_lo, ep_hi, Ppor, Ppst);

rstar = round(rstar);
Ntrial = round(Ntrial);

assert(M > rstar*Ntrial, 'The number of samples is too few to sample from')
%% Define the M models

Mmodels = cell(9, M);

entry.PhiP = 1;
entry.Bp = 2;
entry.Cp = 3;
entry.P = 4;
entry.H = 5;
entry.Ac = 6;
entry.Ax = 7;
entry.b1 = 8;
entry.K_opt = 9;

% [PhiP, Bp, Cp, P, H, Ac, Ax, b1, K_opt] = makeModelandConstraints(y, u, p, params, bnds);
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
    
    [PhiP, Bp, Cp, P, H, Ac, Ax, b1, K_opt] = makeModelandConstraints(y_data, u_data, p, params, bnds);
    % Store the data in the cell array
    Mmodels{entry.PhiP, i} = PhiP;  Mmodels{entry.Bp, i} = Bp; Mmodels{entry.Cp, i} = Cp;
    Mmodels{entry.P, i} = P;  Mmodels{entry.H, i} = H; Mmodels{entry.Ac, i} = Ac;
    Mmodels{entry.Ax, i} = Ax;  Mmodels{entry.b1, i} = b1; Mmodels{entry.K_opt, i} = K_opt;
    
   
end

%% Randomly Sample Rstar * Ntrial Models for performing optimization over

[~, rowMmod] = size(Mmodels);

random_entries = randperm(rowMmod, rstar*Ntrial);
sampleModels = Mmodels(:, random_entries);

%% Calculate the Optimal value of Cstar_i for all the sample Models
% Make constraints Ac*ck <= b1 + Ax*Xp 

% cell array to store optimal references
cStar = cell(2, Ntrial);
% MPC references
% make 
[L2, ~] = chol(Q_bar, 'lower');
Linv = inv(L2);
% options for mpcsolver
options = mpcqpsolverOptions;

for i = 1 : Ntrial
    Ac0 = (i-1)*rstar + 1;
    Ac1 = i*rstar;
    
    Ac = sampleModels(entry.Ac, Ac0: Ac1)';
    Ac = cell2mat(Ac);
    
    Ax = sampleModels(entry.Ax, Ac0: Ac1)';
    Ax = cell2mat(Ax);
    
    b1 = sampleModels(entry.b1, Ac0: Ac1)';
    b1 = cell2mat(b1);
    
    b = b1 + Ax*Xp;
    
    c = mpcqpsolver(Linv, zeros(p, 1), Ac, b, [], zeros(0,1), false(size(b)), options);
    cStar{1, i} = c;
    
end
%% Calcualte the number of violations of all the other inputs
rStarSubSamps = reshape(random_entries, [], rstar);

rStarModels = Mmodels;
% physical bounds on outputs
bndlims = repmat(bnds(1:2, :), (M-rstar), 1);

for i = 1: Ntrial
    
   Xv = zeros(4*(M-rstar), p+1);
   Xv(:, 1) = repmat(Xp, (M-rstar), 1);
   Yv = zeros(2*(M-rstar), p+1);
    
   % px1 vector if input references 
   ck = cStar{1, i}';
   rModels = rStarModels;
   rModels(:, rStarSubSamps(1, :)) = []; 
   % Perform all the models in one go
   bigPhi = sparse(blkdiag(rModels{entry.PhiP, :}));
   bigB = sparse(cell2mat( rModels(entry.Bp, :)' ));
   bigC = blkdiag(rModels{entry.Cp, :});
   
   
   for k = 1: p
       
      Xv(:, k+1) = bigPhi * Xv(:, k) + bigB*ck(1, k);
      Yv(:, k)  = bigC * Xv(:, k);
      
   end
   
   violations = sum((Yv > bndlims), 2);
   cStar{2, i} = sum(violations, 1);
   
end
%% choose the otpimal value for ck 
compare = 0.5*(q_min + q_max);

diff = compare - cell2mat(cStar(2, :));
[~, opt_index] = min(abs(diff));

% choose the optimal input which gives the optimal number of violations
ck = cStar{1, opt_index};

RstarModel = sampleModels(:, (opt_index-1)*rstar + 1 : (opt_index)*rstar - 1);


