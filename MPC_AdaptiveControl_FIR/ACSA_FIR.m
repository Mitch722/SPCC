function [ck, Mmodels, RstarModel, entry, no_coefs, H] = ACSA_FIR(M, y, u, p, params, bnds, Phi, B, C, K_opt)


% check that there is enough data
assert(length(y) >= M*(params.n+params.m), 'There is not enough data to make M models')

% amount of data reqired to create the model
n = params.n;
% the order fo the genereated system
m = params.m;

seglen = n + m + 1;
%%  Run algorithm 1 to gather relative data
% Scenario approach 
ep_lo = 0.10;
ep_hi = 0.20;

Ppor = 0.80;
Ppst = 0.85;

x.dim = 1;
x.sample = zeros(2, M);
x.M = M;

[rstar, ~, Ntrial, q_min, q_max] = algo1(x, ep_lo, ep_hi, Ppor, Ppst);

rstar = round(rstar);
Ntrial = round(Ntrial);

assert(M > rstar*Ntrial, 'The number of samples is too few to sample from')

%% flip the data, simplifies which part of the data to use
y_flip = fliplr(y);
u_flip = fliplr(u);

Mmodels = cell(4, M);

% entry for calling data in a structure
entry.Ffir = 1;
entry.Dc = 2;
entry.Dx = 3;
entry.b0 = 4;
entry.bigFIR = 5;

model.A = Phi;
model.B = B;
model.C = C;
model.K_opt = K_opt;


for i = 1: M
    
    % initialise the last n + m data, seglen = n+m+1
    
    y0 = (i-1)*seglen + 1;
    y1 = i*seglen;
    
    y_data = fliplr(y_flip(:, y0: y1));
    u_data = fliplr(u_flip(:, y0: y1));
    
    % build FIR coefficients
    Fy = buildFIR(y_data(1, :), u_data, n, m);
    Fp = buildFIR(y_data(2, :), u_data, n, m);
    % put the two coefs together
    Ffir = [Fy'; Fp'];
    Ffir = fliplr(Ffir);
    % Store the FIR coefficients cell array
    Mmodels{entry.Ffir, i} = Ffir;
    Mmodels{entry.bigFIR, i} = ConvFIRmat(Ffir, p);
    
    % FIR constraints
    [b0, Dx, Dc, no_coefs, H] = MPC_FIR_ck(model, Ffir, p, bnds);
   
    Mmodels{entry.Dc, i} = Dc;
    Mmodels{entry.Dx, i} = Dx;
    Mmodels{entry.b0, i} = b0;
    
    
end

%% Randomly Sample Rstar * Ntrial Models for performing optimization over

[~, rowMmod] = size(Mmodels);

random_entries = randperm(rowMmod, rstar*Ntrial);
% shuffle models to create a sample 
sampleModels = Mmodels(:, random_entries);

%% Calculate the Optimal value of Cstar_i for all the sample Models
% Make constraints Ac*ck <= b1 + Ax*Xp 

% cell array to store optimal references
cStar = cell(2, Ntrial);
% MPC references
% make 
[L2, ~] = chol(H, 'lower');
Linv = inv(L2);
% options for mpcsolver
options = mpcqpsolverOptions;
% last few steps from reference-input Ck
chat = u(:, end-no_coefs+1:end)';
assert(length(chat(:,1))==no_coefs, 'The above line is incorrect and probs has one extra coeficient')

for i = 1 : Ntrial
    Dc0 = (i-1)*rstar + 1;
    Dc1 = i*rstar;
    
    Dc = sampleModels(entry.Dc, Dc0: Dc1)';
    Dc = cell2mat(Dc);
    
    Dx = sampleModels(entry.Dx, Dc0: Dc1)';
    Dx = cell2mat(Dx);
        
    b1 = sampleModels(entry.b0, Dc0: Dc1)';
    b1 = cell2mat(b1);
    
%     Xp = sampleModels(entry.Xp, Ac0: Ac1)';
%     Xp = cell2mat(Xp);
    
    b = b1 + Dx*chat;
    % solve the quadprog for c
    c = mpcqpsolver(Linv, zeros(p, 1), Dc, b, [], zeros(0,1), false(size(b)), options);
    cStar{1, i} = c;
    
end
%% Calcualte the number of violations of all the other inputs
% Simulation stage
rStarSubSamps = reshape(random_entries, [], rstar);

rStarModels = Mmodels;
% physical bounds on outputs
% bndlims = repmat(bnds(1:2, :), (M-rstar), 1);
totalNoViol = 2*p*(M-rstar);
% a stacked vector of Uk inputs
Ds = vecUk(Phi, B, K_opt, p);
% Uk = Ds * zk;
%% Find the number of violations
for i = 1: Ntrial 
   
   % px1 vector if input references 
   ck = cStar{1, i};
   % initial state plus 
   zk = ck;
         
   rModels = rStarModels;
   rModels(:, rStarSubSamps(1, :)) = []; 
   
   bigF = rModels(entry.bigFIR, :)';
   bigFmat = cell2mat(bigF);
   
   yFir = bigFmat*zk;
   
   bndlims = repmat(bnds, length(yFir)/length(bnds), 1);
   
   violations = sum((abs(yFir) > bndlims), 2);
   violations = sum(violations, 1);
   
   cStar{2, i} = 1 - violations/totalNoViol;
   
end
%% choose the otpimal value for ck 
compare = 0.5*(q_min + q_max)/M;

diff = compare - cell2mat(cStar(2, :));
[~, opt_index] = min(abs(diff));

% choose the optimal input which gives the optimal number of violations
ck = cStar{1, opt_index};

RstarModel = sampleModels(:, (opt_index-1)*rstar + 1 : (opt_index)*rstar - 1);


