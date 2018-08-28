 
load('output_input.mat')
    


y(:, end) = [];
Ck(:, end) = [];

M = 80;
p = 12;

params.m = 4;
params.n = 10;

Q_bar = eye(p);
bnds = [0.8, 1.5, 0.6]';

Ts = 0.01;
%Xp = getState_n_4(Y, Ck, PhiP, Bp, Cp);
Xp = randn(4, 1);

tic
[ck, Mmodels, RstarModel, entry] = ACSA_1(M, y, Ck, p, params, Q_bar, bnds, Xp);
toc

%% test RmodelOutput 
tic

[uk, c] = RmodelOutput(Q_bar, RstarModel, entry, y, Ck, p);


toc
%%
% mod2 = Mmodels(:, 1:3);

[colMmod, rowMmod] = size(Mmodels);

random_entries = randperm(rowMmod, 3);
sampleModels = Mmodels(:, random_entries)';

Ax = sampleModels(:, 6);
Ax = cell2mat(Ax);

sampleModels = sampleModels';

bigPhi = sparse(blkdiag(sampleModels{1, :}));
sparseBigPhi = sparse(bigPhi);

bigXp = repmat(Xp, 3, 1);

dp = sparseBigPhi* bigXp;
dp1 = bigPhi*bigXp;

%% choose optimal Phip
% find the average of all the models

[noItems, noRstar] = size(RstarModel);
aveModel = (1/noRstar)*sum(cat(3, RstarModel{1, :}), 3 );
