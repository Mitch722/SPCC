try 
    
    a  = y;
    
catch
    
    load('output_input.mat')
    
end

M = 80;
p = 12;

params.m = 4;
params.n = 10;

Q_bar = eye(p);
bnds = [0.8, 1.5, 0.6]';

Ts = 0.01;

[~, Mmodels] = ACSA_1(M, y, Ck, p, params, Q_bar, bnds);
%%
% mod2 = Mmodels(:, 1:3);

[colMmod, rowMmod] = size(Mmodels);

random_entries = randperm(rowMmod, 3);
sampleModels = Mmodels(:, random_entries)';

