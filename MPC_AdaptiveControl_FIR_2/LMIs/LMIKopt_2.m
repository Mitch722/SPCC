
% Bias on the Model variance
bias = 0.53;
% noise width for the uniform dist
nWidth = 0.0007;

TsFast = 0.005;
TsObvs = 0.01;

M0 = 1.5;   M = M0;
m0 = 0.2;   m = m0;

[Phi, ~, L, K_opt] = inverted_pen_T(TsFast, M, m);

[Phis, ~, Ls, K_opts] = inverted_pen_T(TsObvs, M, m);

Af = Phi.A(1:4, 1:4) + Phi.B(1:4, :)*K_opt;
Bf = Phi.B(1:4, :);
Cf = Phi.C(:, 1:4);

As = Phis.A(1:4, 1:4) + Phis.B(1:4, :)*K_opts;
Bs = Phis.B(1:4, :);
Cs = Phis.C(:, 1:4);

% A_k = Af + nWidth*rand(4) - bias*nWidth*ones(4);
% A_k = As - [zeros(1, 4); zeros(3, 1), 0.0045*ones(3)];
A_k = As - 0.004*ones(4);
% A_k = Af - 0.0015*ones(4);
% robust set of Phi

Phi = [As-Bs*K_opts-Ls*Cs,   Ls*Cs;
       -Bs*K_opt,        A_k];

% define unknown matrix to be found using LMI
setlmis([])
P = lmivar(1, [size(Phi, 1) 1]);  % structure and size of LMI

% def LMI
lmiterm([1 1 1 P], Phi', Phi);         % Phi' P Phi
lmiterm([1 1 1 P], -1, 1);         % -P
lmiterm([1 1 2 0], 1);            % 0
lmiterm([1 2 2 P], -1, 1);         % -P

LMISYS = getlmis;

[tmin, Psol] = feasp(LMISYS);
P = dec2mat(LMISYS, Psol, P);

eigP = eig(P);
