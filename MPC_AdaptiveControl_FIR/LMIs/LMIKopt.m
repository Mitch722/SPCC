
% Bias on the Model variance
bias = 0.53;
% noise width for the uniform dist
nWidth = 0.05;

TsFast = 0.005;

M0 = 1.5;   M = M0;
m0 = 0.2;   m = m0;

[Phi, ~, L, K_opt] = inverted_pen_T(TsFast, M, m);
Phi = Phi.A;
Phi = Phi(1:4, 1:4) + nWidth*rand(4) - bias*nWidth*ones(4);
% robust set of Phi
Phi_m = Phi + 0.1*ones(4); % + Phi*0.995;

% define unknown matrix to be found using LMI
setlmis([])
P = lmivar(1, [size(Phi, 1) 1]);  % structure and size of LMI

% def LMI
lmiterm([1 1 1 P], -1, 1, 's');         % -P
lmiterm([1 1 2 1], 1, Phi_m');            % Phi'
lmiterm([1 2 2 inv(P)], -1, 1);         % -inv(P)

LMISYS = getlmis;

[tmin, Psol] = feasp(LMISYS);
P = dec2mat(LMISYS, Psol, P);

eigP = eig(P);
