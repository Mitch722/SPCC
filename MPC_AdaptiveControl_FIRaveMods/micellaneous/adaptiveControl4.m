function [c, Ap, Bp, Cp, P] = adaptiveControl4(H, Y, Ck, p, params, bnds, Xp)
n = params.n;
m = params.m;
Ts = params.Ts;

assert(p/length(bnds) == round(p/length(bnds)), 'bounds need to repeat: p must be divisible by 3')
% takes data, generates model and finds optimal input c.

[~, ~, P, ~] = least_squares_params(Y, Ck, n, m);

[Ap, Bp, Cp] = estSS(P, m);

Qp = Cp'*Cp;

K_opt = dlqr(Ap, Bp, Qp, 1, 0);
% Q_bar = dlyap(Ap2, Cp2'*Cp2);
PhiP = Ap - Bp*K_opt;

% find the constraints for optimisation
[Ac, Ax, b1] = make_constraints(PhiP, Bp, Cp, K_opt, p, bnds);

% Q is the Hessian for the otpimization
[L2, ~] = chol(H,'lower');
Linv = inv(L2);

% The acutal correct input
b = b1 + Ax*Xp;

% options for mpcqpsolver:
options = mpcqpsolverOptions;

c = mpcqpsolver(Linv, zeros(p, 1), Ac, b, [], zeros(0,1), false(size(b)), options);
% c = quadprog(Q_bar, f, -Ac, -b, [], [], lb, ub, [], options);

%% The functions inside adaptive control

function [Ac, Ax, b] = make_constraints(PhiP, Bp, C, K_opt, p, bnds)
% Make psi
no_states = length(PhiP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
E = [1, zeros(1, p-1)];

[len_output, ~] = size(C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
F = [C, zeros(len_output, p);
     K_opt, E];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% initialise As
As = F;
for i = 1 : p-1
    
    psi = makePsi(PhiP, Bp, p, i);
    % inefficient but easy to make
    As = [As ; F*psi];
    % A( (i)*ax + 1 : (i + 1)*ax, :) = F*psi^i; 
end
% split As into Ac and Ax
Ac = As(:, (no_states+1):end);
Ax = As(:, 1: no_states);

% repeat A
Ac = [Ac; -Ac];
Ax = [Ax; -Ax];
% Make A negative
Ac = -Ac;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bounds2 = bnds;
for i = 1: length(Ac)/length(bnds) -1
    bounds2 = [bounds2; bnds];
end
    
b = -1*bounds2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end