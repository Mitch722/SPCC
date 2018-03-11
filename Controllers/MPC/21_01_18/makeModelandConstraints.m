function [PhiP, Bp, Cp, P, H, Ac, Ax, b1] = makeModelandConstraints(Y, Ck, p, params, bnds)

% This function is based upon adaptiveControl5.m 
% Functionality:
% 1. Takes input/output data and makes a controllable model
% 2. Stablises Model using LQR 
% 3. Generaters a Hessian, Constraints for the input and for propagating
% States through model
% 

% Define the model parameters 
n = params.n;
m = params.m;

assert(p/length(bnds) == round(p/length(bnds)), 'bounds need to repeat: p must be divisible by 3')
% takes data, generates model and finds optimal input c.

[~, ~, P, ~] = least_squares_params(Y, Ck, n, m);

[Ap, Bp, Cp] = estSS(P, m);

Qp = Cp'*Cp;

K_opt = dlqr(Ap, Bp, Qp, 1, 0);
% Q_bar = dlyap(Ap2, Cp2'*Cp2);
PhiP = Ap - Bp*K_opt;

% find the constraints for optimisation
% [Ac, Ax, b1] = make_constraints(PhiP, Bp, Cp, K_opt, p, bnds);
[H, ~, Ac, Ax, b1, ~, ~, ~] = MPC_vars(PhiP, Bp, Cp, K_opt, 1, p, bnds, 100);


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