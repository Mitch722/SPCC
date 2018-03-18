function [PhiP, Bp, Cp, P, H, Ac, Ax, b1, K_opt] = makeModelandConstraints(Y, Ck, p, params, bnds)

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
%[~, ~, P, ~] = least_squares_params_new(Y, Ck, n, m);

[Ap, Bp, Cp] = estSS(P, m);

Qp = Cp'*Cp;
try
    K_opt = dlqr(Ap, Bp, Qp, 1, 0);
catch
    K_opt = zeros(1, length(Ap));
end
% Q_bar = dlyap(Ap2, Cp2'*Cp2);
PhiP = Ap - Bp*K_opt;

% find the constraints for optimisation
[H, ~, Ac, Ax, b1, ~, ~, ~] = MPC_vars(PhiP, Bp, Cp, K_opt, 1, p, bnds, 100);


end