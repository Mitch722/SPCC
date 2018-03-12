function [H, f, Ac, Ax, b, lb, ub, options] = MPC_vars(A, B, C, K_opt, R, p, bnds, maxF)
% This function makes variables for quadprog function for MPC 


%% begin solution
% A is A + B*K_opt
phi = A;
[~, no_states] = size(A);

K_opt = [K_opt, zeros(1, no_states - length(K_opt))];

M = [zeros(p-1, 1), eye(p-1);
     zeros(1, p)];    

eN = eye(p);
eN = eN(1, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make psi
psi = [phi,                     B*eN;
       zeros(p, length(phi)),   M];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
E = [1, zeros(1, p-1)];

[len_output, ~] = size(C);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
F = [C, zeros(len_output, p);
     K_opt, E];
 
% F = [C, zeros(len_output, p)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
Q = C'*C;

Q_bar2 = [K_opt'*R*K_opt,   K_opt'*R*E;
          E'*R*K_opt,       E'*R*E];

[size_Qx, ~] = size(Q);
Q_bar = Q_bar2;
Q_bar(1: size_Qx, 1: size_Qx) = Q(1:size_Qx, 1: size_Qx) + Q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solve discrete Lyapunov equation:
P_bar = dlyap(psi', Q_bar);
% Construct -As x <= b, the -As part:

% [ax, ay] = size(F);

% initialise As
As = F;
for i = 1 : p-1
    % inefficient but easy to make
    As = [As ; F*psi^i];
    % A( (i)*ax + 1 : (i + 1)*ax, :) = F*psi^i; 
end
% split As into Ac and Ax
Ac = As(:, (no_states+1):end);
Ax = As(:, 1: no_states);

% repeat A
As = [As; -As];
Ac = [Ac; -Ac];
Ax = [Ax; -Ax];
% Make A negative
As = -As;
Ac = -Ac;

%% Construct -Ax <= -b, the -b part:

bounds2 = bnds;
for i = 1: length(As)/length(bnds) -1
    bounds2 = [bounds2; bnds];
end
    
b = -1*bounds2;
%% Find H = Pcc from P_bar 
% Takes bottom right hand corner of P_bar this is the hessian
H = P_bar(no_states+1 : end, no_states+1 : end);

options = optimoptions('quadprog');
options.ConstraintTolerance = 1e-8;

%% Use equality constraints to feed in current position
% f is a vector of zeros
f = zeros(1, p)';

lb = -maxF*ones(p, 1);
ub = maxF*ones(p, 1);
% x = quadprog(H,f,-A,-b,Aeq,beq,lb,ub,x0,options)

end