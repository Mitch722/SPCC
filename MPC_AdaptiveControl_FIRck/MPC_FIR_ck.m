% Aim: [Y; -Y] =< [H; -H] c_k:
% build 
%
%
function [bx, Acp, Acc, no_coefs, H] = MPC_FIR_ck(model, F, p, b)
% 
% find Hessian
A = model.A;
B = model.B;
C = model.C;
K_opt = model.K_opt;

phi = A;
R = 1;

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

Q = C'*C;

Q_bar2 = [K_opt'*R*K_opt,   K_opt'*R*E;
          E'*R*K_opt,       E'*R*E];

[size_Qx, ~] = size(Q);
Q_bar = Q_bar2;
Q_bar(1: size_Qx, 1: size_Qx) = Q(1:size_Qx, 1: size_Qx) + Q;
% solve discrete Lyapunov equation:
P_bar = dlyap(psi', Q_bar);

H = P_bar(no_states+1 : end, no_states+1 : end);

no_outs = length(F(:, 1));
%% bounds
bnds = repmat(b, p, 1);

bx = [bnds; bnds];

%% Ax
% no. of coefficients
no_coefs = length(F(1, :));

Hc = ConvFIRmat(F, p);

if length(Hc(:,1)) >= p*no_outs
    Hc = Hc(1:p*no_outs, :);
else
   
    diff = p*no_outs - length(Hc(:,1));
    Hc = [Hc; zeros(diff, length(Hc(1,:)) )];
end

Hx = ConvFIRmat(F, no_coefs);

if length(Hx(:, 1)) < length(Hc(:, 1))
    
    diff = length(Hc(:, 1)) - length(Hx(:, 1));
    Hx = [Hx; zeros(diff, length(Hx(1, :) ))];
end


%% Split-up the Ax matrix

Acc = [Hc; -Hc];
Acp = [Hx; -Hx];
end
