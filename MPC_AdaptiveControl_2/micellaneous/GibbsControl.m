function [c, Ap, Bp, Cp] = GibbsControl(Q_bar, Y, Ck, p, params, bnds)
n = params.n;
m = params.m;
Ts = params.Ts;

assert(p/length(bnds) == round(p/length(bnds)), 'bounds need to repeat: p must be divisible by 3')
% takes data, generates model and finds optimal input c.

[sys_obv, ~, K_opt] = inverted_pen_T(Ts, M, m);

Ap = sys_obv.Aobv;
Bp = sys_obv.Bov;
Cp = sys_obv.Cobv;

% find the constraints for optimisation
[H, f, Ac, Ax, b, lb, ub, options] = MPC_vars(A, B, C, K_opt, R, p, bnds, maxF)

% Q is the Hessian for the otpimization
[L2, ~] = chol(H,'lower');
Linv = inv(L2);

% generate states from the output for a fourth order system
Xp = getState_n_4(Y, Ck, Ap, Bp, Cp);

% The acutal correct input
b = b1 + Ax*Xp;

% options for mpcqpsolver:
options = mpcqpsolverOptions;

c = mpcqpsolver(Linv, zeros(p, 1), Ac, b, [], zeros(0,1), false(size(b)), options);
% c = quadprog(Q_bar, f, -Ac, -b, [], [], lb, ub, [], options);

%% The functions inside adaptive control

function [Y, D, P, P_expanded] = least_squares_params(y, c, n, m)
% This function takes the input data and makes an input output data model
assert(n >= m, 'the number of rows should be greater than columns, necessary that n > m')
assert(length(y)>= n+m, 'y vector is too small: make n and m smaller'  )

% find y_0 index
y_0_index = length(y) - n - m;
% Y the column vector for data n long
Y = y(:, y_0_index + m + 1: end);
Y = reshape(Y, [], 1);
% make sure y or c is repeated

y2 = y;
c2 = c;

D = NaN*zeros(n*length(y(:,1)), 3*m + 2);
% build up the D matrix
nm_end = y_0_index;

n_point = y_0_index + m - 1;

no_outputs2 = length(y2(:,1));

for i = 1:n
    
   d_inter = fliplr(y2(:, nm_end + i: n_point + i));
   
   c_input = fliplr(c2(:, nm_end + i: n_point + 1 + i)); 
   c_inter = [c_input, zeros(size(c_input))];
   
   c_inter = [c_inter; zeros(size(c_input)), c_input];
   
   d_inter = [d_inter, c_inter]; 
   
   D(i*no_outputs2 - no_outputs2 +1 : i*no_outputs2, :) = d_inter;
    
end

psudo = (D'*D);

P = psudo\D'*Y;

P_expanded = reshape(P, no_outputs2, []);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ap, Bp, Cp] = estSS(P, m)

% take only the parameters for the states:
m = m;
    
Ap = zeros(m);
Ap(2:m, 2:m) = eye(m-1);

last_row_A = fliplr(P(1:m)');

Ap(end, :) = (last_row_A);
% Ap(end-1, :) = last_row_A;

Bp = zeros(m, 1);
% Bp(end-1 : end, 1) = [1; 1];

Bp(end, 1) = 1;

Cp_inter = P(m+2: end, :)';
Cp_inter(m+1) = [];
% flip to put correct coefficients in the correct place
Cp_inter = fliplr(Cp_inter);
% Cp = zeros(2, m);

Cp = [Cp_inter(:, 1:m); Cp_inter(:, m+1:end)];

Cp2 = P(m+1, 1).*[last_row_A(1, 1:m)];
Cp2 = [Cp2; P(2*m+2, 1).*last_row_A(1, 1:m)];

Cp = Cp + Cp2;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ap, Bp, Cp, Xp] = augment_estSS(Ax, Aphi, Bx, Bphi, Cx, Cphi, Yx, Yphi)

Ap = [Ax, zeros(length(Ax));
      zeros(length(Ax)), Aphi];
  
Bp = [Bx; Bphi];

Cp = [Cx,            zeros(1, length(Cx));
      zeros(1, length(Cx)),     Cphi];

Xp = [Yx(end - length(Ax) +1: end, :); Yphi(end - length(Ax) +1: end, :)];

end


function [PhiP, Bp, K_opt] = optmal_gain(Ap, Bp, Q, R, N)
% u = k_opt*x + c

K_opt = -1*dlqr(Ap,Bp,Q,R,N);
PhiP = Ap + Bp*K_opt;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PhiP, Bp, K] = place_poles(Ap, Bp, P)
    K = place(Ap,Bp,P);
    
    PhiP = Ap - Bp*K;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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