function c = adaptiveControl(Q_bar, Y, Ck, p, params, bnds)
n = params.n;
m = params.m;
Ts = params.Ts;

assert(p/length(bnds) == round(p/length(bnds)), 'bounds need to repeat: p must be divisible by 3')
% takes data, generates model and finds optimal input c.

% perform least squares on cart velocity
[Yx, ~, Px, ~] = least_squares_params(Y(1, :), Ck, n, m);
% perform least squares on pendulum velocity
[Yphi, ~, Pphi, ~] = least_squares_params(Y(2, :), Ck, n, m);

[Apx, Bpx] = estSS(Px, m);
[Aphi, Bphi] = estSS(Pphi, m);
% fully augmented matrix
[Ap, Bp, Cp, Xp] = augment_estSS(Apx, Aphi, Bpx, Bphi, Yx, Yphi);
% Perform LQR on the generated matrix

eigen_values = eig(Ap);

[PhiP, Bp, K_opt] = optmal_gain(Ap, Bp, eye(length(Ap)), 0.1, 0);
% PhiP = Ap;
% K_opt = zeros(1, length(PhiP));
% find the constraints for optimisation
[Ac, Ax, b1] = make_constraints(PhiP, Bp, Cp, K_opt, p, bnds);

% Q is the Hessian for the otpimization
[L2, ~] = chol(Q_bar,'lower');
Linv = inv(L2);

b = b1 + Ax*Xp;

% options for mpcqpsolver:
options = mpcqpsolverOptions;

c = mpcqpsolver(Linv, zeros(p, 1), Ac, b, [], zeros(0,1), false(size(b)), options);


%% The functions inside adaptive control

function [Y, D, P, P_expanded] = least_squares_params(y, c, n, m)
% This function takes the input data and makes an input output data model
assert(n >= m, 'the number of rows should be greater than columns, necessary that n > m')
assert(length(y)>= n+m, 'y vector is too small: make n and m smaller'  )

% find y_0 index
y_0_index = length(y) - n - m;
% Y the column vector for data n long
Y = y(:, y_0_index + m + 1: end)';
Y = reshape(Y, [], 1);
% make sure y or c is repeated
if length(y(:,1)) > length(c(:,1))
   
    c2 = repmat(c, length(y(:,1))/length(c(:,1)), 1);
    y2 = y;
    
elseif length(y(:,1)) < length(c(:,1))
   
    y2 = repmat(y, length(c(:,1))/length(y(:,1)), 1);
    c2 = c;
else
    c2 = c;
    y2 = y;
    
end

D = NaN*zeros(n*length(y(:,1)), 2*m);
% build up the D matrix
nm_end = y_0_index + 1;

n_point = y_0_index + m;

no_outputs2 = length(y2(:,1));

for i = 1:n
    
   d_inter = fliplr(y2(:, nm_end + i: n_point + i));
   d_inter = [d_inter, fliplr(c2(:, nm_end + i: n_point + i))]; 
   
   D(i*no_outputs2 - no_outputs2 +1 : i*no_outputs2, :) = d_inter;
    
end

psudo = (D'*D)\D';

P = psudo*Y;

P_expanded = reshape(P, [], no_outputs2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ap, Bp] = estSS(P, m)

Ap = zeros(m);
Ap(2:m, 2:m) = eye(m-1);

last_row_A = flipud(P(1:m));
last_row_A = last_row_A';

Ap(end, :) = last_row_A;

Bp = zeros(m, 1);
Bp(end, 1) = 1;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ap, Bp, Cp, Xp] = augment_estSS(Ax, Aphi, Bx, Bphi, Yx, Yphi)

Ap = [Ax, zeros(length(Ax));
      zeros(length(Ax)), Aphi];
  
Bp = [Bx; Bphi];

Cp = eye(length(Ap));

Xp = [Yx(end - length(Ax) +1: end, :); Yphi(end - length(Ax) +1: end, :)];

end


function [PhiP, Bp, K_opt] = optmal_gain(Ap, Bp, Q, R, N)
% u = k_opt*x + c

K_opt = -1*dlqr(Ap,Bp,Q,R,N);
PhiP = Ap + Bp*K_opt;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ac, Ax, b] = make_constraints(PhiP, Bp, C, K_opt, p, bnds)
% Make psi
no_states = length(PhiP);

if p > 1
    
    eN = eye(p);
    eN = eN(1, :);
    
    M = [zeros(p-1, 1), eye(p-1);
         zeros(1, p)];    

    
    psi = [PhiP,                   Bp*eN;
          zeros(p, length(PhiP)),   M];
elseif p == 1
        
    eN = eye(p);
    eN = eN(1, :);
    
    psi = [PhiP,                   Bp*eN];
    psi = [psi; zeros(size(psi))];
         
end
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
    % inefficient but easy to make
    As = [As ; F*psi^i];
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

end