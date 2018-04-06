% Work out all the functions for the parameter estimation
% do parameter estimation on the velocity components

%% 
try 
    try_statement = A;
    
catch
    MPC_inv_10
    
end


%% Try the least squares 
n = 50;
m = 10;
c = Ck;

y1 = y(1, 1: 1900);
c1 = c(1, 1: 1900);
[Y, D, P, P_expanded] = least_squares_params(y1, c1, n, m);

Y_hat = D*P;

e = Y - Y_hat;

%%
y1 = y(1, 1: 1900);
c1 = c(1, 1: 1900);

Ts = 0.01;

[Y, D, P, P_expanded] = velo_LS(y1, c1, n, m, Ts);

Y_hat = D*P;

e = Y - Y_hat;

%% Build state Space Matrix
[Apx, Bpx] = estSS(P, m);
[Aphi, Bphi] = estSS(P, m);

[Ap, Bp] = augment_estSS(Apx, Aphi, Bpx, Bphi);

% [PhiP, Bp, K_opt] = optmal_gain(Ap, Bp, Q, R, N);

[H, f, Ac, Ax, b, lb, ub, options] = MPC_vars(A, B, C, K_opt, R, p, bnds, maxF);

%% Functions

function [Y, D, P, P_expanded] = least_squares_params(y, c, n, m)
% This function takes the input data and makes an input output data model
assert(n > m, 'the number of rows should be greater than columns, necessary that n > m')
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
%%
function [Y, D, P, P_expanded] = velo_LS(y, c, n, m, Ts)
% This function takes the input data and makes an input output data model
% input args:
%   1. y: Sensor outputs
%   2. c: input sequence
%   3. n: length of filter
%   4. m: outputs to perform least squares
%   5. Ts: sampling time
% returns:
%   1. 
assert(n > m, 'the number of rows should be greater than columns, necessary that n > m')
assert(length(y)>= n+m, 'y vector is too small: make n and m smaller'  )

%% Euler to make output derivatives

H1 = (1/Ts)*eye(length(y));
H1(end, end) = 0;

H2 = (-1/Ts)*eye(length(y) - 1);

H2 = [zeros(length(H2), 1), H2;
      zeros(1, length(y))];
  
DH = H1 + H2;

% change y into derivatives and make it a long vector
y = (DH*y')';
% trim off the end zero
y(:, end) = [];
%%


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

%%
function [Ap, Bp] = estSS(P, m)

Ap = zeros(m);
Ap(2:m, 2:m) = eye(m-1);

last_row_A = flipud(P(1:m));
last_row_A = last_row_A';

Ap(end, :) = last_row_A;

Bp = zeros(m, 1);
Bp(end, 1) = 1;


end
%%
function [Ap, Bp] = augment_estSS(Ax, Aphi, Bx, Bphi)

Ap = [Ax, zeros(length(Ax));
      zeros(length(Ax)), Aphi];
  
Bp = [Bx; Bphi];

end


function [PhiP, Bp, K_opt] = optmal_gain(Ap, Bp, Q, R, N)
% u = k_opt*x + c

K_opt = -1*dlqr(A,B,Q,R,N);
PhiP = Ap + Bp*K_opt;

end

function [Q_bar, psi] = solveLyapunov(PhiP, Bp, Q, K_opt, R, E, p)

if p > 1
    
    eN = eye(p);
    eN = eN(1, :);
    
    M = [zeros(p-1, 1), eye(p-1);
         zeros(1, p)];    

    
    psi = [PhiP,                   Bp*eN;
          zeros(p, length(phi)),   M];
elseif p == 1
        
    eN = eye(p);
    eN = eN(1, :);
    
    psi = [PhiP,                   Bp*eN];
    psi = [psi; zeros(size(psi))];
         
end

Q2 = [Q,                zeros(size(Q));
      zeros(size(Q)),   zeros(size(Q))];
  
Q2 = Q2 + [K_opt'*R*K_opt,  K_opt'*R*eN;
           eN'*R*K_opt,     eN'*R*eN];
% Q_bar is the solution to the discrete Lyapunov equation
Q_bar = dlyap(psi, Q2);

end

function [Ac, Ax, b] = make_constraints(PhiP, Bp, C, K_opt, p, bnds)
% Make psi
if p > 1
    
    eN = eye(p);
    eN = eN(1, :);
    
    M = [zeros(p-1, 1), eye(p-1);
         zeros(1, p)];    

    
    psi = [PhiP,                   Bp*eN;
          zeros(p, length(phi)),   M];
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
for i = 1: length(As)/length(bnds) -1
    bounds2 = [bounds2; bnds];
end
    
b = -1*bounds2;

end
