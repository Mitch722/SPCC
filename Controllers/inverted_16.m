%% Initialise physical model and parameters for controller
% Using QuadProg instead of mpcqpsolver


[sys_obv, L, K_opt] = inverted_pen;

A = sys_obv.A;
B = sys_obv.B;
C = sys_obv.C;
D = sys_obv.D;
Ts = sys_obv.Ts;

[~, no_states] = size(A);
[no_outputs, ~] = size(D);

t = 0: Ts: 20;

% Preallocate data
x = zeros( no_states, length(t)+1 );
y = zeros( no_outputs, length(t)+1 );

x_max = [2; 0.3166; 5; 0.644];
x_min = -x_max;

u_over_bar = 10;
u_under_bar = -10;

% should put some weighting in here
Q = C'*C;
R = 1;

%% Find MPC optimal

% z(k|k) = [x(k|k) ck]'
% z(k+i|k) = psi^i * z(k|k)

% horizon window length
p = 15;
R = 1;
% bounds on 
bounds = [1, 0.1]';
% bounds = [bounds; bounds];

Y = [1; 1];
X = randn(8, 1) + 1;
X(2) = 0.1*X(2);

X = C*X;
% X = 0;
%%
Time_out = 10;
x = zeros(no_states, Time_out/Ts);
y = zeros(no_outputs, Time_out/Ts);

x2 = x;
y2 = y;

Ck = zeros(1, Time_out/Ts);

for k = 1: (Time_out/Ts)-1 
    
    Y = y(:, k);
    [ck, status] = optimal_input(A, B, C, Y, K_opt, R, p, bounds);
    
    c = ck(1);
    
%     if c^2 >= 100^2
%         c = 0;
%     end
    
    w =  0.1*randn(no_states, 1);
    v =  0.01*randn(no_outputs, 1);
    
    x(:, k+1) = A*x(:, k) + B*c + w;
    y(:, k) = C*x(:, k) + v;
    
    x2(:, k+1) = A*x2(:, k) + w;
    y2(:, k) = C*x2(:, k) + v;
    
    
    Ck(k) = c;
end

fprintf('Status = %d \n', status)
%%
figure
plot(y(1, :))
hold on
plot(y(2, :));
hold on

plot(y2(1, :))
hold on
plot(y2(2, :));
figure
plot(Ck)

%% Function for working out constrained MPC optimal input
function [ck, status] = optimal_input(A, B, C, Y, K_opt, R, p, bounds)

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
     K_opt*eye(no_states), E];
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
%% As = zeros(ax*p, ay):

% initialise As
As = F;
for i = 1 : p-1
    % inefficient but easy to make
    As = [As ; F*psi^i];
    % A( (i)*ax + 1 : (i + 1)*ax, :) = F*psi^i; 
end
% repeat A
Ac = As(:, (no_states+1):end);

As = [As; -As];
Ac = [Ac; Ac];
% Make A negative
As = -As;
Ac = -Ac;

% Construct -Ax <= -b, the -b part:
bounds2 = bounds;
for i = 1: length(As)/length(bounds) -1
    bounds2 = [bounds2; bounds];
end
    
b = -1*bounds2;

%% find Linv, the cholsky of P_bar:
N = p+1;
H = zeros((N+1)*length(Q_bar));

[qx, qy] = size(Q_bar);


%% initialise H
H(1:qy, 1:qx) = Q_bar;
for i = 1: N-1
    psi_trans = psi';
    H(i*qy + 1 : (i+1)*qy, i*qx + 1:(i+1)*qx) = (psi_trans^i)*Q_bar*psi^i;
end

H(N*qy +1: end, N*qx +1 : end) = P_bar;

%% Destroy H:

% H = R*eye(p);

[L, ~] = lu(H);
Linv = inv(L);

iA0 = false(size(b));
% iA0(1:2, :) = [true; true];

opt = mpcqpsolverOptions;
opt.IntegrityChecks = false;

%% Use equality constraints to feed in current position

Aeq = [C, zeros(len_output, p)];
%Aeq = [Aeq; zeros(1, (no_states+p))];
% Pad y
% Y = [Y; zeros(length(Aeq) - length(Y),1)];

% Find ck

% quadprog((H,f,A,b,Aeq,beq,lb,ub,x0,options)
ck = quadprog(H,zeros(length(As), 1),As,b,[],[],[],[]);

end
