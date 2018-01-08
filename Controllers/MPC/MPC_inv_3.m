%% Initialise physical model and parameters for controller
% This script is implements the correct Lyapunov hessian weight P_bar
% Uses quadprog instead of mpcqpsolver
% separates MPC hessians from quadprog

[sys_obv, L, K_opt] = inverted_pen;

A = sys_obv.A;
B = sys_obv.B;
C = sys_obv.C;
D = sys_obv.D;
Ts = sys_obv.Ts;

[~, no_states] = size(A);
[no_outputs, ~] = size(D);


%% Find MPC optimal

% z(k|k) = [x(k|k) ck]'
% z(k+i|k) = psi^i * z(k|k)

% horizon window length
p = 2;

Q = C'*C;
R = 1;
% bounds on 
main_bounds = [1, 1]';
% bounds = [bounds; bounds];

%%
Time_out = 20;
x = zeros(no_states, Time_out/Ts);
y = zeros(no_outputs, Time_out/Ts);

x2 = x;
y2 = y;

Ck = zeros(1, Time_out/Ts);

X = x(:, 1);
maxF = 100;

[H, f, Ac, Ax, b1, lb, ub, options] = MPC_vars(A, B, C, K_opt, R, p, main_bounds, maxF);

for k = 1: (Time_out/Ts)-1 
    
    b = b1 + Ax*X;
    ck = quadprog(H, f, -Ac, -b, [], [], lb, ub, [], options);
    
    if isempty(ck)
        c = 0;
    else
    c = ck(1);
    end
    
    varW = 0.01;
    varV = 0.001;
    
    w =  varW*randn(no_states, 1);
    w(3) = w(3)*0.1 + varV*rand(1, 1) - 0.5*varV;
    
    v =  varV*rand(no_outputs, 1);
    v(2) = v(2)*0.1; 
    
    x(:, k+1) = A*x(:, k) + B*c + w;
    y(:, k) = C*x(:, k) + v;
    
    X = x(:, k);
    
    x2(:, k+1) = A*x2(:, k) + w;
    y2(:, k) = C*x2(:, k) + v;
    
    
    Ck(k) = c;
end

%%
figure
plot(y(1, :))
hold on
plot(y2(1, :));
grid on
plot([0, Time_out/Ts], [main_bounds(1), main_bounds(1)], 'k')
plot([0, Time_out/Ts], [-main_bounds(1), -main_bounds(1)], 'k')

plot(y(2, :))
hold on
plot(y2(2, :));

figure
plot(y(2, :))
hold on
plot(y2(2, :));

grid on

figure
stairs(Ck)

grid on
