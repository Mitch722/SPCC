% Run simulation with varying time Dynamics

TsFast = 0.005;
Time_out = 20;

TsObvs = 0.01;

rng default
%% Define the Observer
M0 = 1.5;   M = M0;
m0 = 0.2;   m = m0;

[~, sysdObv, L, K_opt] = inverted_pen_T(TsObvs, M0, m0);

A = sysdObv.A;
B = sysdObv.B;
C = sysdObv.C;
D = sysdObv.D;

p = 24;

Q = C'*C;
R = 1;

% main_bounds = [x, phi, u]
main_bounds = [0.8, 0.15, 0.4]';

ratioTs = TsObvs / TsFast;
%% Initialise the variable masses
x = zeros(4, Time_out/TsFast);
y = zeros(2, Time_out/TsFast);

xhat = zeros(4, Time_out/TsObvs);
yhat = zeros(2, Time_out/TsObvs);

u = x(1, :);

Ck_sequence = zeros(1, Time_out/TsObvs);
c = 0;

maxF = 100;

[H, f, Ac, Ax, b1, lb, ub, opt] = MPC_vars(A-B*K_opt, B, C, K_opt, R, p, main_bounds, maxF);
% cholesky for mpcqpsolver
[L2, ~] = chol(H,'lower');
Linv = inv(L2);

% options for mpcqpsolver:
options = mpcqpsolverOptions;

%%
k0 = 1;
for k = 1 : Time_out/TsFast
        
    % Observer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if k / ratioTs == round(k/ratioTs)
        
        k0 = k/ ratioTs;
        
        varW = 0.01;
        varV = 0.01;
        w =  varW*randn(4, 1);
        w(3) = w(3)*0.1 + varV*rand(1, 1) - 0.5*varV;
        v =  varV*rand(2, 1);
        v(2) = v(2)*0.1; 
        % find the current States
        X = xhat(:, k0);
        b = b1 + Ax*X;
        
        % [x,status] = mpcqpsolver(Linv,f,A,b,Aeq,beq,iA0,options)
        ck = mpcqpsolver(Linv, f, Ac, b, [], zeros(0,1), false(size(b)), options);
        % ck = quadprog(H, f, -Ac, -b, [], [], lb, ub, [], options);

        if isempty(ck) || abs(ck(1)) > 10
            c = 0;
        else
            c = ck(1);
        end
        
        xhat(:, k0 + 1) = (sysdObv.A - sysdObv.B*K_opt - L*sysdObv.C)*xhat(:, k0) + sysdObv.B*c + L*y(:, k) + w;
        yhat(:, k0 + 1) = sysdObv.C*xhat(:, k0) + v;
             
    end       
    % Physics of dyanamics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = M - 0.01*TsFast*M + 0.001*TsFast*randn(1,1);
    m = m - 0.01*TsFast*m + 0.001*TsFast*randn(1,1);
    
    varW = 0.01;
    varV = 0.01;
    w =  varW*randn(4, 1);
    w(3) = w(3)*0.1 + varV*rand(1, 1) - 0.5*varV;
    v =  varV*rand(2, 1);
    v(2) = v(2)*0.1; 
    
    sysd = inverted_pen_T_model(TsFast, M, m);
    % use the MPC output value
    Ck = c;
    
    uk = -K_opt*xhat(:, k0) + Ck; 
    u(k) = uk;
    x(:, k+1) = (sysd.A + 0.05*rand(4) - 0.53*0.05*ones(4))*x(:, k) + sysd.B*uk + w;
    y(:, k+1) = sysd.C*x(:, k+1) + v;
     
end


%%
figure

t1 = linspace(0, Time_out, length(y(1, :)));
t2 = linspace(0, Time_out, length(yhat(1, :)));

plot(t1, y(1, :), 'r')
hold on
% plot(t2, yhat(1, :));
grid on

stairs([0, Time_out], [main_bounds(1), main_bounds(1)], 'k')
stairs([0, Time_out], [-main_bounds(1), -main_bounds(1)], 'k')

title('Cart Position MPC vs Adaptive MPC')
xlabel('Time/s')
ylabel('Cart Position from Centre')

figure
plot(t1, y(2, :), 'r')
hold on
% plot(t2, yhat(2, :));

stairs([0, Time_out], [main_bounds(2), main_bounds(2)], 'k')
stairs([0, Time_out], [-main_bounds(2), -main_bounds(2)], 'k')

grid on
title('Angle of Pendulum phi')
xlabel('Time/s')
ylabel('Angle phi of Pendulum')

figure
stairs(u)

 grid on
 title('Reference Input MPC')