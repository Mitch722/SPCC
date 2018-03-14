% Run simulation with varying time Dynamics

TsFast = 0.001;
Time_out = 10;

TsObvs = 0.01;
% main_bounds = [x, phi, u]
main_bounds = [0.8, 0.15, 0.4]';

ratioTs = TsObvs / TsFast;
%% Initialise the variable masses

M0 = 1.5;   M = M0;
m0 = 0.2;   m = m0;

x = zeros(4, Time_out/TsFast);
y = zeros(2, Time_out/TsFast);

xhat = zeros(4, Time_out/TsObvs);
yhat = zeros(2, Time_out/TsObvs);

[sys_obv, sysdObv, L, K_opt] = inverted_pen_T(TsObvs, M0, m0);
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
        
        xhat(:, k0 + 1) = (sysdObv.A - sysdObv.B*K_opt - L*sysdObv.C)*xhat(:, k0) + L*y(:, k) + w;
        yhat(:, k0 + 1) = sysdObv.C*xhat(:, k0) + v;
        
    end
    % Physics of dyanamics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = M - 0.001*M;
    m = m - 0.0001*m;
    
    varW = 0.01;
    varV = 0.01;
    w =  varW*randn(4, 1);
    w(3) = w(3)*0.1 + varV*rand(1, 1) - 0.5*varV;
    v =  varV*rand(2, 1);
    v(2) = v(2)*0.1; 
    
    sysd = inverted_pen_T_model(TsFast, M, m);
    
    ck = 0;
    
    uk = -K_opt*x(:, k) + ck; 
    x(:, k+1) = sysd.A*x(:, k) + sysd.B*uk + w;
    y(:, k+1) = sysd.C*x(:, k+1) + v;
    
    
end


%%
figure
plot(y(1, :))
hold on
% plot(y2(1, :));
grid on

stairs(main_bounds(1), 'k')
stairs(-main_bounds(1), 'k')

title('Cart Position MPC vs LQR')
xlabel('Time Steps')
ylabel('Cart Position from Centre')

figure
plot(y(2, :))
hold on
% plot(y2(2, :));

% stairs(main_bounds(2), 'k')
% stairs(-main_bounds(2), 'k')

grid on
title('Angle of Pendulum phi')
xlabel('Time Steps')
ylabel('Angle phi of Pendulum')

% figure
% stairs(u)

% grid on
% title('Reference Input MPC')