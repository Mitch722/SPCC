function [y, u, t1, yhat, main_bounds] = AdaptiveMPCsimFIR(bias, Time_out, nWidth, s)

% Run simulation with varying time Dynamics

TsFast = 0.005;

TsObvs = 0.01;

rng(s);
%% Define the Observer
M0 = 1.5;   M = M0;
m0 = 0.2;   m = m0;

[~, sysdObv, L, K_opt] = inverted_pen_T(TsObvs, M0, m0);

A = sysdObv.A;
B = sysdObv.B;
C = sysdObv.C;
D = sysdObv.D;

p = 20;

Q = C'*C;
R = 1;

% main_bounds = [x, phi, u]
main_bounds = [0.8, 0.15, 40]';

ratioTs = TsObvs / TsFast;
%% Initialise the variable masses
x = zeros(4, Time_out/TsFast);
y = zeros(2, Time_out/TsFast);

xhat = zeros(4, Time_out/TsObvs);
yhat = zeros(2, Time_out/TsObvs);

u = x(1, :);
u_adapt = u;
u_lq = u;

Ck = zeros(1, Time_out/TsFast);
c = 0;

maxF = 100;

[H, f, Ac, Ax, b1, lb, ub, opt] = MPC_vars(A-B*K_opt, B, C, K_opt, R, p, main_bounds, maxF);
% cholesky for mpcqpsolver
[L2, ~] = chol(H,'lower');
Linv = inv(L2);

% options for mpcqpsolver:
options = mpcqpsolverOptions;

adaptTime = 1500;
%%
k0 = 1;
for k = 1 : Time_out/TsFast
        
    % Observer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if k / ratioTs == round(k/ratioTs)
        
        k0 = k/ ratioTs;
        
        varW = 0.0051;
        varV = 0.0051;
        w =  varW*randn(4, 1);
        w(3) = w(3)*0.1 + varV*rand(1, 1) - 0.5*varV;
        v =  varV*rand(2, 1);
        v(2) = v(2)*0.1; 
        
        xhat(:, k0 + 1) = (sysdObv.A - sysdObv.B*K_opt - L*sysdObv.C)*xhat(:, k0) + sysdObv.B*c + L*y(:, k) + w;
        yhat(:, k0 + 1) = sysdObv.C*xhat(:, k0) + v;
        
        X = xhat(:, k0);
        
        b = b1 + Ax*X;
        
        if k <= adaptTime - 1
            % [x,status] = mpcqpsolver(Linv,f,A,b,Aeq,beq,iA0,options)
            ck = mpcqpsolver(Linv, f, Ac, b, [], zeros(0,1), false(size(b)), options);
            % ck = quadprog(H, f, -Ac, -b, [], [], lb, ub, [], options);

            if isempty(ck) || abs(ck(1)) > 10
                c = 0;
            else
                c = ck(1);
            end
        end
    end       
    % Initialise Adaptive Control
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % bounds
    bnds = [0.8, 1.5]';
    k0 = round(k/ ratioTs);    % controller step
    X = xhat(:, k0);    % current state estimate
    % parameters for FIR models
    params.n = 10;       % no. of parameters
    params.m = params.n+5;    % no. rows 
    M_samps = 55;       % no. of models in multi-sample
    
    if k == adaptTime - 1
        % Try statement 
        [~, ~, RstarModel, entry] = ACSA_FIR(M_samps, y(:, 1:k-1), u(:, 1:k-1), p, R, params, H, bnds, A-B*K_opt, B, C, K_opt, X);      
        
    end
    % Adaptive Control run at k0 intervals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if k > adaptTime  && k / ratioTs == round(k/ratioTs)
        % feed in the previous state
        Xprev = xhat(:, k0-params.n);
        if k/25 == round(k/25)
            % selects new constraints
            [~, ~, RstarModel, entry] = ACSA_FIR(M_samps, y(:, 1:k-1), u(:, 1:k-1), p, R, params, H, bnds, A-B*K_opt, B, C, K_opt, Xprev);
                       
        end   
        % average the FIR big matrices
        % sum them first 
        % concatonate the matrices into a 3d matrix
        threedFfir = cat(3, RstarModel{entry.bigFIR, :});
        % add up the matrices in the 
        aveFir = sum(threedFfir, 3)/length(threedFfir(1,1,:));
        aveFir = fliplr(aveFir);
        
        [H2, f, Dc2, Dx2, b2, lb, ub, op] = MPC_varsFIR(A-B*K_opt, B, C, K_opt, R, p, bnds, maxF, aveFir);
        [L21, ~] = chol(H2,'lower');
        Linv1 = inv(L21);
        
        % X = xhat(:, k0);
        b = b2 + Dx2*Xprev;
                
        ck = mpcqpsolver(Linv1, f, Dc2, b, [], zeros(0,1), false(size(b)), options);
        if isempty(ck) || abs(ck(1)) > 10
            c = 0;
        else
            c = ck(1);
        end
    end    
    % Physical Model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = M - 0.01*TsFast*M + 0.001*TsFast*randn(1,1);
    m = m - 0.01*TsFast*m + 0.001*TsFast*randn(1,1);
    
    varW = 0.0051;
    varV = 0.0051;
    w =  varW*randn(4, 1);
    w(3) = w(3)*0.1 + varV*rand(1, 1) - 0.5*varV;
    v =  varV*rand(2, 1);
    v(2) = v(2)*0.1; 
    
    sysd = inverted_pen_T_model(TsFast, M, m);
    % use the MPC output value
    Ck(k) = c;
    
    uk = -K_opt*xhat(:, k0) + Ck(k); 
    
    % limit on uk
    if abs(uk) > 100
        uk = 100 * abs(uk)/uk;
    end
    
    u(k) = uk;
    
    u_lq(k) = -K_opt*xhat(:, k0);
    
    Ad = sysd.A + nWidth*rand(4) - bias*nWidth*ones(4);
        
    x(:, k+1) = Ad*x(:, k) + sysd.B*uk + w;
    y(:, k+1) = sysd.C*x(:, k+1) + v;
     
end


%%
figure

t1 = linspace(0, Time_out, length(y(1, :)));
t2 = linspace(0, Time_out, length(yhat(1, :)));

plot(t1, y(1, :), 'b')
hold on
% plot(t2, yhat(1, :));
grid on

stairs([0, Time_out], [main_bounds(1), main_bounds(1)], 'k')
stairs([0, Time_out], [-main_bounds(1), -main_bounds(1)], 'k')

title('Cart Position MPC vs LQR')
xlabel('Time/s')
ylabel('Cart Position from Centre')

figure
plot(t1, y(2, :), 'b')
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
end