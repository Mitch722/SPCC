% Run simulation with varying time Dynamics

TsFast = 0.005;
Time_out = 15;

TsObvs = 0.01;
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
u_adapt = u;

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
        
        varW = 0.01;
        varV = 0.01;
        w =  varW*randn(4, 1);
        w(3) = w(3)*0.1 + varV*rand(1, 1) - 0.5*varV;
        v =  varV*rand(2, 1);
        v(2) = v(2)*0.1; 
        
        xhat(:, k0 + 1) = (sysdObv.A - sysdObv.B*K_opt - L*sysdObv.C)*xhat(:, k0) + sysdObv.B*c + L*y(:, k) + w;
        yhat(:, k0 + 1) = sysdObv.C*xhat(:, k0) + v;
        
        X = xhat(:, k0);
        
        b = b1 + Ax*X;
        
        if k <= 1299
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
    if k == adaptTime - 1
        % Try statement 
        
        params.m = 4;
        params.n = 100 - params.m;
        % params.Ts = Ts;

        [~, ~, P, ~] = least_squares_params(y(:, 1: k-1), u(:, 1:k-1), params.n, params.m);

        [Ap2, Bp2, Cp2] = estSS(P, params.m);

        Kp2 = dlqr(Ap2,Bp2,Cp2'*Cp2,1,0);

        PhiP2 = Ap2 - Bp2*Kp2;
        % generate the states for the parameter estimation
        Xp = getState_n_4(y(:, 1:k-1), u(:, 1:k-1), PhiP2, Bp2, Cp2);

        bnds = [0.8, 1.5, 0.4]';
        [Q_bar, ~, ~, ~, ~, ~, ~, ~] = MPC_vars(PhiP2, Bp2, Cp2, Kp2, R, p, bnds, maxF);

        % This would be threaded used to generate RStarModels
        M_samps = 90;
        
        params.m = 4;
        params.n = 10;
        [~, Mmodels, RstarModel, entry] = ACSA_1(M_samps, y(:, 1:k-1), u(:, 1:k-1), p, params, Q_bar, bnds, Xp);

    end
    % Adaptive Control run at k0 intervals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bnds = [0.8, 1.5, 20]';
    if k > adaptTime && k / ratioTs == round(k/ratioTs)
        if k0/25 == round(k0/25)

            Xp = getState_n_4(y(:, 1:k-1), u(:, 1:k-1), PhiP2, Bp2, Cp2);

            [~, Mmodels, RstarModel, entry] = ACSA_1(M_samps, y(:, 1:k-1), u(:, 1:k-1), p, params, Q_bar, bnds, Xp);
            [noItems, noRstar] = size(RstarModel);
            % recalculate an average model
            PhiP2 = (1/noRstar)*sum(cat(3, RstarModel{entry.PhiP, :}), 3 );
            Bp2 = (1/noRstar)*sum(cat(3, RstarModel{entry.Bp, :}), 3 );
            Cp2 = (1/noRstar)*sum(cat(3, RstarModel{entry.Cp, :}), 3 );
            Kp2 = (1/noRstar)*sum(cat(3, RstarModel{entry.K_opt, :}), 3 );
            
        end
        % genereate the states from set RstarModel
        % generate uk = [Kopt1, ... Kopti, ... KoptR]*[Xp1,...XpR]' + ck
        % use an average model instead of all of the models
        aveModel = cell(3, 1);
        aveModel{entry.PhiP, 1} = PhiP2;
        aveModel{entry.Bp, 1} = Bp2;
        aveModel{entry.Cp, 1} = Cp2;
        [Q_bar, ~, Ac2, Ax2, b12, ~, ~, ~] = MPC_vars(PhiP2, Bp2, Cp2, Kp2, R, p, bnds, maxF);
        aveModel{entry.Ac, 1} = Ac2;
        aveModel{entry.Ax, 1} = Ax2;
        aveModel{entry.K_opt, 1} = Kp2;
        aveModel{entry.b1, 1} = b12;
                       
        [uka, ck] = RmodelOutput(Q_bar, aveModel, entry, y(:, 1:k-1), u(:, 1:k-1), p);

        if abs(ck(1)) > 1000
            ck = 0;
        end
        if abs(uka) > 10000
            uka = 0;
        end
        c = uka;
        u_adapt(k) = uka;
    end    
    % Physical Model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = M - 0.1*TsFast*M + 0.001*TsFast*randn(1,1);
    m = m - 0.1*TsFast*m + 0.001*TsFast*randn(1,1);
    
    varW = 0.01;
    varV = 0.01;
    w =  varW*randn(4, 1);
    w(3) = w(3)*0.1 + varV*rand(1, 1) - 0.5*varV;
    v =  varV*rand(2, 1);
    v(2) = v(2)*0.1; 
    
    sysd = inverted_pen_T_model(TsFast, M, m);
    % use the MPC output value
    Ck(k) = c;
    
    uk = -K_opt*xhat(:, k0) + Ck(k); 
    u(k) = uk;
    x(:, k+1) = (sysd.A + 0.05*rand(4) - 0.53*0.05*ones(4))*x(:, k) + sysd.B*uk + w;
    y(:, k+1) = sysd.C*x(:, k+1) + v;
     
end


%%
figure

t1 = linspace(0, Time_out, length(y(1, :)));
t2 = linspace(0, Time_out, length(yhat(1, :)));

plot(t1, y(1, :))
hold on
grid on

stairs([0, Time_out], [main_bounds(1), main_bounds(1)], 'k')
stairs([0, Time_out], [-main_bounds(1), -main_bounds(1)], 'k')

title('Cart Position MPC vs LQR')
xlabel('Time/s')
ylabel('Cart Position from Centre')

figure
plot(t1, y(2, :))
hold on

stairs([0, Time_out], [main_bounds(2), main_bounds(2)], 'k')
stairs([0, Time_out], [-main_bounds(2), -main_bounds(2)], 'k')

grid on
title('Angle of Pendulum phi')
xlabel('Time/s')
ylabel('Angle phi of Pendulum')

figure
stairs(Ck)

grid on
title('Reference Input MPC')