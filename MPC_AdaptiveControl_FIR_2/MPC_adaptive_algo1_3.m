

Ts = 0.001;
M = 1.5;
m = 0.2;
% create inverted pendulum models
[sys_obv, sysd, L, K_opt] = inverted_pen_T(Ts, M, m);

A = sys_obv.A;
B = sys_obv.B;
C = sys_obv.C;
D = sys_obv.D;

p = 12;

Q = sys_obv.C'*sys_obv.C;
R = 1;
% main_bounds = [x, phi, u]
main_bounds = [0.8, 0.15, 0.4]';

[~, no_states] = size(sys_obv.A);
[no_outputs, ~] = size(sys_obv.D);


%% prepare MPC 
Time_out = 30;
x = zeros(no_states, Time_out/Ts);
y = zeros(no_outputs, Time_out/Ts);

u = x(1, :);
u_adapt = u;

X = x(:, 1);
maxF = 100;

[H, f, Ac, Ax, b1, lb, ub, opt] = MPC_vars(A, B, C, K_opt, R, p, main_bounds, maxF);
% cholesky for mpcqpsolver
[L2, ~] = chol(H,'lower');
Linv = inv(L2);

% options for mpcqpsolver:
options = mpcqpsolverOptions;

for k = 1: (Time_out/Ts)-1 
    
    % progress bar
    if k/100 == round(k/100)
       fprintf('#') 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if k < 1299
    % Run algo 1 to update bound    
    b = b1 + Ax*X;
        
    % [x,status] = mpcqpsolver(Linv,f,A,b,Aeq,beq,iA0,options)
    ck = mpcqpsolver(Linv, zeros(p, 1), Ac, b, [], zeros(0,1), false(size(b)), options);
    % ck = quadprog(H, f, -Ac, -b, [], [], lb, ub, [], options);
    
        if isempty(ck) || abs(ck(1)) > 10
            c = 0;
        else
            c = ck(1);
        end
     
        
    varW = 0.01;
    varV = 0.01;

    w =  varW*randn(no_states, 1);
    w(3) = w(3)*0.1 + varV*rand(1, 1) - 0.5*varV;

    v =  varV*rand(no_outputs, 1);
    v(2) = v(2)*0.1; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x(:, k+1) = sys_obv.A*x(:, k) + sys_obv.B*c + w;
    y(:, k) = sys_obv.C*x(:, k) + v;   
    
    X = x(:, k+1);
    u(1, k) = [K_opt, zeros(1, 4)]*X + c;
    
    % initialise the adaptive control/ input/ output models
        if k == 1298
            % Try statement 

            params.m = 4;
            params.n = 100 - params.m;
            params.Ts = Ts;

            [~, ~, P, ~] = least_squares_params(y(:, 1: k-1), u(:, 1:k-1), params.n, params.m);

            [Ap2, Bp2, Cp2] = estSS(P, params.m);

            Kp2 = dlqr(Ap2,Bp2,Cp2'*Cp2,1,0);
            
            PhiP2 = Ap2 - Bp2*Kp2;
            % generate the states for the parameter estimation
            Xp = getState_n_4(y(:, 1:k-1), u(:, 1:k-1), PhiP2, Bp2, Cp2);
            
            bnds = [0.8, 1.5, 1]';
            [Q_bar, ~, ~, ~, ~, ~, ~, ~] = MPC_vars(PhiP2, Bp2, Cp2, Kp2, R, p, bnds, maxF);
            
            % This would be threaded used to generate RStarModels
            M = 80;
        
            params.m = 4;
            params.n = 10;
            [ck, Mmodels, RstarModel, entry] = ACSA_1(M, y(:, 1:k-1), u(:, 1:k-1), p, params, Q_bar, bnds, Xp);
                        
        end
     
     else
        M = 80;
        
        params.m = 4;
        params.n = 10;

        
        if k/50 == round(k/50)
            
            Xp = getState_n_4(y(:, 1:k-1), u(:, 1:k-1), PhiP2, Bp2, Cp2);
            
            [~, Mmodels, RstarModel, entry] = ACSA_1(M, y(:, 1:k-1), u(:, 1:k-1), p, params, Q_bar, bnds, Xp);
        end
        % genereate the states from set RstarModel
        % generate uk = [Kopt1, ... Kopti, ... KoptR]*[Xp1,...XpR]' + ck
        [uk, ck] = RmodelOutput(Q_bar, RstarModel, entry, y(:, 1:k-1), u(:, 1:k-1), p);
        
        if abs(uk) > 4 || abs(ck(1)) > 4
            uk = 0;
            ck = 0;
        end
        
        c = ck(1);
        u_adapt(k) = uk;
        
        varW = 0.01;
        varV = 0.01;

        w =  varW*randn(4, 1);
        w(3) = w(3)*0.1 + varV*rand(1, 1) - 0.5*varV;

        v =  varV*rand(no_outputs, 1);
        v(2) = v(2)*0.1; 
        
        x(1:4, k+1) = sysd.A*x(1:4, k) + sysd.B*uk + w;
        y(:, k) = sysd.C*x(1:4, k) + v;
        
        u(1, k) = uk + 0.01*randn(1,1);
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
end

fprintf('\n')
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

stairs(main_bounds(2), 'k')
stairs(-main_bounds(2), 'k')

grid on
title('Angle of Pendulum phi')
xlabel('Time Steps')
ylabel('Angle phi of Pendulum')

figure
stairs(u)

grid on
title('Reference Input MPC')

toc


