

[sys_obv, L, K_opt] = inverted_pen;

A = sys_obv.A;
B = sys_obv.B;
C = sys_obv.C;
D = sys_obv.D;
Ts = sys_obv.Ts;

[~, no_states] = size(A);
[no_outputs, ~] = size(D);

T_end = 20;

x = zeros(no_states, T_end/Ts);
y = zeros(no_outputs, length(x));
u = zeros(1, length(x));

for k = 1: length(x) - 1
   
    u(1, k) = randn(1,1);
    
    x(:, k+1) = A*x(:, k) + B*u(1, k);
    
    y(:,k) = C*x(:, k);
    
end

%% Plots

figure
plot(y(1, :))
figure
plot(y(2, :))

%% Parameter estimation
m = 4;

[Y, D, P, P_expanded] = least_squares_params(y(:, 1:1000), u(:, 1:1000), 100, m);

Dcgain = (eye(length(A)) - A)\B;
Dcgain = C*Dcgain;

a_sum = sum(P(1:m, 1), 1);
bx_sum = sum(P(m+1: 2*m+1), 1);
bp_sum = sum(P(m+2:end, 1), 1);

Dcgain2 = [bx_sum/a_sum; bp_sum/a_sum];
%%

[Ap, Bp, Cp] = estSS(P, m);

Q_bar = dlyap(Ap,eye(m));
