% Second order velocity acceleration model 

Ts = 0.001;
runT = 30;

A = [1 0; 
    Ts 1];
B = [1;
     0];
 
C = eye(length(A));

[K,~,~] = dlqr(A,B,C'*C,1,0);

A = A - B*K;

x = zeros(length(A), runT/Ts);
x(1, 1) = 0;
y = zeros(2, runT/Ts);
u = zeros(1, runT/Ts);

for k = 1: length(x)-1
    
    if k < 5000
        uk = 0.1*randn(1, 1);
    
    else
        uk = 0;
    end  
    
    x(:, k+1) = A*x(:, k) + B*uk + 0.1*randn(2,1);
    y(:, k) = C*x(:, k);
    
    u(:, k) = uk;
    
end

n = 100;
m = 2;

[Y, D, P, P_expanded] = least_squares_params(y(:, 1000:2000), u(:, 1000:2000), n, m);

% data = iddata(y(:, 1000:2000)', [u(:, 1000:2000); u(:, 1000:2000)]', Ts);
% sys = armax(data,[2 3 2 1]);

Dcgain = (eye(length(A)) - A)\B;
Dcgain = C*Dcgain;

a_sum = sum(P(1:m, 1), 1);
bx_sum = sum(P(m+1: 2*m+1), 1);
bp_sum = sum(P(m+2:end, 1), 1);

Dcgain2 = [bx_sum/a_sum; bp_sum/a_sum];

%% Plots
figure
plot(y(1, :))
figure
plot(y(2, :))
