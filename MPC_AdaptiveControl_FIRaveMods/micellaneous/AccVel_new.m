% Second order velocity acceleration model 

Ts = 0.001;
runT = 30;

% A = [1 0; 
%     Ts 1];
% B = [1;
%      0];
A = 2*rand(2,2) - 1; B = 2*rand(2,1) - 1;

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
    
    y(:, k) = C*x(:, k);    
    u(:, k) = uk;
    x(:, k+1) = A*x(:, k) + B*uk + 1e-3*randn(2,1);
    
end

n = 100;
m = 2;

[Y, D, P, P_expanded] = least_squares_params_new(y(:, 1000:2000), u(:, 1000:2000), n, m);

% data = iddata(y(:, 1000:2000)', [u(:, 1000:2000); u(:, 1000:2000)]', Ts);
% sys = armax(data,[2 3 2 1]);

Dcgain = (eye(length(A)) - A)\B;
Dcgain = C*Dcgain;

% a_sum = sum(P(1:m, 1), 1);
% bx_sum = sum(P(m+1: 2*m+1), 1);
% bp_sum = sum(P(m+2:end, 1), 1);
% Dcgain2 = [bx_sum/a_sum; bp_sum/a_sum];

p.a = [1,-P(1:m)'];
p.b = reshape(P(m+1:end),m,size(C,1))';
Dcgain2 = sum(p.b')'/sum(p.a);

YY = reshape(Y,2,n);
YYp = reshape(D*P,2,n);
e = YY - YYp;
mse = sum((e.^2)')'/n;

%% Plots
figure;
plot(y');
figure;
plot([YY;YYp]');

sys = ss(A,B,C,0,-1); stf = tf(sys);
[p.a; stf.Denominator{1}], %#ok<*NOPTS>
[[0, p.b(1,:)]; stf.Numerator{1}],
[[0, p.b(2,:)]; stf.Numerator{2}],
