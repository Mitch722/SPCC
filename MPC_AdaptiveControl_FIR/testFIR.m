

uk = 1:10;
Yk = [110:1:200; 410:1:500];

n_c = 2;
m_c = 3;
% build FIR coefficients
F = buildFIR(Yk, uk, n_c, m_c);
%% test FIR coeffs on a 2nd order random problem

A = 0.1*randn(2);
B = [0; 1];

C = eye(2);

n = 100;
x = zeros(2, n);
y = zeros(2, n);

u = 0.001*randn(1, n);

for k = 1:n-1
    
   x(:, k+1) = A*x(:, k) + B*u(:, k);
   y(:, k) = C*x(:, k);
    
end
%% FIR

Fy = buildFIR(y(1,:), u, 20, 50);
Fp = buildFIR(y(2,:), u, 20, 50);

Fcoefs = [Fy'; Fp'];

p = 21;
F_big = ConvFIRmat(Fy', p);

yFIR = F_big*u(end-length(F_big(1,:))+1 :end )';

%% Test the DC gain
% DC gain FIR

dcFir = sum(Fy);
dcSS = C*(eye(2) - A)\B;

%%
time = linspace(1,n, length(y(1, :)));
plot(time, y(1, :), 'b')
hold on
plot(time, y(2, :), 'r')

figure
plot(time(1:p), yFIR')
hold on
plot(time(1:p), y(1, end-p+1 :end ))
