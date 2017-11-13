
[sys, sysd, Ts] = double_integ_0;

A = sysd.A;
B = sysd.B;
C = sysd.C;

L = [1, 1];

n = 2000;
x = zeros(2, n);
y = x;

u = 1;

for i = 2 : n
    
    w = betarnd(2,5);
    x(:,i) = A*x(:, i-1 ) + B*u;
    
    y(1,i) = C * x(:,i);
    
end
