
rng('default')

TsFast = 0.005;

M = zeros(1, 30/TsFast);
m = M;

M0 = 1.5;   M(1) = M0;
m0 = 0.2;   m(1) = m0;

for k = 1 : 30/TsFast-1
    M(k+1) = M(k) - 0.01*TsFast*M(k) + 0.00*TsFast*randn(1,1);
    m(k+1) = m(k) - 0.01*TsFast*m(k) + 0.00*TsFast*randn(1,1);
end
%%

figure
plot(M, 'b');
grid on
figure
plot(m, 'b')
grid on
