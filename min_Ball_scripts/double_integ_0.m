
function [sys, sysd, Ts] = double_integ_0

A = [0 0 ; 1 0];
B = [1; 0];

C = [1, 0];
D = 0;

sys = ss(A, B, C, D);

% bodeplot(sys)
% grid on

Ts = 0.1;
sysd = c2d(sys,Ts);

end


