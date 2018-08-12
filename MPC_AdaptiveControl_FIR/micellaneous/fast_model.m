function x = fast_model(x, u, Ts)


[sys_obv, L, K_opt] = inverted_pen_T(Ts);

A = sys_obv.A;
B = sys_obv.B;
C = sys_obv.C;
D = sys_obv.D;

[~, no_states] = size(A);
[no_outputs, ~] = size(D);

x = A*x + B*u;
y = C*x;

end