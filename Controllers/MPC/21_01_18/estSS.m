function [Ap, Bp, Cp] = estSS(P, m)

% take only the parameters for the states:
    
Ap = zeros(m);
Ap(1:m-1, 2:m) = eye(m-1);

last_row_A = fliplr(P(1:m)');

Ap(end, :) = (last_row_A);
% Ap(end-1, :) = last_row_A;

Bp = zeros(m, 1);
% Bp(end-1 : end, 1) = [1; 1];

Bp(end, 1) = 1;

Cp_inter = P(m+2: end, :)';
Cp_inter(m+1) = [];
% flip to put correct coefficients in the correct place
Cp_inter = fliplr(Cp_inter);
% Cp = zeros(2, m);

Cp = [Cp_inter(:, 1:m); Cp_inter(:, m+1:end)];

Cp2 = P(m+1, 1).*[last_row_A(1, 1:m)];
Cp2 = [Cp2; P(2*m+2, 1).*last_row_A(1, 1:m)];

Cp = Cp + Cp2;


end