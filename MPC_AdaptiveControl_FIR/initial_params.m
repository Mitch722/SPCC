function [PhiP, Bp, Cp, K_opt] = initial_params(Y, Ck, params)

n = params.n;
m = params.m;
Ts = params.Ts;

% takes data, generates model and finds optimal input c.

% perform least squares on complete output
[~, ~, P, ~] = least_squares_params(Y, Ck, n, m);

[PhiP, Bp, Cp] = estSS(P, m);

K_opt = zeros(1, length(PhiP));

%% Additional functions
function [Y, D, P, P_expanded] = least_squares_params(y, c, n, m)
% This function takes the input data and makes an input output data model
assert(n >= m, 'the number of rows should be greater than columns, necessary that n > m')
assert(length(y)>= n+m, 'y vector is too small: make n and m smaller'  )

% find y_0 index
y_0_index = length(y) - n - m;
% Y the column vector for data n long
Y = y(:, y_0_index + m + 1: end);
Y = reshape(Y, [], 1);
% make sure y or c is repeated

y2 = y;
c2 = c;

D = NaN*zeros(n*length(y(:,1)), 3*m);
% build up the D matrix
nm_end = y_0_index;

n_point = y_0_index + m - 1;

no_outputs2 = length(y2(:,1));

for i = 1:n
    
   d_inter = fliplr(y2(:, nm_end + i: n_point + i));
   
   c_input = fliplr(c2(:, nm_end + i: n_point + i)); 
   c_inter = [c_input, zeros(size(c_input))];
   
   c_inter = [c_inter; zeros(size(c_input)), c_input];
   
   d_inter = [d_inter, c_inter]; 
   
   D(i*no_outputs2 - no_outputs2 +1 : i*no_outputs2, :) = d_inter;
    
end

psudo = (D'*D);

P = psudo\D'*Y;

P_expanded = reshape(P, no_outputs2, []);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ap, Bp, Cp] = estSS(P, m)

% take only the parameters for the states:
m = m;
    
Ap = zeros(m);
Ap(2:m, 2:m) = eye(m-1);

last_row_A = P(1:m)';

Ap(end, :) = (last_row_A);

Bp = zeros(m, 1);
Bp(end-1 : end, 1) = [1; 1];

Cp_inter = P(m+1: end, :)';
% flip to put correct coefficients in the correct place
Cp_inter = fliplr(Cp_inter);
Cp = zeros(2, m);

for i = 1: m/2
    
    Cp(1, (i-1)*2 + 1 ) = Cp_inter(1, i);
    Cp(2,  2*i ) = Cp_inter(1, m+i);
    
end

end




end