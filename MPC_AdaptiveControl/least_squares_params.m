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

c2(end) = [];

D = NaN*zeros(n*length(y(:,1)), 3*m + 2);
% build up the D matrix
nm_end = y_0_index;

n_point = y_0_index + m - 1;

no_outputs2 = length(y2(:,1));

for i = 1:n
   
   d1 = nm_end + i;
   d2 = n_point + i;
   
   d_inter = fliplr(y2(:, d1: d2));
   
   c1p = nm_end + i - 1;
   c2p = n_point + i;
   
   c_input = fliplr(c2(:, c1p: c2p)); 
   c_inter = [c_input, zeros(size(c_input))];
   
   c_inter = [c_inter; zeros(size(c_input)), c_input];
   
   d_inter = [d_inter, c_inter]; 
   
   D(i*no_outputs2 - no_outputs2 +1 : i*no_outputs2, :) = d_inter;
    
end

psudo = (D'*D);

P = psudo\D'*Y;

P_expanded = reshape(P, no_outputs2, []);
end
