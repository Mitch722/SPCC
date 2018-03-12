load('output_input.mat')

% This is to test if the Y and D matrices are being made correctly

%% Try the least squares 
n = 50;
m = 10;

y = zeros(2, 100);
c = zeros(1, 100);
for i = 1 : length(y)
   
    y(:, i) = [i; i];
    c(:, i) = i;
end

[Y, D, P, P_expanded] = least_squares_params(y, c, n, m);

Y_hat = D*P;

e = Y - Y_hat;

%% Functions

function [Y, D, P, P_expanded] = least_squares_params(y, c, n, m)
% This function takes the input data and makes an input output data model
assert(n > m, 'the number of rows should be greater than columns, necessary that n > m')
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