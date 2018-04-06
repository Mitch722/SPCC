% Unit test least_squares_params

% Generates an output as 1: 50 over 51: 100
y = [1 : 50;
     51 : 100];
% takes 1:50 as the input sequence 
c = y(1, :);

% runs the code to test that Y = D*P is being genereated properly
[Y, D, P, P_expanded] = least_squares_params(y, c, 3, 2);
