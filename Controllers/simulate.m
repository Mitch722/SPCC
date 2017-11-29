

function [x, y] = simulate(sysd, u, t, v, w)

% Works out the states and output of a digital state space model
% Args: sysd, the ss(A,B,C,D) of the control model in digital time
%       u, the input, a vector same length as t and rows as states

A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

[~, no_states] = size(A);
[rows_d, ~] = size(D);

% Preallocate data
x = zeros( no_states, length(t)+1 );

y = zeros( rows_d, length(t)+1 );

if nargin == 3
    v = zeros(length(x));
    w = zeros(length(y));
    
end

% model_states = no_states/2;

% noise = 0.01*randn(model_states, length(t));
% noise = [zeros(model_states, length(t)); noise];

for i = 1: length(t)
            
    x(:,i+1) = A*x(:,i) + B*u(:,i) + v(:, i);
    y(:,i) = C*x(:,i) + w(:, i);
    
end

x(:,end) = [];
y(:,end) = [];

