
function [sys, sys_d] = pendulum_ss(data, Ts)

% generate Discrete and Continuous State space models for an inverted
% pendulum
%
% Args: data, a struct containing physical data of model
%       Ts, the sampling period
% Returns: 
%       sys, the continuous time model
%       sys_d, the discrete time model

m = data.m;
M = data.M;
I = data.I;
l = data.l;
b = data.b;
g = data.g;

denom = ( I*(M + m) + M*m*l );

a22 = -(I +m*l^2)*b;
a22 = a22 / denom;

a23 = (g)*(m*l)^2;
a23 = a23 / denom;

a42 = -m*l*b;
a42 = a42 / denom;

a43 = m*g*l*(M+m);
a43 = a43 / denom;


b12 = (I + m*l^2);
b12 = b12 / denom;

b14 = m*l;
b14 = b14 / denom;


A = [0  1   0   0;
    
     0 a22 a23  0; 
     
     0   0   0  1;
     
     0  a42 a43 0];
 
B = [ 0;
     b12;
      0;
      b14];
  
C = [1 0 0 0;
     0 0 1 0];
 
D = [0; 0];

states = {'x', 'x_dot', 'phi', 'phi_dot'};
output = {'x', 'phi'};
input = {'u'};

sys = ss(A,B,C,D,'statename',states,'inputname',input,'outputname',output);

sys_d = c2d(sys,Ts,'zoh');

end

  
     

