function [sys_obv, L, K_opt] = inverted_pen_T(Ts, M, m)

%% Make the Inverted Pendulum model
% Physical data stored in a struct
% m = 0.2
data.m = m;
% M = 1.5
data.M = M;
data.I = 0.005;
data.l = 0.5;
data.b = 0.1;
data.g = 9.81;

states = {'x', 'x_dot', 'phi', 'phi_dot'};
output = {'x', 'phi'};
input = {'u'};

[sys, sysd] = pendulum_ss(data, Ts);

%% Controllability

Q_c = ctrb(sysd);
Q_b = obsv(sysd);

CntrlBT = rank(Q_c);
ObsrvBT = rank(Q_b);

assert(CntrlBT == length(sysd.C), 'The model is not Controllable')
assert(ObsrvBT == length(sysd.C), 'The model is not Observable')

%% LQRD

Q = sysd.C'*sys.C;
R = 1;
K_opt = dlqr(sysd.A,sysd.B,Q,R);

Ak = sysd.A - sysd.B*K_opt;
Bk = sysd.B;
Ck = sysd.C;
Dk = sysd.D;

sys_lqr = ss(Ak, Bk,Ck, Dk, Ts, 'statename',states,'inputname',input,'outputname',output);

t = 0 : Ts : 10;
r = 0.2*ones(size(t));

% [y,t,x] = lsim(sys_lqr, r, t);

poles_A = eig(Ak);
%% Find poles 
sysc_lqr = d2c(sys_lqr,'zoh');

polesCA = eig(sysc_lqr.A);

%% Discrete Observer Implementation
% make observer poles 10x faster than LQR poles
multiple = 0.1;

% Put poles near orgin of unit circle in z-domain
poles_lqr = [0.0001, 0.001, 0.002, 0.0002]';
obvs_poles = multiple * real(poles_A);

% make sure poles are negative
if obvs_poles > 0
    obvs_poles = -1 * obvs_poles;
end

if obvs_poles == 0
    error('Poles are at zero')
end

L = place(sys_lqr.A', sys_lqr.C', obvs_poles);
L = L';

% Built the Observer model

Aobv = [sysd.A - sysd.B*K_opt, sysd.B*K_opt;
        zeros(size(sys.A)),    sysd.A - L*Ck];
    
Bobv = [ sysd.B;
         zeros(size( sysd.B ))];
         
Cobv = [Ck, zeros(size(Ck))];

Dobv = [0; 0];

states = {'x' 'x_dot' 'phi' 'phi_dot' 'e1' 'e2' 'e3' 'e4'};
input = {'r'};
output = {'x'; 'phi'};

sys_obv = ss(Aobv,Bobv,Cobv,Dobv,Ts,'statename',states,'inputname',input,'outputname',output);
    
% Step Input and Plot
end