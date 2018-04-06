% test adaptive control

load('output_input.mat')

y(:, end) = [];
Ck(:, end) = [];

%% Run the parameter estimation and check to see if it matches the actual value

p = 12;
Q_bar = eye(12);

k = 500;

params.m = 4;
params.n = 80 - params.m;
params.Ts = 0.01;

bnds = [0.8, 1.5, 0.4]';

[c_inter, PhiP, Bp, Cp, P] = adaptiveControl4(Q_bar, y(:, k-200: k-1), Ck(:, k-201:k-1), p, params, bnds); 

%% find M for the specific case to generate states

M1 = [Cp*PhiP; Cp*PhiP^2];

Minv = inv(M1);

G1 = [Cp*Bp, zeros(2, 1);    
      Cp*PhiP*Bp, Cp*Bp]; 


yk = y(:, k-params.m/2: k-1);
  
Y = reshape(yk, [], 1);  
C_in = fliplr(Ck(:, k-params.m/2: k-1));

C_in = C_in';
  
X = M1 \ (Y - G1*C_in);  

%%

X2 = getState_n_4(y(:, 1 : k-1), Ck(:, 1 : k-1), PhiP, Bp, Cp);

X2k = PhiP*X2 + Bp*Ck(k-1);
Y2k = Cp*X2k;

%% 

X_k1 = PhiP*X + Bp*Ck(k-1);
Y_k1 = Cp*X_k1;

Y_k_actual = y(:, k+1);

Y_k = Cp*X;

%% find the DC gains of the models

[sys_obv, L, K_opt] = inverted_pen;

A = sys_obv.A;
B = sys_obv.B;
C = sys_obv.C;
D = sys_obv.D;
Ts = sys_obv.Ts;

Dcgain = (eye(length(A)) - A)\B;
Dcgain = C*Dcgain;

a_sum = sum(P(1:params.m, 1), 1);
bx_sum = sum(P(params.m+1: 2*params.m+1), 1);
bp_sum = sum(P(params.m+2:end, 1), 1);

Dcgain2 = [bx_sum/a_sum; bp_sum/a_sum];

Dcgain3 = dcgain(sys_obv);

%% 
function X = getStates(Y, A, B, C, n)

no_states = length(A);

[no_outputs, ~] = size(C);

M = zeros(n*no_outputs, no_states);

for i = 1: n
   
    M( (i-1)*no_outputs: i*no_outputs  , :) = C*A^i;
    
end

G = zeros(n*no_outputs, n);

G(no_outputs,1) = C*B;

Co = fliplr(ctrb(A,B));


end


