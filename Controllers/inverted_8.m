%% Initialise physical model and parameters for controller

[sys_obv, L, K_opt] = inverted_pen;

A = sys_obv.A;
B = sys_obv.B;
C = sys_obv.C;
D = sys_obv.D;
Ts = sys_obv.Ts;

[~, no_states] = size(A);
[no_outputs, ~] = size(D);

t = 0: Ts: 20;

% Preallocate data
x = zeros( no_states, length(t)+1 );
y = zeros( no_outputs, length(t)+1 );

x_max = [2; 0.3166; 5; 0.644];
x_min = -x_max;

u_over_bar = 10;
u_under_bar = -10;

%% Find MPC optimal

for k = 1 : length(t)
   
    % work out length of Nc
    for N = 1 : length(t)
    
        u_max = K_opt*A^(N)*x_max;
        u_min = K_opt*A^(N)*x_min;
        
        if u_max <= u_over_bar && u_min >= u_under_bar
            Nc = N;
            break
        end
    end     
    
    M = zeros(no_states, 8*length(t));
    
    for i = 1 : Nc 
        
        
        
    end
    
end

