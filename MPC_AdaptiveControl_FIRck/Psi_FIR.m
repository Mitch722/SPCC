function [Ds] = vecUk(A, B, K_opt, p)
% 

%% begin solution
% A is A + B*K_opt
phi = A;
[~, no_states] = size(A);

K_opt = [K_opt, zeros(1, no_states - length(K_opt))];

M = [zeros(p-1, 1), eye(p-1);
     zeros(1, p)];    

eN = eye(p);
eN = eN(1, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make psi
psi = [phi,                     B*eN;
       zeros(p, length(phi)),   M];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
E = [1, zeros(1, p-1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% update F such that only gives output once multiplied with conv matrix
F = [K_opt, E];
 
% F = [C, zeros(len_output, p)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% initialise As
Ds = F;
for i = 1 : p-1
    % inefficient but easy to make
    Ds = [Ds ; F*psi^i];
    % A( (i)*ax + 1 : (i + 1)*ax, :) = F*psi^i; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end