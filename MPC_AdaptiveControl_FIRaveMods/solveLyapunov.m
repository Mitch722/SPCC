function [Q_bar, psi] = solveLyapunov(PhiP, Bp, Q, K_opt, R, p)

if p > 1
    
    eN = eye(p);
    eN = eN(1, :);
    
    M = eye(p);

    
    psi = [PhiP,                   Bp*eN;
          zeros(p, length(PhiP)),   M];
elseif p == 1
        
    eN = eye(p);
    eN = eN(1, :);
    
    psi = [PhiP,                   Bp*eN];
    psi = [psi; eye(size(psi))];
         
end

Q2 = [Q,                zeros(length(Q), p);
      zeros(p, length(Q)),  R*eye(p)];
  
        
% Q_bar is the solution to the discrete Lyapunov equation
Q_bar = Q2(end-p+1 : end, end-p+1 : end);

[~,d] = chol(Q_bar);
assert(d<=1, 'Q_bar is not +ve semi definite')

end