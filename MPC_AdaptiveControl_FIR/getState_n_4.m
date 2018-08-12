function X = getState_n_4(y, Ck, PhiP, Bp, Cp)

M1 = [Cp*PhiP; Cp*PhiP^2];

G1 = [Cp*Bp, zeros(2, 1);    
      Cp*PhiP*Bp, Cp*Bp]; 


yk = y(:, end-1: end);
  
Y = reshape(yk, [], 1);  
C_in = fliplr(Ck(:, end-1: end));

C_in = C_in';
  
X = M1 \ (Y - G1*C_in);  

end