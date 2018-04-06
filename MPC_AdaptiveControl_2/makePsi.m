function Psi = makePsi(A, B, p, i)

BeN = zeros(length(A), p);

[~, colB] = size(B);

for counter = 1: i
     
   BeN(:, (counter - 1)*colB +1: counter*colB) = A^(i-counter)*B;
    
end

Psi = [A^i,                 BeN;
       zeros(p, length(A)), eye(p)];
end

