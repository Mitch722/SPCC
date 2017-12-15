%% Find H2

Q_bar_2 = Q_bar(end-p:end, end-p:end);
H2 = zeros((N+1)*length(Q_bar_2));

[qx2, qy2] = size(Q_bar_2);

H2(1:qy2, 1:qx2) = Q_bar_2;
for i = 1: N-1
    
    M_trans = M';
    H2(i*qy2 + 1 : (i+1)*qy2, i*qx2 + 1:(i+1)*qx2) = (M_trans^i)*Q_bar_2*M^i;
    
end
H2(N*qy2 +1: end, N*qx2 +1 : end) = P_bar(end-p:end, end-p:end);
