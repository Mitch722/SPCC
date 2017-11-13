
A = [-2, 1;
     0, 1];
 
B = [1 ; 1];

C = [1 1];

Q = C'*C;

K = [2 -1];

R = 1;


Al = (A + B*K)';
Ql = (Q + K'*R*K);

Q_bar = dlyap(Al, Ql);
