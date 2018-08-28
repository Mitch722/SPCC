A = [0 1; -2 -3];

% define unknown matrix to be found using LMI
setlmis([])
P = lmivar(1, [size(A, 1) 1]);  % structure and size of LMI

% def LMI
lmiterm([1 1 1 P], 1 , A, 's');    % A'P+PA
lmiterm([1 1 2 0], 1);             % 0
lmiterm([1 2 2 P], -1, 1);         % P > 0

LMISYS = getlmis;

[tmin, Psol] = feasp(LMISYS);
P = dec2mat(LMISYS, Psol, P);
