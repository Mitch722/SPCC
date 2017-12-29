function [r_star, Ptrail, Ntrail, q_min, q_max] = algo1(x, ep_lo, ep_hi, Ppor, Ppst)

% input args:
%
%   x: a struct containing 2 fields:
%       x.sample : the data contained in x
%       x.M : the size of the data
%       x.dim : the dimension of x
%    
% returns:
%   r_star
%   Ptrial
%   Ntrail
%   

zeta = [2, x.dim + 1];

% Find q underbar and overbar
% Ppst

ep_hat = ep_lo + 0.5*ep_del;

% Ppst = binocdf(qhat, x.M, 1 - ep_hat);
% inputted value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = 1 : x.M;

probmin = binocdf(q-zeta(2), x.M, 1 - ep_hi);
minQ = probmin >= 0.5*( 1 + Ppst);

q_min = min(q(minQ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

probmax = binocdf(q-zeta(1), x.M, 1 - ep_lo);

maxQ = probmin <= 0.5*( 1 - Ppst);
q_max = max(q(maxQ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arrQ = q_max:q_min;

n = 1;
zed = zeta(n);

r = zed : min(arrQ)/10;
r(end) = [];

RQ = zeros(length(r), length(arrQ));

for i = 1:length(r)
        
    A = log(x.M - r(i) + 1);
    coef = betaln(x.M - arrQ + 1, arrQ - r(i) + 1);
    B = betaln(x.M - arrQ + zed, arrQ - zed + 1);
    C = betaln(zed, r(i) - zed + 1);    
    
    logF = -A - coef + B - C;
    
    F = exp(logF);
    
    RQ(i, :) = F;

end

sumQ = sum(RQ, 2);

[~, index] = max(sumQ);
r_star = r(index);

fprintf('r* = %d \n', r_star)

%% P trial

Ptrail = zeros(1, length(arrQ));

A = log(x.M - r_star + 1);
coef = betaln(x.M - arrQ + 1, arrQ - r_star + 1);
B = betaln(x.M - arrQ + zed, arrQ - zed + 1);
C = betaln(zed, r_star - zed + 1);    

logF = -A - coef + B - C;

F = exp(logF);

Ptrail(1, :) = F;
Ptrail = sum(Ptrail, 2);



%% Ntrial

Ntrail = log(1 - (Ppor / Ppst)) / log(1 - Ptrail);

