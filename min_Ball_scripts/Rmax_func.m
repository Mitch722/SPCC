
function [Theo37, Qmax] = Rmax_func(arrBQ, r, x, zed)


Theo37 = zeros(size(arrBQ));

for i = 1: length(arrBQ)
    
    % r = arrBQ(i);
    
    n = x.M - r;
    
    k = arrBQ(i) - r;
    beta_coefln =  betaln(n - k + 1, k+1);
    beta_coefln = log(1 / (n+1)) - beta_coefln  ;
    
    a = x.M - arrBQ(i) + zed;
    b = arrBQ(i) - zed + 1;
    
    c = zed;
    d = r - zed + 1;
    
    if d < 0
        d = 1;
    end
    
    bulk = betaln(a,b) - betaln(c,d);
    assert(bulk <= eps, 'bulk less than or equal to eps')
    
    Theo37_inter = beta_coefln + bulk;
    Theo37(i) = exp(Theo37_inter);
    
end

[~, indexMax] =  max(Theo37);

Qmax = arrBQ(indexMax);

end