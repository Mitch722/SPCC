function val = bigbino(x, n, p)

a = log(n + 1);

b = betaln(n - x + 1, x + 1);
c = x*log(p);
d = (n - x)*log(1 - p);

logsum = 1 - a - b + c + d;

val = exp(logsum);

