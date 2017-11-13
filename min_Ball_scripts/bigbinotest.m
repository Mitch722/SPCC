

x = 1000;
n = 100000;
p = 0.1;

a = binopdf(x, n, p);
b = bigbino(x, n, p);

assert(a == b, 'It does not work dickhead')