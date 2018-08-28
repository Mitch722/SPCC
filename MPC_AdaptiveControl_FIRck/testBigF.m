
F = randn(2, 3);
p = 4;

F_big = ConvFIRmat(F, p);

yhat = F_big*randn(6, 1);