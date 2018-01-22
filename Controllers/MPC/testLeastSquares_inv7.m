% test least squares

try 
    
    ey1.sample = x(7, 1300:end);
    ey1.dim = 2;
    
catch
    
    MPC_inv_7
    ey1.sample = x(7, 1300:end);
    ey1.dim = 2;
      
end

samp_len = 300;

[epsilon, probab] = LS_algo1(ey1, samp_len);
fprintf('epsilon = %f, confidence = %f \n', epsilon, probab)

eps = 0:0.001:1;
m = length(ey1.sample) - samp_len - 1;

q_eps = m*(1 - epsilon);

confidence = binocdf(q_eps - 2, m, 1 - eps);

figure
plot(eps, confidence)
grid on
