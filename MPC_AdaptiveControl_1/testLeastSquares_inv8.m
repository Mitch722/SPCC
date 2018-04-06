% test least squares

try 
    
    ey1.sample = x(5, :);
    ey1.dim = 2;
    
catch
    
    MPC_inv_8
    ey1.sample = x(5, :);
    ey1.dim = 2;
      
end

samp_len = 400;

[epsilon, probab] = LS_algo1(ey1, samp_len);
fprintf('epsilon = %f, confidence = %f \n', epsilon, probab)

eps = 0:0.001:1;
m = length(ey1.sample) - samp_len - 1;

q_eps = m*(1 - epsilon);

confidence = binocdf(q_eps - 1, m, 1 - eps);

figure
plot(eps, confidence)
grid on

hold on 

confidence = binocdf(q_eps - 100, m, 1 - eps);
plot(eps, confidence)
grid on

