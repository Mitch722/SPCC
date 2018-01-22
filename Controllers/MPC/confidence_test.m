

epsilon = 0:0.01:1;
m = 100;

q_eps = 99;

confidence = binocdf(q_eps - 1, m, 1 - epsilon);
plot(epsilon, confidence)
hold on

confidence = binocdf(q_eps - 2, m, 1 - epsilon);
plot(epsilon, confidence)