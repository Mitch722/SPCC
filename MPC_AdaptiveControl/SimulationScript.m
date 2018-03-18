%% Runs Adaptive MPC and MPC 
tic
% Bias on the Model variance
bias = 0.53;
% Time out time
Time_out = 20;

% Run Adaptive Control algorithm
[yAdapt, uAdapt, ~, yhatAdapt, ~] = AdaptiveMPCsim(bias, Time_out);

% Run MPC algorithm
[yMPC, uMPC, t1, yhatMPC, main_bounds] = MPCsim(bias, Time_out);

%% Plot figures 

figure
plot(t1, yAdapt(1, :), 'b')
hold on
plot(t1, yMPC(1, :), 'r')
grid on

stairs([0, Time_out], [main_bounds(1), main_bounds(1)], 'k')
stairs([0, Time_out], [-main_bounds(1), -main_bounds(1)], 'k')

title('Cart Position Adaptive vs MPC')
xlabel('Time/s')
ylabel('Cart Position from Centre')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(t1, yAdapt(2, :), 'b')
hold on
plot(t1, yMPC(2, :), 'r')
grid on
title('Angle phi Adaptive vs MPC')
xlabel('Time/s')
ylabel('Angle phi of Pendulum')

stairs([0, Time_out], [main_bounds(2), main_bounds(2)], 'k')
stairs([0, Time_out], [-main_bounds(2), -main_bounds(2)], 'k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subplot of each individual output
figure
subplot(2,2,1)
plot(t1, yAdapt(1, :), 'b')
grid on
hold on
stairs([0, Time_out], [main_bounds(1), main_bounds(1)], 'k')
stairs([0, Time_out], [-main_bounds(1), -main_bounds(1)], 'k')
title('Cart Position Adaptive MPC')
xlabel('Time/s')
ylabel('Cart Position from Centre')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2)
plot(t1, yMPC(1, :), 'r')
grid on
hold on
stairs([0, Time_out], [main_bounds(1), main_bounds(1)], 'k')
stairs([0, Time_out], [-main_bounds(1), -main_bounds(1)], 'k')
title('Cart Position MPC')
xlabel('Time/s')
ylabel('Cart Position from Centre')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)
plot(t1, yAdapt(2, :), 'b')
grid on
hold on
stairs([0, Time_out], [main_bounds(2), main_bounds(2)], 'k')
stairs([0, Time_out], [-main_bounds(2), -main_bounds(2)], 'k')
title('Angle phi Adaptive')
xlabel('Time/s')
ylabel('Angle phi of Pendulum')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4)
plot(t1, yMPC(2, :), 'r')
grid on
hold on
stairs([0, Time_out], [main_bounds(2), main_bounds(2)], 'k')
stairs([0, Time_out], [-main_bounds(2), -main_bounds(2)], 'k')
title('Angle phi MPC')
xlabel('Time/s')
ylabel('Angle phi of Pendulum')

toc
