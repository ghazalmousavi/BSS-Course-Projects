clc; clear; close all
warning('off', 'all');
load("hw9.mat")
rng(40)

%% Plot:
x = A*S;
x_noisy = A*S + Noise;

figure;
subplot(311)
plot(S(1, :), 'LineWidth', 1.2)
hold on 
plot(S(2, :), 'LineWidth', 1.2)
plot(S(3, :), 'LineWidth', 1.2)
legend('S_1', 'S_2', 'S_3')
title('Sources', 'Interpreter', 'latex')
hold off

subplot(312)
plot(x(1, :), 'LineWidth', 1.2)
hold on 
plot(x(2, :), 'LineWidth', 1.2)
plot(x(3, :), 'LineWidth', 1.2)
legend('X_1', 'X_2', 'X_3')
title('Observation', 'Interpreter', 'latex')
hold off

subplot(313)
plot(x_noisy(1, :), 'LineWidth', 1.2)
hold on 
plot(x_noisy(2, :), 'LineWidth', 1.2)
plot(x_noisy(3, :), 'LineWidth', 1.2)
legend('X_1', 'X_2', 'X_3')
title('Observation with noise', 'Interpreter', 'latex')
hold off
%% Define constants
[N, M] = size(x);
max_itr = 1e3;
lr = 0.01;

B = randn(N, N);
f = zeros(1, max_itr);

%% Q1: Optimization for learning matrix B:
for itr = 1 : max_itr

    y  = B*x_noisy;

    sai = zeros(size(x_noisy));
    for i = 1 : N
        sai(i, :) = mse_score_function(y(i, :));
    end  

    grad =  (sai * x_noisy.') / M   -  inv(B.'); 
    B = B - lr*grad;
    B = B ./ vecnorm(B);

    f(itr) = norm(grad, 'fro')^2;
end

%% Q1: A*B:
permutation_matrix = B*A;
disp('<strong>AB:</strong>')
disp(permutation_matrix)

%% Q2: Source Estimation:

S_hat = B*x_noisy;
S_hat([2 3], :) = S_hat([3 2], :);
alpha =  max(abs(S), [], 2) ./ max(abs(S_hat), [], 2);

for i = 1:N
    if i == 2
        S_hat(i, :) = alpha(i) * S_hat(i, :);
    else
        S_hat(i, :) = -alpha(i) * S_hat(i, :);
    end
end


figure;
subplot(311)
plot(S_hat(1, :), 'r', 'LineWidth', 1.2)
hold on 
plot(S(1, :), 'k', 'LineWidth', 1.2)
legend('S_1hat', 'S_1')
hold off

subplot(312)
plot(S_hat(2, :), 'r', 'LineWidth', 1.2)
hold on 
plot(S(2, :), 'k', 'LineWidth', 1.2)
legend('S_2hat', 'S_2')
hold off

subplot(313)
plot(S_hat(3, :), 'r', 'LineWidth', 1.2)
hold on
plot(S(3, :), 'k', 'LineWidth', 1.2)
legend('S_3hat', 'S_3')
hold off

%% Q2: Error:
err = norm(S_hat - S, 'fro')^2 / norm(S, 'fro')^2;
fprintf('%s\n', repmat('=', 1, 40));
disp(['<strong>Error</strong>: ', num2str(err, '%.6f')])
fprintf('%s\n', repmat('=', 1, 40));

%% Q3: Convergence plot:
figure;
plot(0:max_itr-1, f, 'b', 'LineWidth', 1.2)
title('Convergence over iteration', 'Interpreter', 'latex')
