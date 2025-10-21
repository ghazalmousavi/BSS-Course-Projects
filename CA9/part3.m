clc; clear; close all
warning('off', 'all');

load("hw9.mat")
rng(10)

%% Define constants
X = A*S + Noise;

max_itr = 1e3;
lr = 0.01;
[N, M] = size(X);

B = randn(N, N);
f = zeros(1, max_itr);

%% Q1: Optimization for learning matrix B:

for itr = 1 : max_itr
    y  = B*X;
    sai = zeros(size(X));

    for i = 1 : N
        sai(i, :) = mse_score_function(y(i, :));
    end  

    grad =  (sai * X.') / M   -  inv(B.'); 
    B = (eye(N) - lr * grad * B.')*B;
    B = B ./ vecnorm(B);

    f(itr) = norm(grad, 'fro')^2;
end

%% Q1: BA:
permutation_matrix = B*A;
disp('<strong>Permutation Matrix:</strong>')
disp(permutation_matrix)

%% Q2: Source Estimation:
S_hat = B*X;
S_hat([1 2], :) = S_hat([2 1], :);
S_hat([1 3], :) = S_hat([3 1], :);

alpha =  max(abs(S), [], 2) ./ max(abs(S_hat), [], 2);

for i = 1:N
    S_hat(i, :) = -alpha(i) * S_hat(i, :);
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
plot(1:max_itr, f, 'b', 'LineWidth', 1.2)
title('Convergence over iteration', 'Interpreter', 'latex')
xlim([0 700])