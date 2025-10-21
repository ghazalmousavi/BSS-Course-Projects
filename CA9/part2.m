clc; clear; close all
warning('off', 'all');

load("hw9.mat")
rng(40)

%% Define constants
X = A*S + Noise;
[W, Z] = whitening(X);

max_itr = 1e3;
lr = 0.1;

[N, M] = size(X);
B = randn(N, N);
b1 = B(1, :).';
b2 = B(2, :).';
b3 = B(3, :).';
f = zeros(1, max_itr);

%% Q1: Optimization for learning matrix B:
for itr = 1 : max_itr

    y = B*Z;
    
    sai_1 = mse_score_function(y(1, :));
    b1 = b1 - lr * ((sai_1*Z.')/ M).';
    b1 = b1 ./ vecnorm(b1);

    sai_2 = mse_score_function(y(2, :));
    b2 = b2 - lr * ((sai_2*Z.') / M).';
    b2 = (eye(N) - b1*b1.') * b2;
    b2 = b2 ./ vecnorm(b2);


    sai_3 = mse_score_function(y(3, :));
    b3 = b3 - lr * ((sai_3*Z.') / M).';
    b3 = (eye(N) - [b1,b2]*[b1,b2].') * b3;

    b3 = b3 ./ vecnorm(b3);

    B = [b1.'; b2.'; b3.'];
    grad = [sai_1*Z.'; sai_2*Z.'; sai_3*Z.']/M;
    f(itr) = norm(grad)^2;
end

%% Q1: BWA:
permutation_matrix = B*W*A;
disp('<strong>Permutation Matrix:</strong>')
disp(permutation_matrix)


%% Q2: Source Estimation:
S_hat = B*Z;
S_hat([2 3], :) = S_hat([3 2], :);
alpha =  max(abs(S), [], 2) ./ max(abs(S_hat), [], 2);

for i = 1:N
    if i == 3
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
plot(1:max_itr, f, 'b', 'LineWidth', 1.2)
title('Convergence over iteration', 'Interpreter', 'latex')

%% Functions:
function [W, Z] = whitening(X)
    R = X * X';
    [V, D] = eig(R);
    [sortedD, idx] = sort(diag(D), 'descend');
    D = diag(sortedD);
    V = V(:, idx);
    W = D^(-0.5) * V';
    Z = W * X;
end