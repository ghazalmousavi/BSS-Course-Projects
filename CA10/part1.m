clc; clear; close all
warning('off', 'all');
load("hw10.mat")
rng(1)

%% Plot:
x = A*S;
X = A*S + Noise;

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
plot(X(1, :), 'LineWidth', 1.2)
hold on 
plot(X(2, :), 'LineWidth', 1.2)
plot(X(3, :), 'LineWidth', 1.2)
legend('X_1', 'X_2', 'X_3')
title('Observation with noise', 'Interpreter', 'latex')
hold off

%% 
[W, Z] = whitening(X);

max_itr = 1e2;
lr = 0.01;

[N, M] = size(X);
B = randn(N, N);
B = B ./ vecnorm(B);
b1 = B(1, :).';
b2 = B(2, :).';
b3 = B(3, :).';
f = zeros(1, max_itr);

b = zeros(1, max_itr);
for itr = 1 : max_itr

    y = B*Z;

    [grad1, k1] = kurtosis_gradient(y(1, :), b1, Z);
    b1 = b1 + lr*grad1;
    b1 = b1 / vecnorm(b1);

    [grad2, k2] = kurtosis_gradient(y(2, :), b2, Z);
    b2 = b2 + lr*grad2;
    b2 = (eye(N) - b1*b1')*b2;
    b2 = b2 / vecnorm(b2);
    b(itr) = b2(1);


    [grad3, k3] = kurtosis_gradient(y(3, :), b3, Z);
    b3 = b3 + lr*grad3;
    b3 = (eye(N) - [b1,b2]*[b1,b2].') * b3;
    b3 = b3 / vecnorm(b3);

    B = [b1.'; b2.'; b3.'];
    
    K = [k1, k2, k3];
    f(itr) = mean(K);
end


%% Q1: BWA:
permutation_matrix = B*W*A;
disp('<strong>Permutation Matrix:</strong>')
disp(permutation_matrix)

%% Q2: Source Estimation:

S_hat = B*Z;
alpha = ones(1, N);

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
