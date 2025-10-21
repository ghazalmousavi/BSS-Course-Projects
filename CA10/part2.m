clc; clear; close all
warning('off', 'all');
load("hw10.mat")
rng(26)

%% Plot:
X = A*S + Noise;

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


for itr = 1 : max_itr
    y = B*Z;

    [grad1, k1] = kurtosis_gradient(y(1, :), b1, Z);
    b1 = grad1;
    b1 = b1 / vecnorm(b1);
%end
%for itr = 1 : max_itr
    [grad2, k2] = kurtosis_gradient(y(2, :), b2, Z);
    b2 = grad2;
    b2 = (eye(N) - b1*b1')*b2;
    b2 = b2 / vecnorm(b2);
%end
%for itr = 1 : max_itr
    [grad3, k3] = kurtosis_gradient(y(3, :), b3, Z);
    b3 = grad3;
    b3 = (eye(N) - [b1,b2]*[b1,b2].') * b3;
    b3 = b3 / vecnorm(b3);
%end
    B = [b1.'; b2.'; b3.'];
    G = [k1, k2, k3];
    f(itr) = mean(G);
end
%% Q1: BWA:
permutation_matrix = B*W*A;
disp('<strong>Permutation Matrix:</strong>')
disp(permutation_matrix)


%% Q2: Source Estimation:
S_hat = B*Z;
S_hat([1 2], :) = S_hat([2 1], :);

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
