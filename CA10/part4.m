clc; clear; close all
warning('off', 'all');
load("hw10.mat")
rng(1)

%% Plot:
X = A*S + Noise;

[W, Z] = whitening(X);

max_itr = 1e2;

[N, M] = size(X);
B = randn(N, N);
B = B ./ vecnorm(B);
b1 = B(1, :).';
b2 = B(2, :).';
b3 = B(3, :).';
f = zeros(1, max_itr);

for itr = 1 : max_itr

    y = B*Z;

    [grad1, g1] = compute_gradient(y(1, :), Z);
    b1 = grad1;
    b1 = b1 / vecnorm(b1);
    
    [grad2, g2] = compute_gradient(y(2, :), Z);
    b2 = grad2;
    b2 = (eye(N) - b1*b1')*b2;
    b2 = b2 / vecnorm(b2);

    [grad3, g3] =  compute_gradient(y(3, :), Z);
    b3 = grad3;
    b3 = (eye(N) - [b1,b2]*[b1,b2].') * b3;
    b3 = b3 / vecnorm(b3);

    B = [b1.'; b2.'; b3.'];
    K = [g1, g2, g3];
    f(itr) = norm(g1);
end

%% Q1: BWA
permutation_matrix = B*W*A;
disp('<strong>Permutation Matrix:</strong>')
disp(permutation_matrix)


%% Q2: Source Estimation:
S_hat = B*Z;
 
S_hat(1, :) = - S_hat(1, :);
S_hat(2, :) = - S_hat(2, :);

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

%% Function:

function [gradient_value, expected_Gy] = compute_gradient(y, Z)
    G_y = -exp(-y.^2 / 2);                  
    g_y = 2 .* y .* exp(-y.^2 / 2);
    
    expected_Gy = mean(G_y);                      
    gradient_value = mean(Z .* g_y, 2); 
end

