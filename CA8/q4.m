clc; clear; close all

load('hw8-part4.mat')
rng(2)
%% Part 2:

L = 100;
K = 5;

max_itr = 200;

[s1_hat, tau1, alpha1, x_hat1] = Blind_Deconvolution(x1, L, K);
[s2_hat, tau2, alpha2, x_hat2] = Blind_Deconvolution(x2, L, K);

[~, T] = size(x1);

figure;

subplot(3,1,1)
plot(1:T, x1, 'b', 'LineWidth', 1.5)
hold on 
plot(1:T, x_hat1, 'r', 'LineWidth', 1)
legend('Real $x_1$', 'Estimated $\hat{x}_1$', 'Interpreter', 'latex')
title('Original vs Estimated Signal', 'Interpreter', 'latex')
xlabel('Time Index', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
hold off
grid on

subplot(3,1,2)
plot(s1_hat, 'b', 'LineWidth', 1.5)
xlabel('Sample Index', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
title('Estimated Spike Signal $\hat{s}_1$', 'Interpreter', 'latex')
grid on

subplot(3,1,3)
stem(tau1, alpha1, 'b', 'LineWidth', 1.5) 
xlabel('Sample Index', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
title('Impulse Train with Amplitudes $\alpha_1$', 'Interpreter', 'latex')
xlim([0 2000])
grid on



figure;

subplot(3,1,1)
plot(1:T, x2, 'b', 'LineWidth', 1.5)
hold on 
plot(1:T, x_hat2, 'r', 'LineWidth', 1)
legend('Real $x_2$', 'Estimated $\hat{x}_2$', 'Interpreter', 'latex')
title('Original vs Estimated Signal', 'Interpreter', 'latex')
xlabel('Time Index', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
hold off
grid on

subplot(3,1,2)
plot(s2_hat, 'b', 'LineWidth', 1.5)
xlabel('Sample Index', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
title('Estimated Spike Signal $\hat{s}_2$', 'Interpreter', 'latex')
grid on


subplot(3,1,3)
stem(tau2, alpha2, 'b', 'LineWidth', 1.5) 
xlabel('Sample Index', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
xlim([0 2000])
title('Impulse Train with Amplitudes $\alpha_2$', 'Interpreter', 'latex')
grid on
%% Part 3:

figure;
subplot(211)
plot(1:T, X(1, :), 'b', 'LineWidth', 1.5)
title('Channel 1', 'Interpreter', 'latex')
grid on
subplot(212)
title('Channel 2', 'Interpreter', 'latex')
plot(1:T, X(2, :), 'b', 'LineWidth', 1.5)
grid on

A_hat = randn(2, 2);
A_hat = A_hat ./ vecnorm(A_hat);
s_hat = zeros(2, L);
alpha = zeros(2, K);
tau = zeros(2, K);
for itr = 1:max_itr
    B_hat = pinv(A_hat) * X;
    for row = 1:size(B_hat, 1)
        [s_hat(row, :), tau(row, :), alpha(row, :),B_hat(row, :)] = Blind_Deconvolution(B_hat(row, :), L, K);
    end
    
    A_hat = X * pinv(B_hat);
    A_hat = A_hat ./ vecnorm(A_hat);
end
X_hat = A_hat * B_hat;

% Channel 1:
figure;
subplot(3,1,1)
plot(1:T, X(1, :), 'b', 'LineWidth', 1.5)
hold on 
plot(1:T, X_hat(1, :), 'r', 'LineWidth', 1)
legend('Real $X_1$', 'Estimated $\hat{X}_1$', 'Interpreter', 'latex')
title('Original vs Estimated Signal', 'Interpreter', 'latex')
xlabel('Time Index', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
hold off
grid on

subplot(3,1,2)
plot(s_hat(1, :), 'b', 'LineWidth', 1.5)
xlabel('Sample Index', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
title('Estimated Spike Signal $\hat{s}_1$', 'Interpreter', 'latex')
grid on


subplot(3,1,3)
stem(tau(1, :), alpha(1, :), 'b', 'LineWidth', 1.5) 
xlabel('Sample Index', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
xlim([0 2000])
title('Impulse Train with Amplitudes $\alpha_1$', 'Interpreter', 'latex')
grid on

% Channel 2:
figure;

subplot(3,1,1)
plot(1:T, X(2, :), 'b', 'LineWidth', 1.5)
hold on 
plot(1:T, X_hat(2, :), 'r', 'LineWidth', 1)
legend('Real $X_1$', 'Estimated $\hat{X}_2$', 'Interpreter', 'latex')
title('Original vs Estimated Signal', 'Interpreter', 'latex')
xlabel('Time Index', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
hold off
grid on

subplot(3,1,2)
plot(s_hat(2, :), 'b', 'LineWidth', 1.5)
xlabel('Sample Index', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
title('Estimated Spike Signal $\hat{s}_1$', 'Interpreter', 'latex')
grid on


subplot(3,1,3)
stem(tau(2, :), alpha(2, :), 'b', 'LineWidth', 1.5) 
xlabel('Sample Index', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
xlim([0 2000])
title('Impulse Train with Amplitudes $\alpha_1$', 'Interpreter', 'latex')
grid on

%% Functions:

function [s_hat, tau, alpha, x_hat] = Blind_Deconvolution(x, L, K)

    T = size(x, 2);
    alpha = ones(K, 1);
    tau = round(linspace(1, T - L + 1, K))';
    max_itr = 50;

    y = zeros(L, K);
    
    idx = (1:L)' + (0:T-L-1);  
    z = x(idx); 
    for itr = 1:max_itr
        for k = 1 : K
            y(:, k) = x(tau(k): tau(k) + L -1);
        end
        s_hat = y * alpha;
        s_hat = s_hat ./ norm(s_hat, 'fro');
    
    
        b = s_hat.'*z;
        for k = 1 : K
            [alpha(k), idx] = max(b);
            tau(k) = idx;
            b(max(1, idx-L+1) : min(end, idx+L-1)) = -inf;
        end
    end
   
    x_hat = zeros(1, T);
    for k = 1:K
        x_hat(tau(k):tau(k) + L - 1) = x_hat(tau(k):tau(k) + L - 1) + (s_hat*alpha(k))';
    end


end