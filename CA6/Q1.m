clc; clear; close all
% rng('shuffle')
%% Defining Sources, Mixture Matrix and Observations:

c = [0.2, 0.4, 0.6, -0.1, -0.3];
d = [0.1, 0.3, -0.2, 0.5, -0.3];

A = [0.8, -0.6;
    0.6, 0.8];

K = 5;
fs = 20;
t = 0:1/fs:K-1/fs;

s1 = zeros(fs*K, 1);
s2 = zeros(fs*K, 1);

for i=1:K
    s1((i-1)*fs + 1 : i*fs) = c(i) * sin(2 * pi * t((i-1)*fs + 1 : i*fs)).'; 
    s2((i-1)*fs + 1 : i*fs) = d(i) * sin(4 * pi * t((i-1)*fs + 1 : i*fs)).'; 
end

S = [s1'; s2'];
X = A*S;

%% Part 1: Plotting Sources and Observations

figure(1)

subplot(211)
plot(t, s1, 'Color', [0.5, 0.1, 0.75], 'LineWidth', 1.5)
xlabel('Time', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
ylabel('$s_1$', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
title('$s_1(t)$', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', [0, 0, 0]);
grid on;

subplot(212)
plot(t, s2, 'Color', [0.5, 0.1, 0.75], 'LineWidth', 1.5)
xlabel('Time', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
ylabel('$s_2$', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
title('$s_2(t)$', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', [0, 0, 0]);
grid on
sgtitle('Sources', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');


figure(2)
subplot(211)
plot(t, X(1, :), 'Color', [0.5, 0.1, 0.75], 'LineWidth', 1.5)
xlabel('Time', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
ylabel('$x_1$', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
title('$x_1(t)$', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', [0, 0, 0]);
grid on;

subplot(212)
plot(t, X(2, :), 'Color', [0.5, 0.1, 0.75], 'LineWidth', 1.5)
xlabel('Time', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
ylabel('$x_2$', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
title('$x_2(t)$', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', [0, 0, 0]);
grid on;
sgtitle('Observations', 'Interpreter', 'latex','FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');

%% Part B: Estimating Sources Using Two Windows 

X_1 = X(:, 1:fs);
R_X_1 = X_1*X_1.';

X_2 = X(:, fs + 1:2*fs);
R_X_2 = X_2*X_2.';

[U, D] = eig(pinv(R_X_2)*R_X_1);
B = U.';

S_hat = B*X;

E = norm(S - S_hat, 'fro')^2/ norm(S, 'fro')^2;
fprintf('Error with using two windows: %.2e\n', E);
fprintf('%s\n', repmat('-', 1, 40));

%% Part C: Estimating Sources Using All Windows

max_iterations = 500;
b2 = randn(2,1);
f = zeros(1, max_iterations);

Z = whitening(X);
Rz_cells = compute_Rz(Z, K);

for i = 1: max_iterations

    R1 = zeros(2, 2);
    for k = 1:K
        temp = Rz_cells{k} * b2;
        R1 = R1 + temp*(temp');  
    end
   
    
    [b, d] = eig(R1);
    [~, index] = min(diag(d));     
    b1 = b(:, index);
    b1 = b1 ./ norm(b1);


    R2 = zeros(2, 2);
    for k = 1:K
        temp = Rz_cells{k} * b1;
        R2 = R2 + temp*(temp'); 
    end

    [b, d] = eig(R2);
    [~, index] = min(diag(d));
    b2 = b(:, index);
    b2 = b2 - (b2' * b1) * b1;
    b2 = b2 ./ norm(b2);


    f(i) = obj_func(b1, b2, Rz_cells);
    if i > 1 && abs(f(i) - f(i-1)) < 1e-15
        break;
    end

end

B2 = [b1'; b2'];
S_hat2 = B2*Z;

[error, s_hat] = select_best_estimation(S, S_hat);

fprintf('Error with using all windows: %.2e\n', error);
fprintf('%s\n', repmat('-', 1, 40));

%% Functions:

function Z = whitening(X)
    R = X * X';
    [V, D] = eig(R);
    [sortedD, idx] = sort(diag(D), 'descend');
    D = diag(sortedD);
    V = V(:, idx);
    Z = (D^(-0.5) * V') * X;
end



function f = obj_func(b1, b2, Rz)
    K = length(Rz);
    f = 0;
    for k = 1:K
        Rk = Rz{k};
        f = f + (b1' * Rk * b2)^2 + (b2' * Rk * b1)^2;
    end
      
end


function Rz_cells = compute_Rz(Z, K)
    Rz_cells = cell(1, 5);
    fs = length(Z)/K;

    for k = 1:K
        idx = (k-1)*fs + 1 : k*fs;
        Rz_k = Z(:, idx) * Z(:, idx)';
        Rz_cells{k} = Rz_k;
    end

end


function [error, s_hat] = select_best_estimation(S, S_hat)

    S = (S - mean(S, 2)) ./ std(S, 0, 2);
    S_hat = (S_hat - mean(S_hat, 2)) ./ std(S_hat, 0, 2);

    Error1 = norm(S - S_hat, 'fro')^2 / norm(S, 'fro')^2;

    S_hat_flipped = S_hat([2 1], :);
    S_hat_flipped = (S_hat_flipped - mean(S_hat_flipped, 2)) ./ std(S_hat_flipped, 0, 2);


    Error2 = norm(S - S_hat_flipped, 'fro')^2 / norm(S, 'fro')^2;

    if Error1 <= Error2
        error = Error1;
        s_hat = S_hat;
    else
        error = Error2;
        s_hat = S_hat_flipped;
    end
end
