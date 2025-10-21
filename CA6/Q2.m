clc; clear; close all;

%% Loading Data:
load("hw6-X1.mat")
load("hw6-X2.mat")

%% Defining constants:

fs = 100;
shift = 20;
max_iterations = 100;
taus = [5, 10];
K = 2;

%% Part 1: Estimating and Plotting Sources Using X1 Observations

Z1 = whitening(X1);
Rzk= cross_covariances(Z1, taus);

[B1_t, ~] = eig(Rzk{1}, Rzk{2});

S_hat1 = B1_t'*Z1;
t = 0:1/fs:5-1/fs;

figure(1)
subplot(211)
plot(t, S_hat1(1, :), 'b', 'LineWidth', 0.7)
xlabel('Time', 'Interpreter', 'latex')
ylabel('$\hat{s}_1$', 'Interpreter', 'latex')
grid on

subplot(212)
plot(t, S_hat1(2, :), 'r', 'LineWidth', 0.7)
xlabel('Time', 'Interpreter', 'latex')
ylabel('$\hat{s}_2$', 'Interpreter', 'latex')
grid on

sgtitle('$\hat{S}$', 'Interpreter', 'latex', 'FontSize', 12)


%% Part 2: Plotting Fourier Transform of Sources and Observations(X1):

N = length(t);                 
f = linspace(0, fs, N);

fft_s11 = abs(fft(S_hat1(1, :)));
fft_s21 = abs(fft(S_hat1(2, :)));

figure(2);
subplot(211)
plot(f, fft_s11, 'b', 'LineWidth', 1);
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel('Magnitude', 'Interpreter', 'latex');
grid on

subplot(212)
plot(f, fft_s21, 'r', 'LineWidth', 1);
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel('Magnitude', 'Interpreter', 'latex');
sgtitle('FFT of $\hat{S}$', 'Interpreter', 'latex', 'FontSize', 12);
grid on;

%% Part 3:Estimating and Plotting Sources Using X2 Observations

Z2 = whitening(X2);
Rzk2= cross_covariances(Z2, taus);

[B2_t, D] = eig(Rzk2{1}, Rzk2{2});

S_hat2 = B2_t'*Z2;

figure(3)
subplot(211)
plot(t, S_hat2(1, :), 'b', 'LineWidth', 0.7)
xlabel('Time', 'Interpreter', 'latex')
ylabel('$\hat{s}_1$', 'Interpreter', 'latex')
grid on

subplot(212)
plot(t, S_hat2(2, :), 'r', 'LineWidth', 0.7)
xlabel('Time', 'Interpreter', 'latex')
ylabel('$\hat{s}_2$', 'Interpreter', 'latex')
sgtitle('$\hat{S}$', 'Interpreter', 'latex', 'FontSize', 12)
grid on


%% Part 4: Plotting Fourier Transform of Sources and Observations(X2):

fft_s12 = abs(fft(S_hat2(1, :)));
fft_s22 = abs(fft(S_hat2(2, :)));

figure(4);
subplot(211)
plot(f, fft_s12, 'b', 'LineWidth', 1);
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel('Magnitude', 'Interpreter', 'latex');
grid on

subplot(212)
plot(f, fft_s22, 'r', 'LineWidth', 1);
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel('Magnitude', 'Interpreter', 'latex');
sgtitle('FFT of $\hat{S}$', 'Interpreter', 'latex', 'FontSize', 12);
grid on;


%% Functions:

function Z = whitening(X)
    R = X * X';
    [V, D] = eig(R);
    [sortedD, idx] = sort(diag(D), 'descend');
    D = diag(sortedD);
    V = V(:, idx);
    Z = (D^(-0.5) * V') * X;
end


function Rz = cross_covariances(Z1, taus)

    Rz = cell(1, length(taus));

    for i = 1:length(taus)
        tau = taus(i);
        Z1_t   = Z1(:, (tau+1):end);
        Z1_tau = Z1(:, 1:(end - tau));
        Rz{i} = Z1_t * Z1_tau';
    end
end

