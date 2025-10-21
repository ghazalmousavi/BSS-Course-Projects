clc;
clear;
close all
rng(42);

%% Defining constants:

fc = 150e6;
f1 = 20e3;
f2 = 10e3;
fs = 1e6;
M = 10;
N = 2;
c = 3e8;
T = 10000;
d = 1;
k = 2*pi*fc/c;
theta1 = deg2rad(10);
theta2 = deg2rad(20);
t = 0: 1/fs: 0.01 - 1/fs;
alpha = (0:M-1)';
noise = randn(M, T);

%% Defining Sources & Mixture Function 

s1 = exp(1i*2*pi*f1.*t).';
s2 = exp(1i*2*pi*f2.*t).';


a1 = exp(-1i*k*alpha*d*sin(theta1));
a2 = exp(-1i*k *alpha*d*sin(theta2));

A = [a1, a2];
S = [s1'; s2'];

%% SVD:

Y = A*S + noise;

[U, Q, V_T] = svd(Y);
V = V_T.';

U_sig = U(:, [1, 2]);
U_null = U(:, setdiff(1:size(U, 2), [1, 2]));


V_sig = V(:, [1, 2]);
V_null = V(:, setdiff(1:size(V, 1), [1, 2]));
%% Part C: 

theta = 0:90;
a =  exp(-1i*k*alpha*d*sin(deg2rad(theta)));
F = vecnorm(a'* U_sig, 2, 2);
[~, indices] = maxk(F, 2);

figure(1)
plot(theta, F, 'b', 'LineWidth', 2);
hold on;
scatter(theta(indices), F(indices), 50, 'r', 'filled');
xlabel('$\theta$ (degrees)', 'Interpreter','latex');
ylabel('F', 'Interpreter', 'latex');
title('Beamforming Method', 'Interpreter', 'latex', 'FontSize', 10, 'FontWeight', 'bold');
grid on;

%% Part D: 

G = 1 ./ vecnorm(a'* U_null, 2, 2);
[~, indices] = maxk(G, 2);

figure(2)
plot(theta, G, 'b', 'LineWidth', 2);
hold on;
scatter(theta(indices), G(indices), 50, 'r', 'filled');
xlabel('$\theta$ (degrees)', 'Interpreter','latex');
ylabel('G', 'Interpreter', 'latex');
title('MUSIC Method', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

%% Part E: 

f = 0:1000:5e4;
s =  exp(1i*2*pi*f'*t);

H = vecnorm(s* V_sig, 2, 2);
[~, indices] = maxk(H, 2);

figure(3)
plot(f, H, 'b', 'LineWidth', 2);
hold on;
scatter(f(indices), H(indices), 50, 'r', 'filled');
xlabel('f (Hz)', 'Interpreter','latex');
ylabel('H', 'Interpreter', 'latex');
title('Beamforming Method', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

%% Part F:

I = 1./vecnorm(s* V_null, 2, 2);
[~, indices] = maxk(I, 2);

figure(4)
plot(f, I, 'b', 'LineWidth', 2);
hold on;
scatter(f(indices), I(indices), 50, 'r', 'filled');
xlabel('f (Hz)', 'Interpreter','latex');
ylabel('I', 'Interpreter', 'latex');
title('MUSIC Method', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
grid on;