clc;
clear;
close all
rng(42)

%% Defining constants and Mixture Function:

T = 1000;
A = [1, -2;
     2, -1;
     3, -2];

%% Defining Sources 

s1 = unifrnd(-3, 3, T, 1);
s1 = s1 - mean(s1);

s2 = unifrnd(-2, 2, T, 1);
s2= s2 - mean(s2);

S = [s1'; s2'];

%% Part A:

X = A*S;
Rx = X*(X');

figure(1)
scatter3(X(1, :), X(2, :), X(3, :), 34, [0.5, 0.1, 0.75])
xlabel('$x_1$', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
ylabel('$x_2$', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
zlabel('$x_3$', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
title('Scatter Plot of Observations', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', [0, 0, 0]);
grid on; 


[U, D] = eig(Rx);
[sortedD, idx] = sort(diag(D), 'descend'); 
D = diag(sortedD);
U = U(:, idx);

disp('Matrix D:');
disp(D);
disp('Matrix U:');
disp(U);
disp('=====================================================');


%% Part C:

U_sig = U(:, [1, 2]);
D_reduced = D(:, [1, 2]);
D_reduced(3, :) = [];

B = D_reduced^(-0.5)*U_sig';
disp('Matrix B:');
disp(B)

Z = B * X;
Rz = Z * Z';
disp('Matrix Rz:');
disp(Rz)
disp('=====================================================');


figure(2)
scatter(Z(1, :), Z(2, :),34, [0.5, 0.1, 0.75])
xlabel('$z_1$', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
ylabel('$z_2$', 'FontSize', 10, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
title('\textbf{Scattering Plot After Dimension Reduction}', 'Interpreter', 'latex');

%% Part D:

[Q, G, V_T] = svd(X);

fprintf('Number of non-zero elements in G: %d\n', sum(diag(G) > 1e-10));
disp('=====================================================');

v1 = V_T(:, 1);
v2 = V_T(:, 2);
V_sig = [v1'; v2'];

%% Part E:

F = S * pinv(Z);
disp('Matrix F:');
disp(F)
disp('=====================================================');
