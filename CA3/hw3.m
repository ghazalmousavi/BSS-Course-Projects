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

%% Calculating Observations
X = A*S;

figure(1)
scatter3(X(1, :), X(2, :), X(3, :), 34, [0.5, 0, 0.5])
xlabel('$x_1$', 'FontSize', 14, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
ylabel('$x_2$', 'FontSize', 14, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
zlabel('$x_3$', 'FontSize', 14, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');

title('Scattering Plot', 'FontSize', 16, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
grid on; 
%% Calculating Correlation Matrix:

Rx = X*(X');
[U, D] = eig(Rx);
[sortedD, idx] = sort(diag(D), 'descend'); 

D = diag(sortedD);
disp('Matrix D:');
disp(D);

U = U(:, idx);
disp('Matrix U:');
disp(U);
disp('=================================================================================================');

%% Part B:Steepest Descend:
u1_sd = ones(3, 1);
u2_sd = ones(3, 1);
u3_sd = ones(3, 1);


max_iterations = 50;
mu = 0.01;

[u1_sd, d1_sd] = steepest_descend(u1_sd, Rx, mu, max_iterations);
[u2_sd, d2_sd] = steepest_descend(u2_sd, Rx, mu, max_iterations, u1_sd);
[u3_sd, d3_sd] = steepest_descend(u3_sd, Rx, mu, max_iterations, u1_sd, u2_sd);

Us = [u1_sd, u2_sd, u3_sd];
disp('Matrix U:');
disp(Us)

Ds = diag([d1_sd, d2_sd, d3_sd]);
disp('Matrix D:');
disp(Ds)
disp('=================================================================================================');


%% Part B: Newton

u1_n = ones(3, 1);
u2_n = ones(3, 1);
u3_n = ones(3, 1);
max_iterations = 1500;

[u1_n, d1_n]= newton(u1_n, Rx, max_iterations);

[u2_n, d2_n]= newton(u2_n, Rx, max_iterations, u1_n);

[u3_n, d3_n]= newton(u3_n, Rx, max_iterations, u1_n, u2_n);

Un = [u3_n, u2_n, u1_n];
disp('Matrix U:');
disp(Un)

Dn = diag([d3_n, d2_n, d1_n]);
disp('Matrix D:');
disp(Dn)

disp('=================================================================================================');

%% Part C:

U_reduced = U(:, [1, 2]);
D_reduced = D(:, [1, 2]);
D_reduced(3, :) = [];


C = U_reduced \A;
disp('Matrix C:');
disp(C)
disp('=================================================================================================');

u_3A = U(:, 3)'*A;
%% Part D:
B = D_reduced^(-0.5)*U_reduced';
disp('Matrix B:');
disp(B)

Z = B * X;
Rz = Z * Z';
disp('Matrix Rz:');
disp(Rz)
disp('=================================================================================================');

%%  Part E:
Z_new1 = U_reduced'*X;
fprintf('Energy retention ratio after 2D reduction: %.2f\n', (sum(D(1:2)) / sum(D(1:3)))*100);

figure(2)
scatter(Z_new1(1, :), Z_new1(2, :),34, [0.5, 0, 0.5])
xlabel('$z_1$', 'FontSize', 14, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
ylabel('$z_2$', 'FontSize', 14, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
title('\textbf{Scattering Plot After Dimension Reduction}', 'Interpreter', 'latex');

U_reduced_1d = U(:, 1);
Z_new2 = U_reduced_1d' * X;
fprintf('Energy retention ratio after 1D reduction: %.2f\n', (D(1,1) / sum(diag(D))) * 100);

figure(3)
plot(Z_new2, 34, 'o', 'Color', [0.5, 0, 0.5]); 

xlabel('$z_1$', 'FontSize', 14, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
ylabel('$z_2$', 'FontSize', 14, 'Interpreter', 'latex', 'FontWeight', 'bold', 'Color', 'black');
title('\textbf{Plot After Dimension Reduction}', 'Interpreter', 'latex');

%% Functions:

function [u, d] = steepest_descend(u, Rx, mu, max_iterations, varargin)

    for i = 1:max_iterations

        u(:, i+1) = u(:, i) + 2*mu*Rx*u(:, i);

        if length(varargin) == 1
            u1 = cell2mat(varargin(1));
            u(:, i+1) = u(:, i+1) - (u(:, i+1)'*u1)*u1;
        elseif length(varargin) == 2 

             u1 = cell2mat(varargin(1));
             u2 = cell2mat(varargin(2));

             u(:, i+1) = u(:, i+1) - (u(:, i+1)'*u1)*u1;
             u(:, i+1) = u(:, i+1) - (u(:, i+1)'*u2)*u2;
        end 

        u(:, i+1) = u(:, i+1)./ vecnorm(u(:, i+1));
   
    
        if abs(u(:, i+1)'*Rx*u(:, i+1)- u(:, i)'*Rx*u(:, i)) < 0.0001
   
            fprintf("Steepest Descent converged with mu = %.2f after %d iterations.\n", mu, i);
            u = u(:, end);
            
            break;
        end
    end
    d = u'*Rx*u;
    
end


function [u, d]= newton(u, Rx, max_iterations, varargin)

   for i = 1:max_iterations
        H  = 2*Rx;
        u(:, i+1) = u(:, i) + pinv(H)*(2*Rx*u(:, i));

        if length(varargin) == 1
            u1 = cell2mat(varargin(1));
            u(:, i+1) = u(:, i+1) - (u(:, i+1)'*u1)*u1;
        elseif length(varargin) == 2 
    
             u1 = cell2mat(varargin(1));
             u2 = cell2mat(varargin(2));
    
             u(:, i+1) = u(:, i+1) - (u(:, i+1)'*u1)*u1;
             u(:, i+1) = u(:, i+1) - (u(:, i+1)'*u2)*u2;
        end 
    
        u(:, i+1) = u(:, i+1)./ vecnorm(u(:, i+1));
    
        if abs(u(:, i+1)'*Rx*u(:, i+1)- u(:, i)'*Rx*u(:, i)) < 0.0001 
    
            fprintf("Newton converged after %d iterations.\n", i);
    
            u = u(:, end);
    
            break;
        end
    end   
    d = u'*Rx*u;
end

