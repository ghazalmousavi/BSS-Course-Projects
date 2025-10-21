clc;
clear;
close all;

%% Defining constants and Mixture Function:
T = 10000;
A = [0.6, 0.7071, 0.8;
     0.8, 0.7071, -0.6];
a = -2;
b = 2;

%% Defining Sources 
z1 = unifrnd(a, b, T, 1);
z2 = unifrnd(a, b, T, 1);
z3 = unifrnd(a, b, T, 1);
Z = [z1'; z2'; z3'];

source_idx = randi([1, 3], 1, T);  

S = sparse(source_idx, 1:T, Z(sub2ind(size(Z), source_idx, 1:T)), 3, T);

%% Calculating Observations
X = A*S;
x1 = X(1, :);
x2 = X(2, :);

%% 1, 2.
scattering_plot(x1, x2)
Beta = mixture_function(x1, x2);

%% 3.
R0 = Beta'*X;
[values, indices] = max(R0, [], 1);
Shat = zeros(3, T);
Shat(sub2ind(size(Shat), indices, 1:T)) = values;

%% Functions:

function  [Beta, point1, point2, point3, point4, point5, point6] = mixture_function(x1, x2)

    [x2_min, index1] = min(x2);
    [x2_max, index2] = max(x2);
    [x1_min, index3] = min(x1);
    [x1_max, index4] = max(x1);
    
    
    
    point1 = [x1(index1), x2_min];
    point2 = [x1(index2), x2_max];
    
    point5 = [x1_min, x2(index3)];
    point6 = [x1_max, x2(index4)];
    
    x3 = x1 == x2;
    
    point3 = [min(x2(x3)), min(x2(x3))];
    point4 = [max(x2(x3)), max(x2(x3))];
    
    beta1 = (point2 - point1)/ norm(point2 - point1);
    beta2 = (point4 - point3)/ norm(point4 - point3);
    beta3 = (point6 - point5)/ norm(point6 - point5);
    Beta = [beta1', beta2', beta3'];

end


function  scattering_plot(x1, x2)
    scatter(x1, x2, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0.8, 0.8, 0.8]);
    xlabel('x_1', 'FontSize', 14, 'FontName', 'Times', 'FontWeight', 'bold', 'Color', 'black');
    ylabel('x_2', 'FontSize', 14, 'FontName', 'Times', 'FontWeight', 'bold', 'Color', 'black');
    title('Scattering Plot', 'FontSize', 16, 'FontName', 'Times', 'FontWeight', 'bold', 'Color', 'black');
    grid on; 
    
    [Beta, point1, point2, point3, point4, point5, point6] = mixture_function(x1, x2);
    
    
    figure;
    scatter(x1, x2, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0.8, 0.8, 0.8]);
    hold on;
    plot([point1(1), point2(1)], [point1(2), point2(2)], 'r-', 'LineWidth', 2);
    plot([point3(1), point4(1)], [point3(2), point4(2)], 'r-', 'LineWidth', 2);
    plot([point5(1), point6(1)], [point5(2), point6(2)], 'r-', 'LineWidth', 2);
    hold off
    xlabel('x_1', 'FontSize', 14, 'FontName', 'Times', 'FontWeight', 'bold', 'Color', 'black');
    ylabel('x_2', 'FontSize', 14, 'FontName', 'Times', 'FontWeight', 'bold', 'Color', 'black');
    title('Scattering Plot', 'FontSize', 16, 'FontName', 'Times', 'FontWeight', 'bold', 'Color', 'black');
    grid on; 
    
    disp(table(Beta));
       
end
