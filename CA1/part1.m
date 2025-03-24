clc;
clear;
close all
rng(42)
%% Defining constants and Mixture Function:
T = 1000;
A = [0.6, 0.8;
     0.8, -0.6];

%% Defining Sources 
s1 = unifrnd(-3, 3, T, 1);
s2 = unifrnd(-2, 2, T, 1);
S = [s1'; s2'];

%% Calculating Observations
X = A*S;
x1 = X(1, :);
x2 = X(2, :);

%% 1.
scattering_plot(x1, x2)
Beta = mixture_function(x1, x2);
error = mse(Beta, A);
disp(['MSE:', num2str(error)]);

%% 2. Adding Noise to Observation:
noise = 0.1*randn(1, T);
noisy_x1 = x1 + noise;
noisy_x2 = x2 + noise;

figure;
scattering_plot(noisy_x1, noisy_x2)
noisy_Beta = mixture_function(noisy_x1, noisy_x2);

error = mse(noisy_Beta, A);
disp(['MSE:', num2str(error)]);

%% Histogram:

figure;
h1 = histogram(x1, 'NumBins', 20);

h1.FaceColor = [0.8, 0.8, 0.8]; 
h1.EdgeColor = [0, 0, 0];    
h1.LineWidth = 2;

title('Histogram of x_1', 'FontSize', 12, 'FontName', 'Times', 'FontWeight', 'bold', 'Color', 'black');

figure;
h2 = histogram(x2, 'NumBins', 20);

h2.FaceColor = [0.8, 0.8, 0.8]; 
h2.EdgeColor = [0, 0, 0];    
h2.LineWidth = 2;

title('Histogram of x_2', 'FontSize', 12, 'FontName', 'Times', 'FontWeight', 'bold', 'Color', 'black');



%% Functions:

function  [Beta, point1, point2, point3] = mixture_function(x1, x2)

    [x1_min, index1] = min(x1);
    [x2_min, index2] = min(x2);
    [x1_max, index3] = max(x1);
    
    point1 = [x1_min, x2(index1)];
    point2 = [x1(index2), x2_min];
    point3 = [x1_max, x2(index3)];
    
    beta1 = (point3 - point2)/ norm(point3 - point2);
    beta2 = (point2 - point1)/ norm(point2 - point1);
    Beta = [beta1', beta2'];

end


function  scattering_plot(x1, x2)
    scatter(x1, x2, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0.8, 0.8, 0.8]);
    xlabel('x_1', 'FontSize', 14, 'FontName', 'Times', 'FontWeight', 'bold', 'Color', 'black');
    ylabel('x_2', 'FontSize', 14, 'FontName', 'Times', 'FontWeight', 'bold', 'Color', 'black');
    title('Scattering Plot', 'FontSize', 16, 'FontName', 'Times', 'FontWeight', 'bold', 'Color', 'black');
    grid on; 
    
    [Beta, point1, point2, point3] = mixture_function(x1, x2);
    
    figure;
    scatter(x1, x2, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0.8, 0.8, 0.8]);
    hold on;
    plot([point1(1), point2(1)], [point1(2), point2(2)], 'r-', 'LineWidth', 2);
    plot([point2(1), point3(1)], [point2(2), point3(2)], 'r-', 'LineWidth', 2);
    hold off
    xlabel('x_1', 'FontSize', 14, 'FontName', 'Times', 'FontWeight', 'bold', 'Color', 'black');
    ylabel('x_2', 'FontSize', 14, 'FontName', 'Times', 'FontWeight', 'bold', 'Color', 'black');
    title('Scattering Plot', 'FontSize', 16, 'FontName', 'Times', 'FontWeight', 'bold', 'Color', 'black');
    grid on; 
    
    disp(table(Beta));
       
end


