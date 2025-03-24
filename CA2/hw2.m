clc;
clear;
close all;

%% Q1:

[X1, X2] = meshgrid(-10:.01:10);
F = X1.^2 + X2.^2 -4.*X1 -6.*X2 + 13 + X1.*X2;

figure;
s = mesh(X1,X2, db(F),'FaceAlpha','0.5');
xlabel('$x_1$', 'Interpreter','latex', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('$x_2$', 'Interpreter','latex', 'FontSize', 12, 'FontWeight', 'bold')
zlabel('db(f)', 'Interpreter','latex', 'FontSize', 12, 'FontWeight', 'bold')
title('3D Visualization of the Objective Function \( f(x_1, x_2) \)', ...
    'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');

%% Q2:

figure;
contour(X1,X2, db(F),'ShowText','on')
xlabel('$x_1$', 'Interpreter','latex', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('$x_2$', 'Interpreter','latex', 'FontSize', 12, 'FontWeight', 'bold')
title('Contour Plot of the Objective Function \( f(x_1, x_2) \) using \( db(f) \)', ...
    'Interpreter', 'latex', 'FontSize', 10, 'FontWeight', 'bold');

%% Q3 & Q4:

syms x1 x2
f = x1^2 + x2^2 - 4*x1 - 6*x2 + 13 + x1*x2;

g = gradient(f, [x1, x2]);
h = double(hessian(f, [x1, x2]));

min_point  = solve(g, [x1, x2]);
min_f = double(subs(f, [x1, x2], [min_point.x1, min_point.x2]));

fprintf('Minimum of cost function f: %f at point = (%f, %f)\n', min_f, double(min_point.x1), double(min_point.x2));
disp('=================================================================================================');

%% Q5:

start_point = [6; 6];
max_iterations = 1000;

iteration_points1 = steepest_descend(f, x1, x2, start_point, 0.1, max_iterations);
  

iteration_points2 = steepest_descend(f, x1, x2, start_point, 0.01, max_iterations);

figure;
plot_convergence(f, iteration_points1, x1, x2)
hold on
plot_convergence(f, iteration_points2, x1, x2)
legend('$\mu = 0.1$', '$\mu = 0.01$', 'Interpreter', 'latex');
hold off

%% Q6:

max_iterations = 5;
iteration_points3 = newton(f, x1, x2, start_point, max_iterations);

figure;
plot_convergence(f, iteration_points3, x1, x2)
%% Q7:

max_iterations = 50;
iteration_points4 = alternation_minimization(f, x1, x2, start_point, max_iterations);

figure;
contour(X1, X2, db(F)); 
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
title('Contour Plot of the Objective Function \( f(x_1, x_2) \) using \( db(f) \)', ...
    'Interpreter', 'latex', 'FontSize', 10, 'FontWeight', 'bold');
hold on;
plot(iteration_points4(1, :), iteration_points4(2, :), 'r-o', 'LineWidth', 2, 'MarkerSize', 6);
hold off;

%% Q8:

max_iterations = 30;
iteration_points5 = gradient_projection(f, x1, x2, start_point, 0.1, max_iterations);
figure;
plot_convergence(f, iteration_points5, x1, x2)

%% Functions:

function x = steepest_descend(f, x1, x2, start_point, mu, max_iterations)
    x = zeros(2, max_iterations + 1);  
    x(:, 1) = start_point;
    g = gradient(f, [x1, x2]);

    for i = 1:max_iterations

        x(:, i+1) = x(:, i) - mu * double(subs(g, [x1, x2], x(:, i)'));

        if abs(calc_f_value (x(:, i+1), x1, x2, f) - calc_f_value (x(:, i), x1, x2, f)) < 0.0001

            fprintf("Steepest Descent converged with mu = %.2f after %d iterations.\n", mu, i);
            x = x(:, 1:i+1);
            
            break;
        end
    end

    min_point = x(:, end);   
    min_f = double(subs(f, [x1, x2], [min_point(1), min_point(2)]));
    fprintf('Minimum of cost function f with Steepest Descend (mu = %.2f): %4f at point = (%.4f, %.4f)\n', mu, min_f,  min_point(1), min_point(2));
    disp('=================================================================================================');

end

function f_value = calc_f_value (x, x1, x2, f)
    f_value = double(subs(f, [x1, x2], [x(1), x(2)]));
end


function plot_convergence(f, x, x1, x2)

    f_values = zeros(1, size(x, 2)-1);

    for i=1:size(x, 2)-1

        f_values(i) =  double(subs(f, [x1, x2], x(:, i)'));

    end

    plot(1:size(x, 2)-1, f_values,'LineWidth', 2)

    xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('f(x)', 'Interpreter', 'latex', 'FontSize', 12);
    title('Function Convergence over Iterations', 'Interpreter','latex', 'FontSize', 12, 'FontWeight','bold');
    grid on;
end


function x= newton(f, x1, x2, start_point, max_iterations)
   x = zeros(2,  max_iterations);
   x(:, 1) = start_point;

   g = gradient(f, [x1, x2]);
   h = double(hessian(f, [x1, x2]));

   for i = 1:max_iterations
        x(:, i+1) = x(:, i) -  inv(h)* double(subs(g, [x1, x2], x(:, i)'));

        if abs(calc_f_value (x(:, i+1), x1, x2, f) - calc_f_value (x(:, i), x1, x2, f)) < 0.0001

            fprintf("Newton converged after %d iterations.\n", i);

            x(:, i+1:end) = repmat(x(:, i+1), 1, max_iterations - (i+1) + 1);
%             x = x(:, 1:i+1);
    
            break;
        end
   end
   min_point = x(:, end);
   min_f = double(subs(f, [x1, x2], [min_point(1), min_point(2)]));

   fprintf('Minimum of cost function f with Newton: %4f at point = (%.4f, %.4f)\n', min_f,  min_point(1), min_point(2));
   disp('=================================================================================================');

end

function x = alternation_minimization(f, x1, x2, start_point, max_iterations)
    x = zeros(2,  max_iterations);
    x(:, 1) = start_point;

    for i = 1:2:max_iterations
        
        x(1, i+1) = (4 - x(2, i))/2;
        x(2, i+1) = x(2, i);

        if abs(calc_f_value (x(:, i+1), x1, x2, f) - calc_f_value (x(:, i), x1, x2, f)) < 0.0001

            fprintf("Alternation Minimization converged after %d iterations.\n", i);
            x = x(:, 1:i+1);

            break;
        end
        
        x(2, i+2) = (6 - x(1, i+1))/2;
        x(1, i+2)= x(1, i+1);

   end
   min_point = x(:, end);   
   min_f = double(subs(f, [x1, x2], [min_point(1), min_point(2)]));
   fprintf('Minimum of cost function f with Alternation Minimization: %4f at point = (%.4f, %.4f)\n', min_f,  min_point(1), min_point(2));
   disp('=================================================================================================');



end







function x = gradient_projection(f, x1, x2, start_point, mu, max_iterations)

   x = zeros(2,  max_iterations);
   x(:, 1) = start_point;
   g = gradient(f, [x1, x2]);

   for i = 1:max_iterations
        x(:, i+1) = x(:, i) - mu * double(subs(g, [x1, x2], x(:, i)'));
        x(:, i+1) = x(:, i+1)./ vecnorm(x(:, i+1));

        if abs(calc_f_value (x(:, i+1), x1, x2, f) - calc_f_value (x(:, i), x1, x2, f)) < 0.0001

            fprintf("Gradient Projection converged after %d iterations.\n", i);
            x = x(:, 1:i+1);
  
            break;
        end
   end
   min_point = x(:, end);
   min_f = double(subs(f, [x1, x2], [min_point(1), min_point(2)]));

   fprintf('Minimum of cost function f with Gradient Projection: %4f at point = (%.4f, %.4f)\n', min_f,  min_point(1), min_point(2));
   disp('=================================================================================================');
end

