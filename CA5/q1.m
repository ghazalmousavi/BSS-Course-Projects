clc;
clear;
close all
%%  Part A:

load('hw5-1')

TrainData_class1 = TrainData_class1 - mean(mean(TrainData_class1,3),2);
TrainData_class2 = TrainData_class2 - mean(mean(TrainData_class2,3),2);

Rx_class1 = mean(pagemtimes(TrainData_class1, 'none', TrainData_class1, 'ctranspose'), 3);
Rx_class2 = mean(pagemtimes(TrainData_class2, 'none', TrainData_class2, 'ctranspose'), 3);

[W, D] = eig(Rx_class1, Rx_class2);
W = W ./ vecnorm(W, 2, 1);
[sortedD, idx] = sort(diag(D), 'descend'); 

D = diag(sortedD);
W = W(:, idx);

W_csp1 = W(:, 1);
W_csp30 = W(:, 30);

X1 = TrainData_class1(:, :, 49);
X2 = TrainData_class2(:, :, 49);
t = 0:255;

figure(1)
subplot(211)
plot(t, W_csp1.' * X1, 'r', 'LineWidth', 2);
hold on
plot(t, W_csp1.' * X2, 'b', 'LineWidth', 2);
hold off
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Variance', 'Interpreter', 'latex');
title('Variance of Experiment 49 Using $W_{csp1}$ Over Time', 'Interpreter', 'latex');
xlim([0 256])
legend('class1', 'class2')

subplot(212)
plot(t, W_csp30.' * X1, 'r', 'LineWidth', 2);
hold on
plot(t, W_csp30.' * X2, 'b', 'LineWidth', 2);
hold off
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Variance', 'Interpreter', 'latex');
title('Variance of Experiment 49 Using $W_{csp30}$ Over Time', 'Interpreter', 'latex');
xlim([0 256])
legend('class1', 'class2')

fprintf("Variance of class1 using W_csp1: %.4f\n", var(W_csp1.' * X1));
fprintf("Variance of class2 using W_csp1: %.4f\n", var(W_csp1.' * X2));

disp('========================================');

fprintf("Variance of class1 using W_csp30: %.4f\n", var(W_csp30.' * X1));
fprintf("Variance of class2 using W_csp30: %.4f\n", var(W_csp30.' * X2));

disp('========================================');


%% Part B:

channels = 1:30;
figure(2)
scatter(channels, abs(W_csp1), 40, 'r', 'filled')
xlabel('Channels', 'Interpreter', 'latex');
ylabel('$W_1$', 'Interpreter', 'latex');
title('Absolute Value of $W_{csp1}$ Across Each Channel', 'Interpreter', 'latex', 'FontWeight', 'bold');
grid('on')


figure(3)
scatter(channels, abs(W_csp30), 40, 'b', 'filled')
xlabel('Channels', 'Interpreter', 'latex');
ylabel('$W_2$', 'Interpreter', 'latex');
title('Absolute Value of $W_{csp30}$ Across Each Channel', 'Interpreter', 'latex', 'FontWeight', 'bold');
grid('on')

%% Part C: 

W_csp = [W(:, 1:7), W(:, 24:30)];
number_of_features = size(W_csp, 2);

class1_features = zeros(size(W_csp, 2), size(TrainData_class2, 3));
class2_features = zeros(size(W_csp, 2), size(TrainData_class2, 3));
N =60;

for i=1:size(TrainData_class2, 3)
    x1 = TrainData_class1(:, :, i);
    x2 = TrainData_class2(:, :, i);
    class1_features(:, i) = var(W_csp' * x1, 0, 2);
    class2_features(:, i) = var(W_csp' * x2, 0, 2);
end

mu1 = mean(class1_features, 2);
mu2 = mean(class2_features, 2);

sigma1 = ((class1_features - mu1)*(class1_features - mu1)')./ N;
sigma2 = ((class2_features - mu2)*(class2_features - mu2)')./ N;

M = (mu1 - mu2)*(mu1 - mu2)';
S = sigma1 + sigma2;
[W_LDA, G] = eig(M, S);
[~, index] = max(diag(G));
W_LDA = W_LDA(:, index);
W_LDA = W_LDA ./ vecnorm(W_LDA, 2, 1);

c = (W_LDA'*mu1 + W_LDA'*mu2)/2;

fprintf("Decision boundary: %.4f\n", c);

disp('========================================');

%% Part D & E:

TestData = TestData - mean(mean(TestData,3),2);
test_features = zeros(size(W_csp, 2), size(TestData, 3));

for i=1:size(TestData, 3)
    test_features(:, i) = var(W_csp' * TestData(:, :, i), 0, 2);

    predicted_label = 2 - (W_LDA' * test_features > c);

end


accuracy = sum(TestLabel == predicted_label) / numel(TestLabel);
fprintf('Accuracy: %.2f%%\n', accuracy * 100);

figure(4);
n = 1:size(TestData, 3);
scatter(n, TestLabel, 60, 'blue', 'LineWidth', 1.5)
hold on
scatter(n, predicted_label, 50, 'red', 'filled')
legend('True Labels', 'Predicted Labels')
xlabel('Sample Index', 'Interpreter','latex')
ylabel('Class Labels', 'Interpreter', 'latex')
title('True vs Predicted Labels','Interpreter','latex', 'FontWeight', 'bold')
ylim([0 3])
hold off




