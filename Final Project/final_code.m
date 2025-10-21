clc; clear; close all
load('dataset\subj_15.mat')
 
N  = size(data, 2);
max_itr = 100;

n_trials1 = size(data{1}, 3);
n_trials2 = size(data{2}, 3);
n_trials3 = size(data{3}, 3);
n_trials4 = size(data{4}, 3);

n_features = 6;
channels   = [39, 15, 44, 48, 20, 24, 49, 53, 25, 57, 21, 50, 54, 29, 31, 26, 60, 22, 55, 51, 27, 30, 32, 59, 58, 23, 47, 43, 19, 52, 28];


train_cm   = zeros(N, N);
test_cm    = zeros(N, N);
bands = [[12, 45]; [5, 30]; [2 , 35]];


for itr = 1 : max_itr
    thresholds = zeros(1, N - 1);
    W_ldas = cell(1, N - 1);
    W_csps = cell(1, N - 1);

    r = [randi(n_trials1), randi(n_trials2), randi(n_trials3), randi(n_trials4)];
    
    
    train_data_class1 = data{1}(channels, :, [1:r(1)-1, r(1)+1:n_trials1]);
    train_data_class2 = data{2}(channels, :, [1:r(2)-1, r(2)+1:n_trials2]);
    train_data_class3 = data{3}(channels, :, [1:r(3)-1, r(3)+1:n_trials3]);
    train_data_class4 = data{4}(channels, :, [1:r(4)-1, r(4)+1:n_trials4]);

    % Combine class 1, class 2, class 3
    train_data_class123 = cat(3, train_data_class1, train_data_class2, train_data_class3);
    train_data_class123_filtered = filtering(train_data_class123, bands(1, :));
    train_data_class4_filtered = filtering(train_data_class4, bands(1, :));
    
    % Classifier 1: class4 vs class1, class2 and class3
    N1 = n_trials1 + n_trials2 + n_trials3;
    [thresholds(1), W_ldas{1}, W_csps{1}] = train_binary_eeg_classifier(train_data_class4_filtered, train_data_class123_filtered, n_features, n_trials4, N1);
    
    % Classifier 2: class3 vs class1 and class2
    N2 = n_trials1 + n_trials2;
    train_data_class12 = cat(3, train_data_class1, train_data_class2);
    train_data_class12_filtered = filtering(train_data_class12, bands(2, :));
    train_data_class3_filtered = filtering(train_data_class3, bands(2, :));

    [thresholds(2), W_ldas{2}, W_csps{2}] = train_binary_eeg_classifier(train_data_class3_filtered, train_data_class12_filtered, n_features, n_trials3, N2);
    
    % Classifier 3: class2 vs class1
    train_data_class2_filtered = filtering(train_data_class2, bands(3, :));
    train_data_class1_filtered = filtering(train_data_class1, bands(3, :));

    [thresholds(3), W_ldas{3}, W_csps{3}] = train_binary_eeg_classifier(train_data_class2_filtered, train_data_class1_filtered, n_features, n_trials1, n_trials2);
    
    % Confusion Matrix (Training Data):
    predicted_labels_class1 = zeros(1, n_trials1 - 1);
    for k = 1 : n_trials1 - 1
        sample = train_data_class1(:, :, k); 
        predicted_labels_class1(k) = predict_label(sample, W_csps, W_ldas, thresholds, bands);
    end
    train_cm = confusion_matrix(train_cm, ones(1, n_trials1 - 1), predicted_labels_class1);

    predicted_labels_class2 = zeros(1, n_trials2 - 1);
    for k = 1 : n_trials2 - 1
        sample = train_data_class2(:, :, k); 
        predicted_labels_class2(k) = predict_label(sample, W_csps, W_ldas, thresholds, bands);
    end 
    train_cm = confusion_matrix(train_cm, 2.*ones(1, n_trials2 - 1), predicted_labels_class2);

    predicted_labels_class3 = zeros(1, n_trials3 - 1);
    for k = 1 : n_trials3 - 1
        sample = train_data_class3(:, :, k); 
        predicted_labels_class3(k) = predict_label(sample, W_csps, W_ldas, thresholds, bands);
    end
    train_cm = confusion_matrix(train_cm, 3.*ones(1, n_trials3 - 1), predicted_labels_class3);

    predicted_labels_class4 = zeros(1, n_trials4 - 1);
    for k = 1 : n_trials4 - 1
        sample = train_data_class4(:, :, k); 
        predicted_labels_class4(k) = predict_label(sample, W_csps, W_ldas, thresholds, bands);
    end
    train_cm = confusion_matrix(train_cm, 4.*ones(1, n_trials4 - 1), predicted_labels_class4);

    % Confusion Matrix (Validation Data):
    valid_dataset = cat(4, data{1}(channels, :, r(1)), data{2}(channels, :, r(2)), data{3}(channels, :, r(3)), data{4}(channels, :, r(4)));
    val_predicted_labels = zeros(1, N);
    
    for i = 1 : size(valid_dataset, 4)
        sample = valid_dataset(:, :, :, i); 
        val_predicted_labels(i) = predict_label(sample, W_csps, W_ldas, thresholds, bands);
    end
    
    true_labels = 1:N;    
    test_cm = confusion_matrix(test_cm, true_labels, val_predicted_labels);
end


%% Compute Accuracy:
 
fprintf('<strong>Validation Dataset Accuracy:</strong>\n');

class_accuracies_test = diag(test_cm) ./ sum(test_cm, 2);
for i = 1 : length(class_accuracies_test)
    fprintf('Class %d: %.2f%%\n', i, class_accuracies_test(i) * 100);
end

fprintf('%s\n', repmat('=', 1, 25));

fprintf('<strong>Train Dataset Accuracy:</strong>\n');  

class_accuracies_train = diag(train_cm) ./ sum(train_cm, 2);
for i = 1 : length(class_accuracies_test)
    fprintf('Class %d: %.2f%%\n', i, class_accuracies_train(i) * 100);
end

