function label = predict_label(sample, W_csps, W_ldas, thresholds, bands)

    % Classifier 1: Class 4 vs Classes 1,2,3
    sample1 = filtering(sample, bands(1,:));
    pred = classify_test_data(sample1, W_csps{1}, W_ldas{1}, thresholds(1));
    if pred == 1
        label = 4;
        return;
    end
    
    % Classifier 2: Class 3 vs Classes 1,2
    sample2 = filtering(sample, bands(2,:));
    pred = classify_test_data(sample2, W_csps{2}, W_ldas{2}, thresholds(2));
    if pred == 1
        label = 3;
        return;
    end
    
    % Classifier 3: Class 2 vs Class 1
    sample3 = filtering(sample, bands(3,:));
    pred = classify_test_data(sample3, W_csps{3}, W_ldas{3}, thresholds(3));
    if pred == 1
        label = 2;
    else
        label = 1;
    end
    
end
