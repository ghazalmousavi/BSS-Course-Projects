function [threshold, W_LDA, W_csp] = train_binary_eeg_classifier(train_dataset_class1, train_dataset_class2, n_features, N1, N2)

    x1 = gpuArray(train_dataset_class1);              
    x2 = gpuArray(train_dataset_class2);               

    Rx1 = mean(pagemtimes(x1, 'none', x1, 'ctranspose'), 3);
    Rx2 = mean(pagemtimes(x2, 'none', x2, 'ctranspose'), 3);
    
    
    W_csp = extract_features_CSP(gather(Rx1), gather(Rx2), n_features);
    W_csp = gpuArray(W_csp);
    W_csp_T = permute(W_csp.', [1, 2]);  
    
    projected1 = pagefun(@mtimes, W_csp_T, x1); 
    class1_features = var(projected1, 0, 2);
    class1_features = squeeze(class1_features); 
    
    
    projected2 = pagefun(@mtimes, W_csp_T, x2); 
    class2_features = var(projected2, 0, 2);
    class2_features = squeeze(class2_features); 
    
    class1_features = gather(class1_features);
    class2_features = gather(class2_features);
    
    [threshold, W_LDA] = linear_discriminant_analysis(class1_features, class2_features, N1, N2);
     
end