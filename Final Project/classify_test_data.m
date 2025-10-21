function predicted_label = classify_test_data(test_dataset, W_csp, W_LDA, threshold)

    test_dataset = gpuArray(test_dataset);
    W_csp = gpuArray(W_csp);
    W_LDA = gpuArray(W_LDA);
    threshold = gpuArray(threshold);

    projected = pagefun(@mtimes, W_csp.', test_dataset);  

    test_features = var(projected, 0, 2);  
    test_features = squeeze(test_features); 

    scores = W_LDA.' * test_features;   

    predicted_label = gather(1 + (scores > threshold));  

end