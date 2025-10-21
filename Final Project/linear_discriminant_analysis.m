function [threshold, W_LDA] = linear_discriminant_analysis(class1_features, class2_features, N1, N2)

    mu1 = mean(class1_features, 2);
    mu2 = mean(class2_features, 2);
    
    sigma1 = ((class1_features - mu1)*(class1_features - mu1)')./ N1;
    sigma2 = ((class2_features - mu2)*(class2_features - mu2)')./ N2;
    
    M = (mu1 - mu2)*(mu1 - mu2)';
    S = sigma1 + sigma2;
    [W_LDA, G] = eig(M, S);
    [~, index] = max(diag(G));
    W_LDA = W_LDA(:, index);
    W_LDA = W_LDA ./ vecnorm(W_LDA, 2, 1);
    
    threshold = (W_LDA'*mu1 + W_LDA'*mu2)/2;

end