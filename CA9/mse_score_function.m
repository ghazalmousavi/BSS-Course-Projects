function sai = mse_score_function(ym)
    
    n = length(ym);

    kernel = @(y) [ones(1, n); y; y.^2; y.^3; y.^4; y.^5];
    kernel_prime = @(y) [zeros(1, n); ones(1, n); 2*y; 3*y.^2; 4*y.^3; 5*y.^4];

    k = kernel(ym);
    K = k*k.';
    kp = kernel_prime(ym);

    theta = inv(K / n) * mean(kp, 2);
    sai = theta.' * k;
end