function W_csp = extract_features_CSP(R1, R2, num_features)

    [U, D] = eig(R1, R2);
    U = U ./ vecnorm(U, 2, 1);
    [~, idx] = sort( diag(D), 'descend');
    U = U(:, idx);

    half = num_features / 2;
    W_csp = [U(:, 1:half), U(:, end - half + 1:end)];
end