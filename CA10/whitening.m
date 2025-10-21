function [W, Z] = whitening(X)
    R = X * X';
    [V, D] = eig(R);
    [sortedD, idx] = sort(diag(D), 'descend');
    D = diag(sortedD);
    V = V(:, idx);
    W = D^(-0.5) * V';
    Z = W * X;
end
