clc; clear; close all;

load("hw8-part3.mat")
rng(8)
%% Parameters;
N0 = 2;
[M, N] = size(D);
T = size(X, 2);

%% Part 1: Mutual Coherence

G = D.'*D;
G(logical(eye(size(G)))) = 0;
mu = max(abs(G(:)));
fprintf('Mutual Coherence: %f\n', mu);
fprintf('%s\n', repmat('-', 1, 60));

%% Part 2: MOD

max_itr = 100;
threshold = 1e-6;

D_MOD = randn(M, N);
D_MOD = D_MOD ./ vecnorm(D_MOD);   
S_MOD = zeros(N, T);
Error_MOD = zeros(1, max_itr);

tic;
for itr = 1:max_itr
    % --- Sparse Recovery Step (OMP) ---
    for t = 1:T
        S_MOD(:, t) = OMP(D_MOD, X(:, t), N0, N);
    end

    % --- Dictionary Update Step ---
    D_MOD = X * pinv(S_MOD);
    D_MOD = D_MOD ./ vecnorm(D_MOD); 

    Error_MOD(itr) = norm(X - D_MOD * S_MOD, 'fro') ^ 2 / (norm(X, 'fro')^2);

    if itr > 1 && abs(Error_MOD(itr) - Error_MOD(itr - 1)) < threshold
        Error_MOD(itr: end) = Error_MOD(itr); 
        fprintf('MOD Converged after %d iterations in %.4f seconds\n', itr, toc);
        break;
    end

    if itr == max_itr
        fprintf('MOD Converged after %d iterations in %.4f seconds\n', max_itr, toc);
    end

end

%% Part 2: K-SVD  

D_KSVD = randn(M, N);
D_KSVD = D_KSVD ./ vecnorm(D_KSVD);   
S_KSVD = randn(N, T);
Error_KSVD = zeros(1, max_itr);

tic;
for itr = 1:max_itr

    % --- Sparse Recovery Step using OMP ---
    for t = 1:T
        S_KSVD(:, t) = OMP(D_KSVD, X(:, t), N0, N);
    end

    % --- Dictionary Update Step ---
    for k = 1:N
        idx = true(N, 1);   
        idx(k) = false;     

        % Residual without atom k
        X_r = X - D_KSVD(:, idx) * S_KSVD(idx, :);

        nonzero_cols = S_KSVD(k, :) ~= 0;
        X_r = X_r(:, nonzero_cols);

        [U, Sigma, V] = svd(X_r, 'econ');

        D_KSVD(:, k) = U(:, 1);
        D_KSVD = D_KSVD ./ vecnorm(D_KSVD);   

        S_KSVD(k, nonzero_cols) = Sigma(1, 1) * V(:, 1)';
    end

    Error_KSVD(itr) = norm(X - D_KSVD * S_KSVD, 'fro') ^ 2/ norm(X, 'fro')^2;

    if itr > 1 && abs(Error_KSVD(itr) - Error_KSVD(itr-1)) < threshold
        Error_KSVD(itr: end) = Error_KSVD(itr); 
        fprintf('K-SVD Converged after %d iterations in %.4f seconds\n', itr, toc);
        break;
    end

    if itr == max_itr
        fprintf('K-SVD Converged after %d iterations in %.4f seconds\n', max_itr, toc);
    end

end

%% Part 3: Convergence Plot:
figure;
plot(1:max_itr, Error_MOD, 'LineWidth', 1.5)
hold on
plot(1:max_itr, Error_KSVD, 'LineWidth', 1.5)
title('Reconstruction Error', 'Interpreter', 'latex')
xlabel('Iteration', 'Interpreter', 'latex')
ylabel('Normalized Error', 'Interpreter', 'latex')
legend('MOD', 'K-SVD', 'Interpreter', 'latex')
grid on
hold off
 
%% Part 5: Successful Recovery Rate:

fprintf('%s\n', repmat('-', 1, 60));

srr_MOD = successful_recovery_rate(D_MOD, D);
fprintf('MOD Successful Recovery Rate: %.2f%%\n', srr_MOD * 100)

srr_KSVD = successful_recovery_rate(D_KSVD, D);
fprintf('K-SVD Successful Recovery Rate: %.2f%%\n', srr_KSVD * 100)

fprintf('%s\n', repmat('-', 1, 60));

%% Part 6: 

S_hat_MOD = align_using_correlation(S, S_MOD);
error_MOD = norm(S_hat_MOD - S, 'fro')^2 / norm(S_hat_MOD, 'fro')^2;
fprintf('MOD Error:   %.6f\n', error_MOD);

S_hat_KSVD = align_using_correlation(S, S_KSVD);
error_KSVD = norm(S_hat_KSVD - S, 'fro')^2 / norm(S_hat_KSVD, 'fro')^2;
fprintf('K-SVD Error: %.6f\n', error_KSVD);

fprintf('%s\n', repmat('-', 1, 60));


%% Functions: 

function srr = successful_recovery_rate(Dhat, D)
    D_temp = D;
    total_matched = 0;
    N = size(D, 2);
    for i = 1:N
        corr_vals = abs(corr(Dhat(:, i), D_temp));
        idx = find(corr_vals > 0.98);
        total_matched = total_matched + numel(idx);
    
        D_temp(:, idx) = []; 
    
    end
    srr = total_matched / N;
end


function S_hat_aligned = align_using_correlation(S, S_hat)
    [~, T] = size(S);
    S_temp = S;
    used_indices = false(1, T);
    S_hat_aligned = zeros(size(S));

    for i = 1:T
        corr_vals = abs(corr(S_hat(:, i), S_temp));
        [~, best_idx] = max(corr_vals);
        real_idx = find(~used_indices);
        matched_col = real_idx(best_idx);
        used_indices(matched_col) = true;
        scale = (S(:, matched_col)' * S_hat(:, i)) / (S_hat(:, i)' * S_hat(:, i));
        S_hat_aligned(:, matched_col) = S_hat(:, i) * scale;
        S_temp(:, best_idx) = [];
    end
end


function S_hat= OMP(D, x, N0, N)
    xr = x;
    subsets = zeros(1, N0);
    S_hat = zeros(N, 1);
    D_selected = [];
    
    for i = 1:N0
        projections = D' * xr;
        [~, idx] = max(abs(projections));
        subsets(i) = idx;
        D_selected = [D_selected, D(:, idx)];
        coeffs = pinv(D_selected) * x;
        xr = x - D_selected * coeffs;
    end

    S_hat(subsets) = pinv(D(:, subsets)) * x;
end

