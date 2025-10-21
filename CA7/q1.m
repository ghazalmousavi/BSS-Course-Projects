clc; clear; close all
load('hw7.mat')

N0 = 3;
N = size(D, 2);
%% Part a:
S_subset = SubsetSelection(D, x, N0, N);

%% Part b:
S_norm2 = pinv(D)*x;

%% Part c:
S_MP1 = MatchingPursuit_KnownSparsityLevel(D, x, N0, N, 1);
S_MP2 = MatchingPursuit(D, x, N);

%% Part d:
S_OMP1 = OMP_KnownSparsityLevel(D, x, N0, N, 1);
S_OMP2 = OrthogonalMatchingPursuit(D, x, N);

%% Part e:
S_BP = BasisPursuit(D, x, N);

%% Part f:
S_IRLS = IRLS(D, x, N);

%%  Functions:

function S_hat = SubsetSelection(D, X, N0, N)
 
    tic;
    vals = 1 : N;
    combs = nchoosek(vals, N0);
    
    error = zeros(size(combs, 1), 1);
    S_hat = zeros(N, 1);

    for i = 1: length(combs)
        idx = combs(i, :);
        s_hat = pinv(D(:, idx))*X;   
        error(i) = norm(D(:, idx)*s_hat - X, 'fro');
    end

    [~, pos] = min(error);
    subset = combs(pos, :);

    D_subset = D(:, subset);      
    S_hat(subset) = pinv(D_subset)* X; 
    
    elapsedTime = toc;  
    fprintf('%s\n', repmat('=', 1, 60));
    fprintf('<strong>Subst Selection:</strong>\n');

    
    PrintResults(subset, S_hat, error(pos), elapsedTime)

end

function [S_hat, subsets] = MatchingPursuit_KnownSparsityLevel(D, x, N0, N, flag)
    tic;
    xr = x;
    subsets = zeros(1, N0);
    S_hat = zeros(N, 1);

    for i = 1:N0

        projections = D' * xr;

        [~, idx] = max(abs(projections));

        subsets(i) = idx;

        d_star = D(:, idx);

        a = xr'*d_star;
        xr = xr - a*d_star;

        S_hat(idx) = projections(idx);
    end
    
    error = norm(D*S_hat - x);
    elapsedTime = toc;
    if flag == 1
        fprintf('<strong>Matching Pursuit Kown N0 = 3:</strong>\n');
        PrintResults(subsets, S_hat, error, elapsedTime)
    end

end


function S_hat = MatchingPursuit(D, x, N)
    tic; 
    S_hat = zeros(N, 1);
    threshold = 0.05;
    for i = 1 : N 
        [S_hat, subsets] = MatchingPursuit_KnownSparsityLevel(D, x, i, N, 0);
        if norm(x - D*S_hat)/norm(x) < threshold
            elapsedTime = toc;
            error = norm(x - D*S_hat);
            fprintf('<strong>Matching Pursuit</strong>\n');
            fprintf('Sparsity Level is: N0 = %d.\n', i);
            PrintResults(subsets, S_hat, error, elapsedTime)
            break
        end
    end
end


function [S_hat, subsets] = OMP_KnownSparsityLevel(D, x, N0, N, flag)
    tic;
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

    
    S_hat(subsets, :) = pinv(D(:, subsets))*x;
    if flag == 1
        elapsedTime = toc;
        error = norm(x - D*S_hat);
        fprintf('<strong>Orthogonal Matching Pursuit Kown N0 = 3 </strong>\n');
        PrintResults(subsets, S_hat, error, elapsedTime)
    end

end


function S_hat = OrthogonalMatchingPursuit(D, x, N)
    tic;
    S_hat = zeros(N, 1);
    threshold = 0.05;
    for i = 1 : N 
        [S_hat, subsets] = OMP_KnownSparsityLevel(D, x, i, N, 0);
        if norm(x - D*S_hat)/norm(x) < threshold
            elapsedTime = toc;
            error = norm(x - D*S_hat);
            fprintf('<strong>Orthogonal Matching Pursuit</strong>\n');
            fprintf('Sparsity Level is: N0 = %d.\n', i);
            PrintResults(subsets, S_hat, error, elapsedTime)
            break
        end
    end
end

function S_hat = BasisPursuit(D, x, N)
    f = ones(2*N, 1);
    Aeq = [D, -D];
    beq = x; 
    lb = zeros(2*N, 1);
    tic;
    S_hat = linprog(f, [], [], Aeq, beq, lb, []);
    S_hat = S_hat(1:N) - S_hat(N + 1:2*N);
    
    elapsedTime = toc;
    error = norm(D*S_hat -x, 'fro');
    subset = find(S_hat);
    fprintf('%s\n', repmat('=', 1, 60));
    fprintf('<strong>Basis Pursuit:</strong>\n'); 
    PrintResults(subset, S_hat, error, elapsedTime)
end

function S_hat = IRLS(D, x, N)
    w = ones(N, 1);
    S_hat = zeros(N, 1);
    max_itr = 15;
    tic;
    for i=1:max_itr
        w = diag(w);
        S_hat = (w^(-1)*D.')*(D*w^(-1)*D')^(-1)*x;
        w = 1 ./ abs(S_hat);
    end
    fprintf('<strong>IRLS:</strong>\n'); 
    elapsedTime = toc;
    S_hat(abs(S_hat) < 1e-10) = 0;
    subset = find(abs(S_hat) > 0.1);
    error = norm(D*S_hat - x);
    PrintResults(subset, S_hat, error, elapsedTime)
end

function PrintResults(subset, s_hat, error, elapsedTime)

    fprintf('Runtime: %.3e\n', elapsedTime);
    fprintf('Error is: %.2e\n', error);
    fprintf('%s\n', repmat('-', 1, 30));

    
    fprintf('%10s | %15s\n', 'Index', 'S_hat(Value)');
    fprintf('%s\n', repmat('-', 1, 30));
    for i = 1:length(subset)
        fprintf('%10d | %15.6f\n', subset(i), s_hat(subset(i)));
    end
    fprintf('%s\n', repmat('-', 1, 30));
    
    fprintf('%s\n', repmat('=', 1, 60));

end


