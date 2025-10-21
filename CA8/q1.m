clc; clear; close all;

load('hw8-part1.mat')
N = size(D, 2);
N0 = 2;

%% Part 1:
f = ones(2 * N, 1);
Aeq = [D, -D];
beq = x; 
lb = zeros(2 * N, 1);
s_hat_bp = linprog(f, [], [], Aeq, beq, lb, []);
s_hat_bp = s_hat_bp(1 : N) - s_hat_bp(N + 1 : 2 * N);

%% Part 2:
s_hat_mp = MP(D, x, N0, N);

%% Part 3:
beq_noisy = x_noisy; 
lb = zeros(2 * N, 1);
s_hat_bp_noisy = linprog(f, [], [], Aeq, beq_noisy, lb, []);
s_hat_bp_noisy = s_hat_bp_noisy(1 : N) - s_hat_bp_noisy(N + 1 : 2 * N);


%% Part 4:
s_hat_lasso = rand(N, 1);
lambda = 0.7;
max_itr = 200;

for itr = 1 : max_itr
    s_hat_lasso = LASSO(x_noisy, D, N, lambda, s_hat_lasso);
end


%% LASSO Function:
function s_hat = LASSO(x, D, N, lambda, s_hat)
    for i = 1:N
        idx = true(N,1); 
    
        idx(i) = false; 
    
        r = x - D(:, idx)*s_hat(idx);
    
        rho = r.' * D(:, i);
        
        if rho > lambda / 2
            s_hat(i) = rho - lambda / 2;
        elseif rho < -lambda / 2
            s_hat(i) = rho + lambda / 2;
        else
            s_hat(i) = 0;
        end
    
    end
end

%% MP Function:

function [S_hat, subsets] = MP(D, x, N0, N)
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

end
    