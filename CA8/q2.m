clc; clear; close all

%% Frame 2 * N

N = 2:50;
mus = zeros(1, 49);
for n = 1 : length(N)
    theta = (0:N(n)-1) * pi/N(n); 
    D = [cos(theta); sin(theta)];
    mus(n) = mutual_coherence(D);
end

figure;
plot(N, mus, 'b', 'LineWidth', 1.5)
xlabel('N', 'Interpreter', 'latex')
ylabel('Mutual Coherence', 'Interpreter', 'latex')
title('Mutual Coherence vs Number of Vectors $N$', 'Interpreter', 'latex')
grid on


%%  Frame 3*N:

N_values = 2:50;
mus = zeros(1, length(N_values));

for n = 1:length(N_values)
    N = N_values(n);
    if N == 2
        D = eye(3, 2);
    elseif N == 3
        D = eye(3);
    else
        indices = 0:N-1;
        phi_angle = pi * (3 - sqrt(5));
        theta = indices * phi_angle;
        z = linspace(1 - 1/(2*N), 1/(2*N) - 1, N);
        radius = sqrt(1 - z.^2);
        x = radius .* cos(theta);
        y = radius .* sin(theta);
        D = [x; y; z];
        D = D ./ vecnorm(D);
    end
    mus(n) = mutual_coherence(D);
end

figure;
plot(N_values, mus, 'r', 'LineWidth', 1.5)
xlabel('N', 'Interpreter', 'latex')
ylabel('Mutual Coherence', 'Interpreter', 'latex')
title('Mutual Coherence vs Number of Vectors $N$ in 3D', 'Interpreter', 'latex')
grid on

%% Functions:

function mu = mutual_coherence(D)
    G = abs(D' * D);
    G(logical(eye(size(G)))) = 0;
    mu = max(G(:));
end

