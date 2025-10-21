clc;
clear;
close all;
load('hw5-2')

%% Defining constants:

[~, T, N] = size(data);
fs = 250;
t = 0:1/fs:(T-1)/fs;
ro = zeros(N, length(freq));
predicted_freq = zeros(size(label));

%% Template Matrices:

X = cell(1, length(freq));  

for i = 1:length(freq)
    X{i} = create_template_matrix(freq(i), t);
end

%% Feature Extraction and Prediction

for i = 1: N
    for j = 1:length(freq)
        ro(i, j) = extract_features(X{j}, data(:, :, i));  
    end

    [~, idx] = max(ro(i, :));     
    predicted_freq(i) = freq(idx);   
end


%% Functions:

function ro = extract_features(X, Y)

    Rx =X*X.';
    Ry = Y*Y.';
    Rxy = X*Y.';
    Ryx = Y*X.';

    sigma1 = Rx^(-0.5)*Rxy*Ry^(-1)*Ryx*Rx^(-0.5);
    sigma2 = Ry^(-0.5)*Ryx*Rx^(-1)*Rxy*Ry^(-0.5);
    
    [c, D] = eig(sigma1);
    [~, index1] = max(diag(D));
    c = c(:, index1);
    
    a = Rx^(-0.5)*c;
    

    [d, G] = eig(sigma2);
    [~, index2] = max(diag(G));
    d = d(:, index2);

    b = Ry^(-0.5)*d;

    ro = abs(a'*X*Y.'*b)/((a'*Rx*a)^0.5*(b'*Ry*b)^0.5);

end



function template_matrix = create_template_matrix(f, t)
    
    max_harmonic_freq = 40;
    n_harmonics = floor(max_harmonic_freq / f);
    template_matrix = zeros(2 * n_harmonics, size(t, 2));

    for n = 1:n_harmonics
        template_matrix(2 * n - 1, :) = sin(2 * pi * n * f * t);
        template_matrix(2 * n, :) = cos(2 * pi * n * f * t);
    end
end



