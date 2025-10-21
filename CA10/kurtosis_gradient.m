function [gradient_value, kurt_val] = kurtosis_gradient(y, b, Z)
    kurt_val = kurtosis(y);
    m = mean(Z .* y.^3, 2);  
    gradient_value = sign(kurt_val) * (m - 3 * b);
end
