function k = kurtosis(y)
    m4 = mean(y.^4);     % 4th moment
    m2 = mean(y.^2);     % 2nd moment
    k = m4 - 3*m2^2;
end
