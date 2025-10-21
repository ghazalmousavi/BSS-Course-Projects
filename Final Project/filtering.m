function filtered_data = filtering(data, band)

    fs = 2.4e3;
    gpu_data = gpuArray(data);
    [~, T, ~] = size(data);

    f = (0:T-1) * fs / T;  
    filter = (f >= band(1) & f <= band(2));  
    filter = reshape(filter, [1, T, 1]);  
    
    % FFT along time dimension
    data_fft = fft(gpu_data, [], 2);  
    
    data_fft_filtered = data_fft .* gpuArray(filter);  
    
    filtered_data = real(ifft(data_fft_filtered, [], 2));  
    
    filtered_data = gather(filtered_data);

end