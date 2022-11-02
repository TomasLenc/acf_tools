function tests = test_fft_normalization
    tests = functiontests(localfunctions);
end

function test_noise_subtraction(test_case)
    mX = zeros(1, 10);  
    idx_ones = [5]; 
    mX(idx_ones) = 1; 
    out = subtract_noise_bins(mX, 2, 5); 
    out = out(idx_ones); 
    expected = 1; 
    verifyEqual(test_case, out, expected)
end

function test_fft_interpolation(test_case)
    mX = zeros(1, 20);  
    idx_ones = [3, 17]; 
    mX(idx_ones) = 1; 
    expected = zeros(1, 20);  
    out = interp_noise_bins(mX, idx_ones, 2, 5); 
    verifyEqual(test_case, out, expected)
end

function test_noise_subtraction_multidim(test_case)
    mX = zeros(3, 4, 20);  
    idx_ones = {2, 3, [3, 13]}; 
    mX(idx_ones{:}) = 1; 
    out = subtract_noise_bins(mX, 2, 5); 
    out = squeeze(out(2, 3, :))'; 
    expected = subtract_noise_bins(squeeze(mX(2, 3, :))', 2, 5);
    verifyEqual(test_case, out, expected)
end

function test_fft_interpolation_multidim(test_case)
    mX = zeros(3, 4, 20);  
    idx_ones = {2, 3, [3, 13]}; 
    mX(idx_ones{:}) = 1; 
    out = interp_noise_bins(mX, [3, 13], 2, 5); 
    expected = zeros(size(mX));  
    verifyEqual(test_case, out, expected)
end