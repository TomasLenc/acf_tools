function tests = test_acf
    tests = functiontests(localfunctions);
end


function test_ap_robust(test_case)
    fs = 100; 
    N = fs * 60; 
    hN = floor(N/2)+1; 
    freq = [0 : hN-1] / N * fs; 
    freq = freq(2:end); 
    mX = aperiodic([3, 1.5], freq); 
    log_pow = log10(mX .^ 2); 
    
    theta1 = fit_aperiodic(freq, log_pow, 'robust', true); 
    theta2 = fit_aperiodic(freq, log_pow, 'robust', false); 
    
    assert(all(abs(theta1 - theta2) < 0.01))
    
%     figure
%     plot(freq, mX)
%     hold on
%     
%     ap = aperiodic(theta1, freq); 
%     ap_pow = 10.^ap; 
%     ap_linear = sqrt(ap_pow); 
%     plot(freq, ap_linear, 'linew', 3)
%     
%     ap = aperiodic(theta2, freq); 
%     ap_pow = 10.^ap; 
%     ap_linear = sqrt(ap_pow); 
%     plot(freq, ap_linear, 'linew', 3)
end


function test_acf_multidim(test_case)

    fs = 100; 
    x = pinknoise(fs * 60)'; 
    x = x - min(x); 
    x_multi_dim = []; 
    x_multi_dim(:, 1, 1, 1, 1, :) = [x; 4.4 + 2 * x]; 
    [acf, lags, ap, mX, freq, ap_par] = get_acf(x_multi_dim, fs, ...
                                               'rm_ap', true, ...
                                               'f0_to_ignore', 1/2.4 ...
                                               ); 
    expected_size = size(x_multi_dim); 
    expected_size(end) = length(lags); 
    assert(all(size(acf) == expected_size))
    verifyNotEqual(test_case, squeeze(mX(1,1,1,1,1,:)), squeeze(mX(2,1,1,1,1,:)))
    verifyNotEqual(test_case, squeeze(acf(1,1,1,1,1,:)), squeeze(acf(2,1,1,1,1,:)))

    [acf_1, ~, ap_1, mX_1, ~, ap_par]  = get_acf(squeeze(x_multi_dim(1,1,1,1,1,:))', fs, ...
                   'rm_ap', true, ...
                   'f0_to_ignore', 1/2.4 ...
                   ); 
    [acf_2, ~, ap_2, mX_2]  = get_acf(squeeze(x_multi_dim(2,1,1,1,1,:))', fs, ...
                   'rm_ap', true, ...
                   'f0_to_ignore', 1/2.4 ...
                   ); 
               
    verifyEqual(test_case, squeeze(acf(1,1,1,1,1,:)), acf_1'); 
    verifyEqual(test_case, squeeze(acf(2,1,1,1,1,:)), acf_2'); 
    
    verifyEqual(test_case, squeeze(ap(1,1,1,1,1,:)), ap_1'); 
    verifyEqual(test_case, squeeze(ap(2,1,1,1,1,:)), ap_2'); 
    
    verifyEqual(test_case, squeeze(mX(1,1,1,1,1,:)), mX_1'); 
    verifyEqual(test_case, squeeze(mX(2,1,1,1,1,:)), mX_2'); 

end


function test_mX_multidim(test_case)

    fs = 100; 
    x = pinknoise(fs * 60)'; 
    x = x - min(x); 
    x_multi_dim = []; 
    x_multi_dim(:, 1, 1, 1, 1, :) = [x; 4.4 + 2 * x]; 
    [~, ~, ~, mX] = get_acf(x_multi_dim, fs, ...
                                               'rm_ap', true, ...
                                               'f0_to_ignore', 1/2.4, ...
                                               'fit_knee', true); 
                                           
    [~, ~, ~, mX_1]  = get_acf(squeeze(x_multi_dim(1,1,1,1,1,:))', fs, ...
                   'rm_ap', true, ...
                   'f0_to_ignore', 1/2.4, ...
                   'fit_knee', true ...
                   ); 
               
    verifyEqual(test_case, squeeze(mX(1,1,1,1,1,:)), mX_1'); 
end


function test_mX_equal_knee_no_knee(test_case)

    fs = 100; 
    x = pinknoise(fs * 60)'; 
    x = x - min(x); 
    [~, ~, ~, mX_1]  = get_acf(x, fs, ...
                               'rm_ap', true, ...
                               'f0_to_ignore', 1/2.4, ...
                               'fit_knee', true ...
                               ); 
    [~, ~, ~, mX_2]  = get_acf(x, fs, ...
                               'rm_ap', true, ...
                               'f0_to_ignore', 1/2.4, ...
                               'fit_knee', false ...
                               ); 
    verifyEqual(test_case, mX_1, mX_2); 
end


function test_acf_plot(test_case)

    f = figure(); 
    ax = axes(f); 
    plot_acf(ax, ...
         randn(1, 100), ...
         [1:100], ...
         'features', struct('z', 0.123456), ...
         'lags_meter_rel', [20, 60], ...
         'lags_meter_unrel', [10, 30, 70], ...
         'min_lag', 2, ...
         'max_lag', 99, ...
         'prec', 1000); 
    close(f); 
    
end







