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




function test_acf_f0_ignore(test_case)

%     fs = 100; 
%     trial_dur = 60; 
%     noise = pinknoise(fs * trial_dur)'; 
%     noise = zscore(noise); 
%     
%     p0 = 0.8; 
%     s = zeros(size(noise)); 
%     event = ones(1, round(0.2 * fs)); 
%     onsets = [p0 : p0 : trial_dur-p0];
%     for i=1:length(onsets)
%         idx = round(onsets(i) * fs); 
%         s(idx+1 : idx+length(event)) = event; 
%     end
%     x = s + noise;
    
    load('test_acf_f0_ignore_data.mat'); 
    
    % wihtout 1/f subtraction 
    [acf_raw, lags, ap, mX, freq, ap_par] = get_acf(...
                                       x, fs, ...
                                       'rm_ap', false ...
                                       );     
                                   
    % with 1/f subtraction but keeping all frequencies when calculating ACF
    [acf_subtr, lags, ap, mX, freq, ap_par] = get_acf(...
                                       x, fs, ...
                                       'rm_ap', true, ...
                                       'f0_to_ignore', 1/p0, ...
                                       'only_use_f0_harmonics', false ...
                                       ); 

    % with 1/f subtraction but only keeping f0 harmonics when calculating ACF
    [acf_subtr_onlyF0, lags, ap, mX, freq, ap_par] = get_acf(...
                                       x, fs, ...
                                       'rm_ap', true, ...
                                       'f0_to_ignore', 1/p0, ...
                                       'only_use_f0_harmonics', true ...
                                       );                                   
                                   
    acf_raw = zscore(acf_raw); 
    acf_subtr = zscore(acf_subtr); 
    acf_subtr_onlyF0 = zscore(acf_subtr_onlyF0); 
    
    idx_half = floor(length(acf_raw) / 2); 
    
    slope_index_raw = ...
        mean(acf_raw(1 : idx_half)) - mean(acf_raw(idx_half : end));
    
    slope_index_subtr = ...
        mean(acf_raw(1 : idx_half)) - mean(acf_subtr(idx_half : end));
    
    slope_index_subtr_onlyF0 = ...
        mean(acf_subtr_onlyF0(1 : idx_half)) - mean(acf_subtr_onlyF0(idx_half : end));
    
    assert(slope_index_raw > slope_index_subtr); 
    assert(slope_index_subtr > slope_index_subtr_onlyF0); 
    
%     figure
%     plot(acf_raw); 
%     hold on 
%     plot(acf_subtr); 
%     plot(acf_subtr_onlyF0); 

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


function test_single_input_ap_fit_error(test_case)

    fs = 128; 
    x = single(randn(1, fs * 22)); 
                                                          
    verifyError(test_case, ...
        @() get_acf(x, fs, 'rm_ap', true),...
        'fit_aperiodic:inputNotDouble');
end


function test_ap_flims(test_case)

    % create 1/f signal 
    fs = 128; 
    N = fs * 100; 
    hN = floor(N/2)+1; 
    x = get_colored_noise(N, fs, -1); 
    x = zscore(x); 
    
    % add power at high frequencies 
    freq = [0 : hN-1] / N * fs; 
    idx = [dsearchn(freq', 30) : dsearchn(freq', fs/2)];
    X = fft(x); 
    X(idx) = X(idx) * 3; 
    X_new = [X(1:hN), flip(X(2:hN-1))];
    x = real(ifft(X_new));
                                                          
    [acf, lags, ap1, mX, freq] = get_acf(x, fs, 'rm_ap', true, ...
                                         'ap_fit_flims', [0.01, fs/2]);  
                                     
    [acf, lags, ap2, mX, freq] = get_acf(x, fs, 'rm_ap', true, ...
                                         'ap_fit_flims', [0.01, 29]);  
    
%     figure
%     plot(freq, mX); 
%     hold on 
%     plot(freq, ap1, ':', 'linew', 2); 
%     plot(freq, ap2, ':', 'linew', 2); 
%     
    % taking into the account the higher frequencies results in a flatter 1/f
    % estimate
    assert(all(ap1(2:10) < ap2(2:10))); 
    
end


function test_acf_flims(test_case)

    % create 1/f signal 
    fs = 200; 
    N = fs * 100; 
    hN = floor(N/2)+1; 
    x = get_colored_noise(N, fs, -1); 
    x = zscore(x); 
    
    % add power at a frequency
    t = [0 : N-1] / fs; 
    s = cos(2 * pi * t * 10); 
    x = x + 2*s; 
                                                          
    [acf1, lags] = get_acf(x, fs, 'rm_ap', true, ...
                                         'acf_flims', [0, 10-fs/N]);  
                                     
    [acf2, lags] = get_acf(x, fs, 'rm_ap', true, ...
                                     'acf_flims', [0, 10]);  
    
    [acf3, lags] = get_acf(x, fs, 'rm_ap', true, ...
                                     'acf_flims', [0, fs/2]);  

    [acf4, lags] = get_acf(x, fs, 'rm_ap', true);  

%     figure
%     plot(lags, acf1); 
%     hold on 
%     plot(lags, acf2); 
%     plot(lags, acf3); 
%     xlim([0, 2]); 
    
    assert(max(abs(acf1 - acf2)) > 1e4 * eps * min(abs([acf1 - acf2]))); 
    assert(max(abs(acf3 - acf4)) <= 1e4 * eps * min(abs([acf3 - acf4]))); 
    
    
end





