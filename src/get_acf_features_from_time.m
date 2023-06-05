function feat = get_acf_features_from_time(x, fs, lags_meter_rel, lags_meter_unrel)
% Calculate prominence of meter-related lags from time-domain input. This is an
% alternative to get_acf_features() function that requires pre-calculating the
% whole autocorrelatino function. The problem is that we can't get normalized
% correlation coefficients from FFT-implementation of autocorrelation. Given
% that we're only interested in the value of autocorrelation at a few apriori
% specified lags, we don't need the whole function and can simply calculate
% pearson correlation with time-shifted versions of the input signal. 
% 
% Parameters
% ----------
% x : array_like, shape=[..., time]
%     Input x with time as the last dimension. 
% fs : int
%     Sampling rate. 
% lags_meter_rel : array_like
%     1-D array of time lags that are meter related. 
% lags_meter_unrel : array_like
%     1-D array of time lags that are meter unrelated. 
%
% Returns
% -------
% feat : struct
%     Structure with calculated features. 
%     

shape = size(x); 

% get indices for lags of interest
lags_all = [lags_meter_rel, lags_meter_unrel];
lags_all_idx = round(lags_all * fs); 

% get pearson correlation for each lag
r = nan([shape(1:end-1), length(lags_all)]); 
for i_lag=1:length(lags_all)
    x_shift = circshift(x, lags_all_idx(i_lag), ndims(x));
    index = repmat({':'}, 1, length(shape)); 
    index{end} = i_lag; 
    r(index{:}) = correlation(x, x_shift); 
end

% calculate mean acf value for lags of interest
index = repmat({':'}, 1, length(shape)); 

index{end} = [1 : length(lags_meter_rel)]; 
acf_mean_meter_rel = mean(r(index{:}), ndims(x)); 

index{end} = [length(lags_meter_rel)+1 : length(lags_all)]; 
acf_mean_meter_unrel = mean(r(index{:}), ndims(x)); 

% z-score
z = zscore(r, [], ndims(x)); 
index{end} = [1 : length(lags_meter_rel)]; 
z_meter_rel = mean(z(index{:}), ndims(x)); 

% acf ratio 
ratio_meter_rel = acf_mean_meter_rel ./ acf_mean_meter_unrel;

% acf contrast 
contrast_meter_rel = (acf_mean_meter_rel-acf_mean_meter_unrel) ./ ...
                     (acf_mean_meter_rel+acf_mean_meter_unrel); 

feat = []; 
feat.mean_meter_rel = acf_mean_meter_rel; 
feat.z_meter_rel = z_meter_rel; 
feat.ratio_meter_rel = ratio_meter_rel; 
feat.contrast_meter_rel = contrast_meter_rel; 

