## acf_tools

**acf_tools** is a collection of MATLAB functions that can be used to work with autocorrelation in the context of rhythmic signal analysis. 

The library provides functions to compute autocorrelation function from a time-domain signal, while minimizing the bias from a 1/f noise that might be mixed in with the signal of interest. 

Additional tools allow using the resulting autocorrelation function to quantify periodic recurrence of a signal at particular rate. This can be used, for example, to estimate how prominent a particular periodic beat is represented within a signal. 

![summary diagram of the method](media/summary_fig.png)

## installation 

Clone the project directory from github and add the directory (with all subdirectories) on MATLAB path. 

_**Dependencies**_: You also need to have the [rnb_tools package](https://github.com/TomasLenc/rnb_tools) on your MATLAB path. 


## quickstart



```matlab

clear all
close all

% add the necessary library folders to MATLAB path (make 
% sure you change these to your local paths)
addpath(genpath('/datadisk/projects_git_dl/rnb_tools'));
addpath(genpath('/datadisk/projects_git_dl/acf_tools'));

% define sampling rate
sampling_rate = 200; 

% define the shortest, unitary inter-onset interval for the rhythm
grid_ioi = 0.200; 

% define structure of the rhythmic pattern (1 = event, 0 = nothing)
pat = [1 1 1 0 1 1 1 0 1 1 0 0]; 

% create impulse-response kernel 
ir = get_square_kernel(sampling_rate, ...
        'duration', 0.100, ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 

% simulate a "signal": repeating rhythmic pattern 
[x_clean, t] = get_s(pat, grid_ioi, sampling_rate, ...
					  'n_cycles', 16, 'ir', ir); 

% simulate 1/f noise
noise_exponent = -1.5; 
noise = get_colored_noise2([1, length(x_clean)], sampling_rate, noise_exponent); 

% scale the noise to the requested SNR 
snr = 1; 
x = add_signal_noise(x_clean, noise, snr);

% ==============================================================

% when fitting the 1/f noise, we will ignore peaks at frequencies 
% corresponding to the pattern repetition frequency and 
% harmonics. 
pat_f0_to_ignore = length(pat) * grid_ioi; 

% compute ACF after accounting for 1/f noise  
[acf_norm, lags, aperiodic_estimate, mX_raw, freqs] = get_acf(x, ...
								 sampling_rate, ...
                                 'rm_ap', true, ...
                                 'f0_to_ignore', pat_f0_to_ignore, ...
                                 'min_freq', 0.1, ...
                                 'max_freq', 9);
                             
% define the closest and furthest neighbouring bin (on each side) that is going
% to be used to subtract the 1/f noise from magnitude spectrum
noise_bins = [3, 13];
% normalize FFT magnitude spectrum by subtracting mean magnitude at 
% neighouring bins 
mX_norm = subtract_noise_bins(mX_raw, noise_bins(1),  noise_bins(2)); 
                             
% compute ACF WITHOUT accounting for 1/f noise  
[acf_raw, lags] = get_acf(x, sampling_rate, 'rm_ap', false);


% ==============================================================

% Let's measure the prominence of periodicity at 0.8 s. In order to normalize
% the meaure, we will contrast to periodicities at 0.6 and 1.0 s. 

% we will use all valid lags, i.e. up to half signal duration 
max_lag = length(x) / sampling_rate / 2; 

% get beat-related lags 
% ------------------

% Make sure there's no overlap with muiltiples of beat-unrelated lags
lags_beat_related = get_lag_harmonics(...
                            0.8, ...
                            max_lag,...
                            'lag_harm_to_exclude', [0.6] ...
                            ); 

% beat-unrelated lags 
% --------------------

% Make sure there's no overlap with muiltiples of beat-related lags 
lags_beat_unrel_1 = get_lag_harmonics(...
                            0.6, ...
                            max_lag,...
                            'lag_harm_to_exclude', [0.8] ...
                            ); 

lags_beat_unrel_2 = get_lag_harmonics(...
                           1.0, ...
                            max_lag,...
                            'lag_harm_to_exclude', [0.8] ...
                            ); 

lags_beat_unrelated = uniquetol([lags_beat_unrel_1, lags_beat_unrel_2], 1e-8); 

% make sure one more time that there's no overlap between beat-rel and -unrel !!! 
assert(~any( min(abs(bsxfun(@minus, lags_beat_related', lags_beat_unrelated))) < 1e-9 ))


% compute features from the autocorrelation function 
features_acf = get_acf_features(acf_norm, ...
                                lags, ...
                                lags_beat_related, ...
                                lags_beat_unrelated);

fprintf('\nThe mean z-score at beat-related lags = %.2f\n\n', features_acf.z_meter_rel);


% ==============================================================
% summary figure 

f = figure('color', 'white');

fontsize = 9; 

ax = subplot(3, 3, 1:3);
plot(ax, t, x); 
ax.XLabel.String = 'time (s)'; 
ax.YLabel.String = 'amplitude'; 
ax.FontSize = fontsize; 
box(ax, 'off'); 

% plot raw ACF
ax = subplot(3, 3, 4:5); 
plot_acf(ax, acf_raw, lags, 'prec', 1e6, 'linew', 0.7, ...
         'lags_meter_rel', lags_beat_related, ...
         'lags_meter_unrel', lags_beat_unrelated, ...
		  'opacity_lagz', 0.4); 
ax.XLabel.String = 'lag (s)'; 
ax.YLabel.String = 'autocorrelation'; 
ax.FontSize = fontsize; 
box(ax, 'off'); 

% plot raw FFT
ax = subplot(3, 3, 6); 
max_freq_to_plot = 10; 
signal_frequencies = [1/pat_f0_to_ignore : 1/pat_f0_to_ignore : max_freq_to_plot]; 
plot_fft(freqs, mX_raw, ...
         'ax', ax, ....
         'maxfreqlim', max_freq_to_plot, ...
         'frex_meter_rel', signal_frequencies, ...
         'fontsize', fontsize); 
ax.YTick = []; 
ax.XAxis.Visible = 'on'; 
ax.XLabel.String = 'frequency (Hz)'; 
ax.YLabel.String = 'magnitude'; 

% plot 1/f component
hold(ax, 'on');
plot(ax, freqs, aperiodic_estimate, '--', 'color', 'k', 'linew', 2);    

% plot noise-subtracted ACF
ax = subplot(3, 3, 7:8); 
plot_acf(ax, acf_norm, lags, 'prec', 1e6, 'linew', 0.7, ...
         'lags_meter_rel', lags_beat_related, ...
         'lags_meter_unrel', lags_beat_unrelated, ...
		  'opacity_lagz', 0.4); 
ax.XLabel.String = 'lag (s)'; 
ax.YLabel.String = 'autocorrelation'; 
ax.FontSize = fontsize; 
box(ax, 'off'); 

% plot noise-subtracted FFT
ax = subplot(3, 3, 9); 
plot_fft(freqs, mX_norm, 'ax', ax, 'maxfreqlim', 10, 'fontsize', fontsize, ...
         'frex_meter_rel', signal_frequencies); 
ax.YTick = []; 
ax.XAxis.Visible = 'on'; 
ax.XLabel.String = 'frequency (Hz)'; 
ax.YLabel.String = 'magnitude'; 


```



















































