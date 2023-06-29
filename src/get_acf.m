function [acf, lags, ap_linear, mX, freq, ap_par, x_norm, ap_optim_exitflag] = ...
                                                get_acf(x, fs, varargin)
% Get autocorrelation of time-domain x. Optionally removes 1/f component. 
% 
% Parameters
% ----------
% x : array_like, shape=[..., time]
%     Input x with time as the last dimension. 
% fs : int
%     Sampling rate. 
% rm_ap : bool, default=false
%     Whether to fit and remove the aperiodic component (1/f) from acf. 
% only_use_f0_harmonics : bool, default=true
%     If true, and the `f0_to_ignore` parameter is provided, after removing the
%     estimated 1/f component, only complex values at harmonics of f0 will be
%     retained. Everything else will be set to zero, and ACF will be computed
%     on the resulting spectrum. This option should be used when we have a
%     perfectly periodic signal, and we know the exact fundamnetal frequency at
%     which the signal will project. This is powerful information that can be
%     used to separate even better siganl from noise. 
% f0_to_ignore : float, optional
%     Fundamental frequency that will be ignored during 1/f fitting, along with 
%     all its harmonics. 
% bins : [int, int], default=[2, 5]
%     Minimum and maximum frequency bin (on both sides) used to ignore 
%     frequencies during 1/f fitting.
% ap_fit_flims : [float, float], default=[0.1, fs/2]
%     The lowest and highest frequency that will be considered when fitting 
%     the 1/f noies (default 0.1 Hz and nyquist). 
% acf_flims : [float, float], default=[0, fs]
%     The lowest and highest frequency that will be used to calculate the ACF.
%     All frequencies outside of this range will be zeroed out in the complex
%     spectrum before calculating the ACF. (default 0.1 Hz and sampling rate). 
% plot_diagnostic : bool, default=false
%     If true, a diagnostic plot for debugging In case of multidimensional 
%     input `x`, only the first element is plotted. E.g. if the input has 
%     dimensions [subject x time], the only the data for the first subject 
%     are plotted in the diagnostic plot. 
% get_x_norm : bool, default=false
%     Whether to also return 1/f-subtracted signal in the time domain. 
% fit_knee : bool, default=false
%     Whether to use the knee parameter when fitting 1/f. 
% robust : bool, optional, default=false
%     If true, the 1/f estimate will be computed using a "robust fit" procedure
%     described by Donoghue et al. 2020. 
% normalize_x : bool, default=true
%     Normalize time-domain input to mean 0 and SD 1 (i.e. zscore). 
% normalize_acf_to_1 : bool, default=false
%     Divide the resulting full ACF by its maximum value. 
% normalize_acf_z : bool, default=false
%     Normalize the restulting full ACF to mean 0 and SD 1 (i.e. zscore). 
% verbose : int {0, 1, 2}, optional, default=0
%     Verbosity level for the 1/f parameter fitting routine. The higher the 
%     number the more verbose the optimizer will be. 
% max_iter : int, optional, default=5000
%     Maximum number of iterations during 1/f parameter fitting. 
% 
% Returns 
% -------
% acf : array_like, shape=[..., lags]
%     Autocorrelation function of the input (lags are on the last dimension). 
% lags : array_like
%     Lags in seconds for the autocorrelation function. 
% ap_linear : array_like, shape=[..., frequency]
%     Fitted 1/f component for the FFT in linear scale (same shape as mX). 
% mX : array_like, shape=[..., frequency]
%     Magnitude spectra (frequencies are on the last dimension). 
% freq : array_like
%     Frequencies for the FFT. 
% ap_par : cell
%     Parameters of the fittted 1/f component. 
% x_norm : array_like
%     1/f-normalized spectrum converted back to time domain (same shape as x).
% ap_optim_exitflag : array_like
%     Optimalization exitflag for each aperiodic fit. 
% 

parser = inputParser; 
addParameter(parser, 'fit_ap', false)
addParameter(parser, 'rm_ap', false)
addParameter(parser, 'normalize_x', true)
addParameter(parser, 'force_x_positive', false)
addParameter(parser, 'normalize_acf_to_1', false)
addParameter(parser, 'normalize_acf_z', false)
addParameter(parser, 'get_x_norm', false)
addParameter(parser, 'fit_knee', false)
addParameter(parser, 'robust', false)
addParameter(parser, 'verbose', 0)
addParameter(parser, 'max_iter', 5000)
addParameter(parser, 'ap_fit_flims', [0.1, fs/2])
addParameter(parser, 'acf_flims', [0, Inf])
addParameter(parser, 'f0_to_ignore', [])
addParameter(parser, 'only_use_f0_harmonics', true)
addParameter(parser, 'bins', [2, 5])
addParameter(parser, 'plot_diagnostic', false)

parse(parser, varargin{:})

rm_ap               = parser.Results.rm_ap; 
fit_ap              = parser.Results.fit_ap | rm_ap; % if rm_ap requested we have to fit anyway
force_x_positive    = parser.Results.force_x_positive; 
normalize_x         = parser.Results.normalize_x; 
normalize_acf_to_1  = parser.Results.normalize_acf_to_1; 
normalize_acf_z     = parser.Results.normalize_acf_z; 
get_x_norm          = parser.Results.get_x_norm; 
fit_knee            = parser.Results.fit_knee; 
robust              = parser.Results.robust; 
verbose             = parser.Results.verbose; 
max_iter            = parser.Results.max_iter; 
ap_fit_flims        = parser.Results.ap_fit_flims; 
acf_flims           = parser.Results.acf_flims; 
f0_to_ignore        = parser.Results.f0_to_ignore; 
only_use_f0_harmonics = parser.Results.only_use_f0_harmonics; 
bins                = parser.Results.bins; 
plot_diagnostic     = parser.Results.plot_diagnostic; 


if isrow(f0_to_ignore)
    f0_to_ignore = f0_to_ignore'; 
end

% check if x is a row vector (common mistake, the default warnings are cryptic...)
if iscolumn(x)
    warning('input must be a row (got column instead)...transposing...')
    x = x'; 
end

% normalize input to mean 0 var 1
if normalize_x
    x = zscore(x, 1, ndims(x)); 
end

% check if x is strictly positive
if any(x(:) < 0) && force_x_positive
    warning(strcat('Forcing x to be positive...')); 
    x = x - min(x, [], ndims(x));    
end

% allocate
acf = []; 
lags = []; 
ap_linear = []; 
mX = []; 
freq = []; 
ap_linear = [];
ap_par = []; 
x_norm = []; 
freq_to_ignore = [];
freq_to_ignore_idx = []; 

% get whole-trial FFT
nyq = fs/2; 
N = size(x, ndims(x)); 
hN = floor(N / 2) + 1; 
freq = [0 : hN-1]' / N * fs; 
mX = abs(fft(x, [], ndims(x))) / N * 2; 
index = repmat({':'}, 1, ndims(x));
index{end} = [1 : hN];
mX = mX(index{:}); 

% get lags for ACF
lags = [0 : hN-1] / fs; 


% fit 1/f if requested
% --------------------

if fit_ap
    
    % check if DC is zero - warn the user 
    index = cell(1, ndims(x));
    index(:) = {':'};
    index{end} = 1;
    if any(mX(index{:}) < 1e4*eps(min(mX(:)))) && ap_fit_flims(1) <= 0 
        warning(sprintf([...
            'The magnitude at 0 Hz (DC) = 0. \n', ...
            'This will be a problem for 1/f fitting. \n', ...
            'Perhaps set the fitting range to start higher than 0?'])); 
    end
    
    % find frequencies which will be considered for the 1/f fit
    min_freq_idx = dsearchn(freq, ap_fit_flims(1)); 
    max_freq_idx = dsearchn(freq, ap_fit_flims(2)); 
    freq_to_fit = freq(min_freq_idx : max_freq_idx); 
    
    % ignore all harmonics of f0 up to nyquist frequency
    freq_to_ignore = [f0_to_ignore : f0_to_ignore : nyq]'; 
    freq_to_ignore_idx = dsearchn(freq, freq_to_ignore); 

    % for 1/f fitting, replace harmonics of f0 with mean of the bins around
    mX_to_fit = mX; 
    for i_f=1:length(freq_to_ignore)
        idx_1 = max(freq_to_ignore_idx(i_f) - bins(2), 1); 
        idx_2 = max(freq_to_ignore_idx(i_f) - bins(1), 1); 
        idx_3 = min(freq_to_ignore_idx(i_f) + bins(1), hN); 
        idx_4 = min(freq_to_ignore_idx(i_f) + bins(2), hN); 
        
        index = cell(1, ndims(x));
        index(:) = {':'};
        index{end} = [idx_1:idx_2, idx_3:idx_4];
        mean_around_bin = mean(mX(index{:}), ndims(x)); 
        
        index = cell(1, ndims(x));
        index(:) = {':'};
        index{end} = freq_to_ignore_idx(i_f);
        mX_to_fit(index{:}) = mean_around_bin; 
    end
            
    % allocate
    shape = size(x);
    shape(end) = length(freq); 
    ap_linear = zeros(shape);
    ap_par = cell([shape(1:end-1), 1]); 
    ap_optim_exitflag = nan([shape(1:end-1), 1]); 

    % time is on the last dimension - we will loop over everything else
    
    % prepare index vector as cell, e.g. {1,1,1,1,':'} 
    nv = ndims(x) - 1;  % exclude last dimension
    idx_while_loop = [repmat({1}, 1, nv), {':'}]; 
    max_size_per_dim = size(x); % size of each dimension
    ready = false; 
    
    while ~ready
      
        % get log power spectra
        log_pow = log10(squeeze(mX_to_fit(idx_while_loop{:})) .^ 2); 

        % restrict only to frequency range for 1/f fitting
        log_pow_to_fit = log_pow(min_freq_idx : max_freq_idx); 

        if ~iscolumn(log_pow_to_fit)
            log_pow_to_fit = log_pow_to_fit'; 
        end
        
        % fit aperiodic component
        [ap_par_iter, ap_optim_exitflag_iter] = fit_aperiodic(...
                            freq_to_fit, log_pow_to_fit, ...
                            'robust', robust, ...
                            'fit_knee', fit_knee,...
                            'verbose', verbose,...
                            'max_iter', max_iter...
                            ); 
                        
        ap_par{idx_while_loop{:}} = ap_par_iter;
        ap_optim_exitflag(idx_while_loop{:}) = ap_optim_exitflag_iter;

        % generate aperiodic with the estimated parameters across the whole
        % frequency range and convert from loq-power space to linear magnitude 
        % space 
        ap_log = aperiodic(ap_par{idx_while_loop{:}}, freq); 
        ap_pow = 10.^ap_log; 
        ap_linear(idx_while_loop{:}) = sqrt(ap_pow); 

        % Update the index vector:
        % Assume that the WHILE loop is ready
        ready = true;       
         % Loop through dimensions
        for k = 1:nv        
            % Increase current index by 1
            idx_while_loop{k} = idx_while_loop{k} + 1;   
            % Does it exceed the length of current dim?
            if idx_while_loop{k} <= max_size_per_dim(k) 
               % No, WHILE loop is not ready now
               ready = false;  
               % idx_while_loop(k) increased successfully, leave "for k" loop
               break;         
            end
            % Reset idx_while_loop{k}, proceed to next k
            idx_while_loop{k} = 1;          
        end
    end    

end


% get ACF
% -------

% get full complex spectra
X = fft(x, [], ndims(x)) ./ N .* 2; 

if rm_ap

    % If the first frequency bin is zero, the value of estimated aperiodic
    % value will be Inf. As we're going to divide by AP bin by bin, let's set
    % this value to 1 (instead of Inf). This will keep the DC magnitude
    % untouched. 
    ap_for_norm = ap_linear; 
    if freq(1) == 0
        index = repmat({':'}, 1, ndims(x)); 
        index{end} = 1; 
        ap_for_norm(index{:}) = 1; 
    end

    % mirror the aperiodic compoent so we also have it for negative frequencies 
    index = cell(1, ndims(x));
    index(:) = {':'};
    if mod(N, 2) == 0
        % even number of frequency bins => we have bin at pi/2 => make sure you
        % don't copy this bin twice! 
        index{end} = [2 : (size(ap_for_norm, ndims(ap_for_norm)) - 1)];      
    else
        % odd number of frequency bins => no bin at pi/2 => just skip DC and mirror
        index{end} = [2 : size(ap_for_norm, ndims(ap_for_norm))];      
    end
    ap_whole_spect = cat(...
                         ndims(x), ...
                         ap_for_norm, ...
                         flip(ap_for_norm(index{:}), ndims(x))...
                         ); 
        
    % Prepare vectors with same direction as X, and unit length (make sure we
    % don't divide by 0)
    mX_whole_spect = abs(X);
    mX_whole_spect(mX_whole_spect==0) = 1; 
    X_norm_vecs = X ./ mX_whole_spect; 
    
    % Scale them with the estimated 1/f magnitudes 
    X_ap_vecs = X_norm_vecs .* ap_whole_spect; 

    % prepare X normalized 
    X_norm = X; 

    % set X bins where the signal is BELOW noise amplitude to have amplitude of
    % zero and don't touch them again
    mask_dont_touch = abs(X) - abs(X_ap_vecs) < 0; 
    X_norm(mask_dont_touch) = 0; 
    
    % now wherever the magnitude of the observed spectrum is ABOVE the
    % estimated 1/f noise magnitude, shorten the vector at that frequency bin
    % by adding a vector with opposite direction, and magnitude equal to the
    % estimated 1/f noise magnitude at that frequency. 
    X_norm(~mask_dont_touch) = X(~mask_dont_touch) - X_ap_vecs(~mask_dont_touch); 
        
    % retain this if diagnostic plots are requested
    X_norm_all_frex = X_norm; 
    
    % If we know which frequency bins the signal is going to project to, we can
    % simply ONLY RETAIN SIGNAL FREQUENCIES and set the complex numbers at all
    % other frequency bins to zero. 
    if ~isempty(freq_to_ignore_idx) && only_use_f0_harmonics

        freq_to_keep_idx = [freq_to_ignore_idx; N - freq_to_ignore_idx + 2]; 
        index = repmat({':'}, 1, ndims(x)); 
        index{end} = freq_to_keep_idx; 

        X_norm_frex_only = zeros(size(X_norm)); 
        X_norm_frex_only(index{:}) = X_norm(index{:}); 

        X_norm = X_norm_frex_only; 
    end
        
else
    
    % without 1/f normalization
    X_norm = X;
        
end

% If requested, zero out frequencies outside of requested range before 
% computing the ACF
if acf_flims(1) == -inf
    acf_flims(1) = min(freq); 
end
if acf_flims(2) == inf
    acf_flims(2) = max(freq); 
end
flims_idx = dsearchn(ensure_col(freq), ensure_col(acf_flims)); 
flims_mask = zeros(1, hN, 'logical'); 
flims_mask(flims_idx(1) : flims_idx(2)) = 1; 
if mod(N, 2)
    % odd 
    flims_mask_all = [flims_mask, flip(flims_mask(2:end))]; 
else
    % even
    flims_mask_all = [flims_mask, flip(flims_mask(2:end-1))]; 
end
assert(length(flims_mask_all) == size(X_norm, ndims(x)))
index = repmat({':'}, 1, ndims(x)); 
index{end} = find(flims_mask_all); 

X_norm_flims = zeros(size(X_norm)); 
X_norm_flims(index{:}) = X_norm(index{:}); 
X_norm = X_norm_flims; 

% Convert 1/f-normalized spectrum back to time domain. 
if get_x_norm
    x_norm = real(ifft(X_norm, [], ndims(x))); 
end

% get the autocorrelation frunciton from the complex spectrum 
acf = real(ifft(X_norm .* conj(X_norm), [], ndims(x))); 

% scale ACF to +-1
if normalize_acf_to_1
    acf = acf ./ max(abs(acf), [], ndims(acf)); 
end

% zscore ACF
if normalize_acf_z
    acf = zscore(acf, 1, ndims(acf)); 
end

% only output acf lags up to N/2 
index = cell(1, ndims(x));
index(:) = {':'};
index{end} = 1:hN;   
acf = acf(index{:}); 

% set acf at lag 0 to the value of lag 1 (it's just variance anyway)
if lags(1) == 0
    index_0 = cell(1, ndims(x));
    index_0(:) = {':'};
    index_0{end} = 1;   
    index_1 = cell(1, ndims(x));
    index_1(:) = {':'};
    index_1{end} = 2;       
    acf(index_0{:}) = acf(index_1{:}); 
end

% set DC to 0 for output (it's just the mean)
index = cell(1, ndims(x));
index(:) = {':'};
index{end} = 1;
mX(index{:}) = 0; 

% make sure frequencies are a row vector
if ~isrow(freq)
    freq = freq'; 
end

%% diagnostic plots 

if plot_diagnostic
        
    idx_to_plot = repmat({1}, 1, ndims(x)-1); 
    idx_to_plot{ndims(x)} = ':'; 

    if ~isrow(x) && ~iscolumn(x)
        warning('Multidimensional input: diagnistic plot for x(%s)', ...
            strjoin(cellfun(@(x)num2str(x), idx_to_plot, 'uni', false), ', '));  
    end

    col_freq_lims = [48, 201, 209]/255; 

    f = figure('color', 'white', 'Position', [318 1299 1602 779]);

    pnl = panel(f);

    pnl.pack('v', [5, 95]);
    pnl(2).pack('h', [30, 70]);
    pnl(2, 1).pack('v', 2);

    % pnl.select('all');

    % ===========================================================================
    % time-domain 
    % ===========================================================================

    x_to_plot = ensure_row(squeeze(x(idx_to_plot{:}))); 

    ax = pnl(1).select();

    plot_erp(x_to_plot, 'fs', fs, 'col', [0, 0, 0], 'linew', 1, 'ax', ax);

    ax.XLim = [0, N/fs];
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'on';
    ax.TickLength = [0, 0]; 
    ax.XTick = [ax.XTick(1), ax.XTick(end)]; 
    pnl(1).xlabel('time');

    % ===========================================================================
    % magnitude spectrum for 1/f fitting
    % ===========================================================================

    % plot raw FFT
    mX_to_plot = ensure_row(squeeze(mX(idx_to_plot{:}))); 

    ax = pnl(2, 1, 1).select(); 
    title(ax, 'raw mX (up to nyquist)'); 
    hold(ax, 'on');

    plot_fft(freq, mX_to_plot, ...
             'ax', ax, ...
             'frex_meter_rel', freq_to_ignore, ...
             'maxfreqlim', nyq, ...
             'linew', 1); 
    ax.YAxis.Visible = 'on';
    ax.XAxis.Visible = 'on';
    ax.XTick = [0, nyq];
    if ~isempty(freq_to_ignore)
        ax.YLim = [0, max(mX_to_plot(freq_to_ignore_idx))];
    end
    
    
    if rm_ap
            
        % plot the 1/f estimate over the spectrum 
        ap_to_plot = ensure_row(squeeze(ap_linear(idx_to_plot{:}))); 
        plot(ax, freq, ap_to_plot, '--', 'color', 'k', 'linew', 2);    
        
        % overlay frequency limits used to fit 1/f 
        plot(ax, [ap_fit_flims(1), ap_fit_flims(1)], [0, ax.YLim(2)], ':', ...
            'linew', 3, 'color', col_freq_lims); 
        plot(ax, [ap_fit_flims(2), ap_fit_flims(2)], [0, ax.YLim(2)], ':', ...
            'linew', 3, 'color', col_freq_lims); 

        %--------
        
        % plot spectrum used to fit 1/f component
        mX_to_fit_plot = ensure_row(squeeze(mX_to_fit(idx_to_plot{:}))); 

        ax = pnl(2, 1, 2).select(); 
        title(ax, 'mX to estimate 1/f (up to nyquist)'); 
        hold(ax, 'on');

        plot_fft(freq, mX_to_fit_plot, ...
                 'ax', ax, ...
                 'frex_meter_rel', freq_to_ignore, ...
                 'maxfreqlim', nyq, ...
                 'linew', 1); 
        ax.YAxis.Visible = 'off';
        if ~isempty(freq_to_ignore)
            ax.YLim = [0, max(mX_to_plot(freq_to_ignore_idx))];
        end
        ax.XAxis.Visible = 'on';
        ax.XTick = [0, nyq];

        % plot the 1/f estimate over the spectrum 
        ap_to_plot = ensure_row(squeeze(ap_linear(idx_to_plot{:}))); 
        plot(ax, freq, ap_to_plot, '--', 'color', 'k', 'linew', 2);    
        
        % overlay frequency limits used to fit 1/f 
        plot(ax, [ap_fit_flims(1), ap_fit_flims(1)], [0, ax.YLim(2)], ':', ...
            'linew', 3, 'color', col_freq_lims); 
        plot(ax, [ap_fit_flims(2), ap_fit_flims(2)], [0, ax.YLim(2)], ':', ...
            'linew', 3, 'color', col_freq_lims); 
        
    end


    % ===========================================================================
    % full spectrum 
    % ===========================================================================

    freq_all = [0 : N-1] / N * fs;
    freq_to_keep_idx = [freq_to_ignore_idx; N - freq_to_ignore_idx + 2]; 

    pnl(2, 2).pack('v', 4);

    mX_full_to_plot = ensure_row(abs(squeeze(X(idx_to_plot{:})))); 

    ax = pnl(2, 2, 1).select(); 
    title(ax, 'mX (up to fs)'); 
    plot_fft(freq_all, mX_full_to_plot, ...
             'frex_meter_rel', freq_all(freq_to_keep_idx), ...
             'ax', ax, ...
             'linew', 0.2); 
    ax.XLim = [0, fs];
    ax.YAxis.Visible = 'off';

    if rm_ap                
        % add 1/f component to the plot
        ap_whole_to_plot = ensure_row(squeeze(ap_whole_spect(idx_to_plot{:}))); 
        
        ax.YLim = [0, max(mX_full_to_plot(freq_to_keep_idx))];
        
        hold(ax, 'on');
        plot(ax, freq_all, ap_whole_to_plot, '--', 'color', 'k', 'linew', 2);    
        plot(ax, [ap_fit_flims(1), ap_fit_flims(1)], [0, ax.YLim(2)], ':', ...
            'linew', 3, 'color', col_freq_lims); 
        plot(ax, [ap_fit_flims(2), ap_fit_flims(2)], [0, ax.YLim(2)], ':', ...
            'linew', 3, 'color', col_freq_lims); 
        
        % plot 1/f-normalized spectrum
        X_norm_to_plot = ensure_row(squeeze(X_norm_all_frex(idx_to_plot{:}))); 
        
        ax = pnl(2, 2, 2).select(); 
        title(ax, '1/f-normalized mX (up to fs)'); 
        plot_fft(freq_all, abs(X_norm_to_plot), ...
                 'frex_meter_rel', freq_all(freq_to_keep_idx), ...
                 'ax', ax, ...
                 'linew', 0.2); 
        ax.XLim = [0, fs];
        ax.YAxis.Visible = 'off';
    end

    if exist('X_norm_frex_only', 'var')
        % plot full spectrum after only harmonics of f0 are retained 
        ax = pnl(2, 2, 3).select(); 
        title(ax, '1/f-normalized mX, only f0 harmonics (up to fs)'); 
        plot_fft(freq_all, ...
             abs(ensure_row(squeeze(X_norm_frex_only(idx_to_plot{:})))), ...
             'frex_meter_rel', freq_all(freq_to_keep_idx), ...
             'ax', ax, ...
             'linew', 0.2); 
        ax.XAxis.Visible = 'on'; 
        ax.XLim = [0, fs];
        ax.XTick = [0, fs]; 
        ax.XTickLabel = ax.XTick;
        ax.TickLength = [0, 0]; 
        ax.YAxis.Visible = 'off';

    end
    
    if exist('X_norm_flims', 'var')
        % plot full spectrum after only frequencies within the requested limits
        % are retained 
        mX_norm_flims_to_plot = ...
                    abs(ensure_row(squeeze(X_norm_flims(idx_to_plot{:})))); 
        ax = pnl(2, 2, 4).select(); 
        title(ax, '1/f-normalized mX, only f0 harmonics, flims applied (up to fs)'); 
        hold(ax, 'on');
 
        clusts = bwconncomp(flims_mask_all); 
        
        for i_clust=1:length(clusts.PixelIdxList)
            
            fill(ax, ...
                [freq_all(clusts.PixelIdxList{i_clust}(1)), ...
                 freq_all(clusts.PixelIdxList{i_clust}(end)), ...
                 freq_all(clusts.PixelIdxList{i_clust}(end)), ...
                 freq_all(clusts.PixelIdxList{i_clust}(1)), ...
                 ], ...
                [min(mX_norm_flims_to_plot), min(mX_norm_flims_to_plot), ...
                 max(mX_norm_flims_to_plot), max(mX_norm_flims_to_plot)], ...
                [255, 187, 0] / 255, 'LineStyle', 'none', 'FaceAlpha', 0.1); 
        end
        
        plot_fft(freq_all, ...
             mX_norm_flims_to_plot, ...
             'frex_meter_rel', freq_all(freq_to_keep_idx), ...
             'ax', ax, ...
             'linew', 0.2); 
        ax.XAxis.Visible = 'on'; 
        ax.XLim = [0, fs];
        ax.XTick = [0, fs]; 
        ax.XTickLabel = ax.XTick;
        ax.TickLength = [0, 0]; 
        ax.YAxis.Visible = 'off';

    end
    
    pnl(2, 1).xlabel('frequency');
    pnl(2, 2).xlabel('frequency');

    pnl(2, 2).de.margin = [10, 1, 3, 12]; 
    pnl(2).margintop = 15; 
    pnl.margin = [10, 10, 3, 1]; 
    
            
end













