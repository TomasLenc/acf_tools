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
% fit_knee : bool, default=false
%     Whether to use the knee parameter when fitting 1/f. 
% min_freq : float, default=0.1
%     The lowest frequency that will be considered when fitting 1/f. 
% max_freq : float, default=fs/2
%     The highest frequency that will be considered when fitting 1/f (default 
%     is nyquist). 
% f0_to_ignore : float
%     Fundamental frequency that will be ignored during 1/f fitting, along with 
%     all its harmonics. 
% bins : [int, int], default=[2, 5]
%     Minimum and maximum frequency bin (on both sides) used to ignore 
%     frequencies during 1/f fitting.
% get_x_norm : bool, default=false
%     Whether to also return 1/f-subtracted signal in the time domain. 
% normalize_x : bool, default=true
%     Normalize time-domain input to mean 0 and SD 1 (i.e. zscore). 
% normalize_acf_to_1 : bool, default=false
%     Divide the resulting ACF by its maximum value. 
% normalize_acf_z : bool, default=false
%     Normalize the restulting ACF to mean 0 and SD 1 (i.e. zscore). 
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
addParameter(parser, 'min_freq', 0.1)
addParameter(parser, 'max_freq', fs/2)
addParameter(parser, 'f0_to_ignore', [])
addParameter(parser, 'bins', [2, 5])

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
min_freq            = parser.Results.min_freq; 
max_freq            = parser.Results.max_freq; 
f0_to_ignore        = parser.Results.f0_to_ignore; 
bins                = parser.Results.bins; 


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
ap_all = [];
ap_linear = [];
ap_par = []; 
x_norm = []; 

% get whole-trial FFT
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
    
    % find frequencies which will be considered for the 1/f fit
    min_freq_idx = dsearchn(freq, min_freq); 
    max_freq_idx = dsearchn(freq, max_freq); 
    freq_to_fit = freq(min_freq_idx : max_freq_idx); 
    
    % ignore all harmonics of f0 up to nyquist frequency
    nyq = fs/2; 
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
    ap_all = zeros(shape);
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
        % frequency range 
        ap = aperiodic(ap_par{idx_while_loop{:}}, freq); 
        ap_all(idx_while_loop{:}) = ap; 

        % convert from loq-power space to linear magnitude space 
        ap_pow = 10.^ap; 
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

if rm_ap

    % get full complex spectra
    X = fft(x, [], ndims(x)) ./ N .* 2; 
    
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

    % Divide each frequency bin by the estimated 1/f magnitude at that
    % frequency -> if the spectrum precisely follows the 1/f component, all
    % magnitudes shold be normalized to 0. 
    X_norm = X ./ ap_whole_spect; 
    
    % Convert 1/f-normalized spectrum back to time domain. 
    if get_x_norm
        x_norm = real(ifft(X_norm, [], ndims(x))); 
    end
    
    % get the autocorrelation frunciton from the normalized spectrum 
    acf = real(ifft(X_norm .* conj(X_norm), [], ndims(x))); 

else
    
    % get raw autocorrelation (without 1/f normalization)
    X = fft(x, [], ndims(x)) / N * 2; 
    acf = (real(ifft(X .* conj(X), [], ndims(x)))); 
        
end

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











