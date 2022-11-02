function [acf, lags, ap_linear, mX, freq, ap_par] = get_acf(x, fs, varargin)
% Get autocorrelation of time-domain x. Optionally removes 1/f component. 
% 
% Parameters
% ----------
% x : array_like, shape=[..., time]
%     Input x with time as the last dimension
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
% 

parser = inputParser; 
addParameter(parser, 'rm_ap', false)
addParameter(parser, 'fit_knee', false)
addParameter(parser, 'min_freq', 0.1)
addParameter(parser, 'max_freq', fs/2)
addParameter(parser, 'f0_to_ignore', [])
addParameter(parser, 'bins', [2, 5])

parse(parser, varargin{:})

rm_ap               = parser.Results.rm_ap; 
fit_knee            = parser.Results.fit_knee; 
min_freq            = parser.Results.min_freq; 
max_freq            = parser.Results.max_freq; 
f0_to_ignore        = parser.Results.f0_to_ignore; 
bins                = parser.Results.bins; 


if isrow(f0_to_ignore)
    f0_to_ignore = f0_to_ignore'; 
end

% check if x is strictly positive
if any(x < 0)
    warning(...
        strcat('Time-domain input signal is not strictly positive. ', ...
               'This may lead to biased measures of periodicity from autocorelation.')...
           ); 
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

N = size(x, ndims(x)); 
hN = floor( (N + 1) / 2); 
lags = [0 : hN-1] / fs; 


% whole-trial FFT
N = size(x, ndims(x)); 
hN = ceil(N / 2); 
freq = [0 : hN-1]' / N * fs; 

mX = abs(fft(x, [], ndims(x))) / N * 2; 

index = cell(1, ndims(x));
index(:) = {':'};
index{end} = [1:hN];
mX = mX(index{:}); 


if ~rm_ap
    
    % get autocorrelation from FFT
    X = fft(x, [], ndims(x)) / N * 2; 
    acf = (real(ifft(X .* conj(X), [], ndims(x)))); 

else
    
    min_freq_idx = dsearchn(freq, min_freq); 
    max_freq_idx = dsearchn(freq, max_freq); 
    freq_to_fit = freq(min_freq_idx : max_freq_idx); 
    
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
        
    % fit aperiodic 
    % -------------
    
    shape = size(x);
    shape(end) = length(freq); 
    ap_all = zeros(shape);
    ap_linear = zeros(shape);
    ap_par = cell([shape(1:end-1), 1]); 

    % time is on the last dimension - we will loop over everything else
    
    % prepare index vector as cell, e.g. {1,2,3,4,':'} 
    nv = ndims(x) - 1;  % Exclude last dimension
    idx_while_loop = [repmat({1}, 1, nv), {':'}]; 
    max_size_per_dim = size(x); % size of each dimension
    ready = false; 
    
    while ~ready
      
        log_pow = log10(squeeze(mX_to_fit(idx_while_loop{:})) .^ 2); 

        log_pow_to_fit = log_pow(min_freq_idx : max_freq_idx); 

        if ~iscolumn(log_pow_to_fit)
            log_pow_to_fit = log_pow_to_fit'; 
        end
        
        % fit aperiodic component
        ap_par{idx_while_loop{:}} = fit_aperiodic(...
                            freq_to_fit, log_pow_to_fit, ...
                            'fit_knee', fit_knee); 

        % generate aperiodic with the estimated parameters
        ap = aperiodic(ap_par{idx_while_loop{:}}, freq); 
        ap_all(idx_while_loop{:}) = ap; 

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

    % use the estimated aperiodic to get normalized ACF 
    X = fft(x, [], ndims(x)) ./ N .* 2; 

    if mod(N, 2) == 0
        ap_whole_spect = cat(...
                             ndims(x),...
                             ap_linear, ...
                             flip(ap_linear, ndims(x))...
                             ); 
    else
        index = cell(1, ndims(x));
        index(:) = {':'};
        index{end} = [2 : size(ap_linear, ndims(ap_linear))];      
        ap_whole_spect = cat(...
                             ndims(x), ...
                             ap_linear, ...
                             flip(ap_linear(index{:}), ndims(x))...
                             ); 
    end

    X_norm = X ./ ap_whole_spect; 
    
    acf = real(ifft(X_norm .* conj(X_norm), [], ndims(x))); 

end


index = cell(1, ndims(x));
index(:) = {':'};
index{end} = 1:hN;   
acf = acf(index{:}); 

% set acf at lag 0 to the value of lag 1
if lags(1) == 0
    index_0 = cell(1, ndims(x));
    index_0(:) = {':'};
    index_0{end} = 1;   
    index_1 = cell(1, ndims(x));
    index_1(:) = {':'};
    index_1{end} = 2;       
    acf(index_0{:}) = acf(index_1{:}); 
end


if ~isrow(freq)
    freq = freq'; 
end

% set DC to 0
index = cell(1, ndims(x));
index(:) = {':'};
index{end} = 1;
mX(index{:}) = 0; 










