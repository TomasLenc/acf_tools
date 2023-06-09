function feat = get_fft_features(mX, freq, freq_meter_rel, freq_meter_unrel, varargin)
% Calculate SNR as zscore across harmonics. 
% 
% Parameters
% ----------
% mX : array_like, shape=[..., frequency]
%     Raw (not noise-subtracted!!!) magnitude spectra with frequency as the 
%     last dimension. 
% freq : array_like
%     Frequencies for the FFT. 
% freq_meter_rel : array_like
%     1-D array of frequencies (in Hz) that are meter related. 
% freq_meter_unrel : array_like
%     1-D array of frequencies (in Hz) that are meter unrelated. 
% tol : float, optional, default=1e-8
%     Tolarance for finding lags of interest in the array of lags. If there is
%     no lag closer than `tol`, the requested lag will be ommitted. 
% 
% Returns 
% -------
% feat : struct
%     Structure with calculated features. 

parser = inputParser; 

addParameter(parser, 'tol', 1e-3)

parse(parser, varargin{:})

tol = parser.Results.tol; 


if ~isrow(freq)
    freq = freq';
end
if ~isrow(freq_meter_rel)
    freq_meter_rel = freq_meter_rel';
end
if ~isrow(freq_meter_unrel)
    freq_meter_unrel = freq_meter_unrel';
end

idx_meter_rel = ensure_row(find_idx_tol(freq, freq_meter_rel, 'tol', tol)); 
idx_meter_unrel = ensure_row(find_idx_tol(freq, freq_meter_unrel, 'tol', tol)); 

index = cell(1, ndims(mX));
index(:) = {':'};
index{end} = [idx_meter_rel, idx_meter_unrel];
amps = mX(index{:});

z = zscore(amps, [], ndims(mX)); 

index = cell(1, ndims(mX));
index(:) = {':'};
index{end} = [1 : length(idx_meter_rel)];
z_meter_rel = mean(z(index{:}), ndims(mX)); 

feat = []; 
feat.z_meter_rel = z_meter_rel; 
