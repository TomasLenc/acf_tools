%% Separate the spectra of fractal component and oscillatory component from mixed time series
%   amri_sig_fractal()
%
% Usage
%   spec = amri_sig_fractal(sig,srate,...)
% 
% Inputs
%   sig   - a time-series vector. If sig is a matrix, then separate spectra for each column  
%   srate - sampling rate
%
% Outputs
%   spec  - spectrum
%           .freq = a vector of frequency points
%           .srate= sample rate;
%           .mixd = spectrum of mixed time series
%           .frac = spectrum of fractal component
%           .osci = spectrum of oscillatory component
%
% Keywords
%   frange  - [fmin fmax](default [0 srate/4]), the output frequency range.
%   detrend - 1 or 0 (default 1): 1 means detrending data before fft, otherwise 0
%   filter  - 1 or 0 (default 1): 1 means filtering before downsampling to avoid aliasing.
%   hset    - (default 1.1:0.05:1.9) an array containing scaling factors (> 1).
%
% See also
%   amri_sig_genfrac
%   amri_sig_plawfit
%
% Version
%   0.10
% 
% Reference
%   -Wen H. and Liu Z. (2015) Separating Fractal and Oscillatory Components in 
%    the Power Spectrum of Neurophysiological Signals

%% History
% 0.01 - HGWEN - 12/15/2013 - Use resample function instead of interp1
%                           - upsample signal before extracting fractals
% 0.02 - HGWEN - 01/15/2014 - add a new method 'nivcgsa'
%                           - fit power-law line after resampling data in euqal space
%                           - set nfft=2^nextpow2(2*Ndata);
% 0.03 - HGWEN - 02/27/2014 - Change the name 'nivcgsa' to "IRASA". 
% 0.04 - HGWEN - 03/01/2014 - In IRASA, use median instead of min operator in the final step
% 0.05 - HGWEN - 05/14/2014 - If sig is a matrix, separate spectra for each column
% 0.06 - HGWEN - 08/20/2014 - Use multiple h values in CGSA
%                           - remove the power-law fitting section, and add a new function "amri_sig_plawfit"
%                           - Only return freq, srate, mixd, and frac.
% 0.07 - HGWEN - 10/11/2014 - Add a keyword "upsample"
% 0.08 - HGWEN - 10/19/2014 - Reorganized the code
% 0.09 - HCWEN - 04/11/2015 - Added keywords 'hset' and 'filter', and removed the keyword 'upsample'
% 0.10 - HGWEN - 09/26/2015 - Reorganized the structure.

%%

function spec = amri_sig_fractal(sig,srate,varargin)

if nargin<2
    eval('help amri_sig_fractal');
    return
end

%% defaults
flag_detrend = 1;
fmin = 0; 
fmax = srate/4; 
flag_filter = 1;
hset=1.1:0.05:1.9;

%% Keywords
for i = 1:2:size(varargin,2) 
    Keyword = varargin{i};
    Value   = varargin{i+1};
    if strcmpi(Keyword,'frange')
        fmin = max(Value(1),fmin);
        fmax = min(Value(2),fmax);
    elseif strcmpi(Keyword,'detrend')
        flag_detrend = Value;
    elseif strcmpi(Keyword,'filter')
        flag_filter = Value;
    elseif strcmpi(Keyword,'hset')
        hset = Value;
    else
        warning(['amri_sig_fractal(): unknown keyword ' Keyword]);
    end
end

%% preprocessing
sig = double(sig);
if isvector(sig)
   sig = sig(:); 
end

% detrend signal
if flag_detrend >= 1
    sig = detrend(sig,'linear');
end

%% apply IRASA method to separate fractal and oscillatory components
[Smixd, Sfrac, freq] = irasa(sig,srate,hset,flag_filter);  

%% only keep the given frequency range
ff = (freq>=fmin & freq<=fmax & freq>0);
freq = freq(ff); 
Smixd = Smixd(ff,:); 
Sfrac = Sfrac(ff,:); 

%% outputs
spec.freq  = freq;
spec.srate = srate;
spec.mixd  = Smixd;
spec.frac  = Sfrac;
spec.osci  = Smixd - Sfrac;

end

%% IRASA Irregular-Resampling Auto-Spectral Analysis

function [Smixd, Sfrac, freq] = irasa(sig,srate,hset,flag_filter)
% Given a discrete time series (sig) of length (Ntotal)
Ntotal = size(sig,1);
dim = size(sig,2);

% Ndata is the power of 2 that does not exceed 90% of Ntotal.
Ndata = 2^floor(log2(Ntotal*0.9));

% Nsubset is fixed to 15
Nsubset = 15;

% compute the auto-power spectrum of the originally sampled time series
L = floor((Ntotal-Ndata)/(Nsubset-1));

% set nfft greater than ceil(hset(end))*Ndata, asure that do fft without truncating
nfft = 2^nextpow2(ceil(hset(end))*Ndata);

% set output data length Nfrac
Nfrac = nfft/2 + 1;
freq = srate/2*linspace(0,1,Nfrac); freq = freq(:);

% compute the spectrum of mixed data
Smixd = zeros(Nfrac,dim);
taper = gettaper([Ndata dim]);
for k = 0:Nsubset-1
    i0 = L*k+1;
    x1 = sig(i0:1:i0+Ndata-1,:);
    p1 = fft(x1.*taper,nfft)/min(nfft,size(x1,1));
    p1(2:end,:) = p1(2:end,:)*2;
    Smixd = Smixd+abs(p1(1:Nfrac,:)).^2;
end
Smixd = Smixd/Nsubset;

% filter the input signal to avoid alising when downsampling
if flag_filter == 1
    sig_filtered = sig;
    for i = 1 : size(sig,2)
        sig_filtered(:,i) = amri_sig_filtfft(sig(:,i),srate,0,srate/(2*ceil(hset(end))));
    end
end

% compute fractal component.
Sfrac = zeros(Nfrac,dim,length(hset));
for ih = 1:length(hset)
    % compute the auto-power spectrum of xh
    h = hset(ih);
    [n, d] = rat(h); % n > d
    Sh = zeros(Nfrac,dim);
    for k = 0 : Nsubset-1
        i0 = L*k + 1;
        x1 = sig(i0:i0+Ndata-1,:);
        xh = myresample(x1, n, d); 
        taperh = gettaper(size(xh));
        ph = fft(xh.*taperh,nfft)/min(nfft,size(xh,1));
        ph(2:end,:) = ph(2:end,:)*2;
        tmp = (abs(ph)).^2;
        Sh = Sh + tmp(1:Nfrac,:);
    end
    Sh = Sh / Nsubset;
    
    % compute the auto-power spectrum of X1h
    S1h = zeros(Nfrac, dim);
    for k = 0 : Nsubset - 1
        i0 = L*k + 1;
        if (flag_filter==1)
            x1 = sig_filtered(i0:1:i0+Ndata-1,:);
        else
            x1 = sig(i0:1:i0+Ndata-1,:);
        end
        x1h = myresample(x1,d,n);
        taper1h = gettaper(size(x1h));
        p1h = fft(x1h.*taper1h,nfft)/min(nfft,size(x1h,1));
        p1h(2:end,:) = p1h(2:end,:)*2;
        tmp = (abs(p1h)).^2;
        S1h = S1h + tmp(1:Nfrac,:);
    end
    S1h = S1h / Nsubset;
    Sfrac(:,:,ih)= sqrt(Sh.*S1h);
end

% taking median
Sfrac = median(Sfrac,3);
end


%% subfunctions
function taper = gettaper(S)
% get a tapering function for power spectrum density calculation
if license('test','signal_toolbox')
    taper = hann(S(1),'periodic');    
else
    taper = 0.5*(1-cos(2*pi*(1:S(1))/(S(1)-1)));   
end
    taper = taper(:);
    taper = repmat(taper,1,S(2));
end

function data_out = myresample(data,L,D)
% resample signal with upsample L and downsample D
if license('test','signal_toolbox')
    data_out = resample(data,L,D);
else
    N = size(data,1);
    x0 = linspace(0,1,N);
    x1 = linspace(0,1,round(N*L/D));
    data_out = interp1(x0,data,x1);    
end
end
