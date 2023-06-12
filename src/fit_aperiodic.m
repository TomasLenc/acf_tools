function [theta_opt, exitflag] = fit_aperiodic(x, y, varargin)
% Fit aperiodic 1/f component to data. The model is equivalent to the one
% described in: 
% T. Donoghue et al., “Parameterizing neural power spectra into periodic 
%   and aperiodic components,” Nature Neuroscience, vol. 23, no. 12,
%   pp. 1655–1665, 2020, doi: 10.1038/s41593-020-00744-x.
% 
% 
% Parameters
% ----------
% x : array_like
%     1-D array of frequencies (linear!). 
% y : array_like
%     1-D array of log-power. 
% fit_knee : bool, optional, default=false
%     Whether to use the knee parameter when fitting 1/f. 
% method : {'fminunc', 'lsq'}, optional, default='fminunc'
%     Which matlab fitting method implementation to use. 
% robust : bool, optional, default=false
%     If true, the 1/f estimate will be computed using a "robust fit" procedure
%     described by Donoghue et al. 2020.
% init : array_like, optional
%     [offset, knee, exponent] vector of initialization parameters if fitting
%     with knee, else [offset, exponent]
% verbose : int {0, 1, 2}, optional, default=0
%     Verbosity level for the optimization routine. The higher the number the
%     more verbose the optimizer will be. We will keep it silent by default. 
% max_iter : int, optional, default=5000
%     Maximum number of iterations during parameter fitting. 
%
% Returns
% -------
% theta_opt : array_like
%     Array of 1/f parametersestimated by optimizing least-squares cost. 
% exitflag: int
%     Integer code output from the optimizer, indicating whether there were any
%     issues. 
% 

parser = inputParser; 

addParameter(parser, 'fit_knee', false)
addParameter(parser, 'method', 'fminunc')
addParameter(parser, 'init', [])
addParameter(parser, 'robust', false)
addParameter(parser, 'verbose', 0)
addParameter(parser, 'max_iter', 5000)

parse(parser, varargin{:})

fit_knee = parser.Results.fit_knee; 
method = parser.Results.method; 
theta_init = parser.Results.init; 
robust = parser.Results.robust; 
verbose = parser.Results.verbose; 
max_iter = parser.Results.max_iter; 
          
% check input types
if ~isa(x, 'double') 
    error('fit_aperiodic:inputNotDouble', ...
        'x must be double float (got %s instead)', class(x)); 
end
if ~isa(y, 'double') 
    error('fit_aperiodic:inputNotDouble', ...
        'y must be double float (got %s instead)', class(y)); 
end

% check input shapes 
if ~all(size(x) == size(y))
    error('x and y must have the same shape')
end
if ndims(x) > 2 
    error('x must be a column or row vector, got %d-dim instead', ndims(x))
end
if ndims(y) > 2 
    error('y must be a column or row vector, got %d-dim instead', ndims(y))
end
if ~isrow(x) & ~iscolumn(x)
    error('x must be a column or row vector (got matrix instead)')
end
if ~isrow(y) & ~iscolumn(y)
    error('y must be a column or row vector (got matrix instead)')
end

if verbose == 0
    display_opt = 'off';
elseif verbose == 1
    display_opt = 'notify';
elseif verbose > 1
    display_opt = 'final';
end

% guess params for initialization
offset_guess = 0; 
knee_guess = 1; 
exponent_guess = 1;

if isempty(theta_init)
    if fit_knee
        theta_init = [offset_guess, knee_guess, exponent_guess]; 
    else
        theta_init = [offset_guess, exponent_guess]; 
    end
end

% These two methods (lsqnonlin and fminunc) should be equivalent 
if strcmp(method, 'lsq')
    
    options = optimoptions('lsqnonlin', ...
                           'MaxFunctionEvaluations', max_iter, ...
                           'Display', display_opt); 
    % first fit 
    if fit_knee
        lsq_fun = @(theta) ...
            (theta(1) - log10(theta(2) + x .^ theta(3))) - y;
    else
        lsq_fun = @(theta) ...
            (theta(1) - log10( x .^ theta(end))) - y;
    end
    [theta_opt, ~, exitflag, output] = lsqnonlin(...
                            lsq_fun, theta_init, [], [], options...
                            ); 
    
    % if robust fit
    if robust
        % flatten specrum based on initial fit 
        y_flat = y - aperiodic(theta_opt, x); 
        % remove negative residual values 
        y_flat(y_flat < 0) = 0; 
        % ignore 2.5% of datapoints with max error 
        perc_thr = prctile(y_flat, 0.025); 
        mask = y_flat <= perc_thr; 
        x_masked = x(mask); 
        y_masked = y(mask); 
        % fit again, using the raw first fit as initial parameters
        if fit_knee
            lsq_fun = @(theta) ...
                (theta(1) - log10(theta(2) + x_masked .^ theta(3))) - y_masked;
        else
            lsq_fun = @(theta) ...
                (theta(1) - log10( x_masked .^ theta(end))) - y_masked;
        end
        [theta_opt, ~, exitflag, output] = lsqnonlin(...
                                lsq_fun, theta_opt, [], [], options...
                                ); 
    end
    
elseif strcmp(method, 'fminunc')
    
    options = optimoptions('fminunc', ...
                           'MaxFunctionEvaluations', max_iter, ...
                           'Display', display_opt); 
    try
        [theta_opt, ~, exitflag, output] = fminunc(...
                            @(theta) sum((y - aperiodic(theta, x)).^2), ...
                            theta_init, ...
                            options...
                            ); 
        % if robust fit
        if robust
            % flatten specrum based on initial fit 
            y_flat = y - aperiodic(theta_opt, x); 
            % remove negative residual values 
            y_flat(y_flat < 0) = 0; 
            % ignore 2.5% of datapoints with max error 
            perc_thr = prctile(y_flat, (1 - 0.025) * 100); 
            mask = y_flat <= perc_thr; 
            x_masked = x(mask); 
            y_masked = y(mask); 
            % fit again, using the raw first fit as initial parameters
            fun = @(theta) sum((y_masked - aperiodic(theta, x_masked)).^2); 
            [theta_opt, ~, exitflag, output] = fminunc(fun, theta_opt, options); 
        end
        
    catch ME
        
        if strcmp(ME.identifier, 'optim:fminusub:UsrObjUndefAtX0')
            warning('Seems like objective function is undefined at initial point. This can be because you are trying to fit clean signal without any noise?')
            theta_opt = zeros(size(theta_init)); 
        else
            rethrow(ME);
        end
        exitflag = 999;
    end
    
else
    error('method "%s" not implemented', method)
end

if exitflag == 0
    warning('aperiodic fit didnt converge...')
elseif exitflag < 0
    warning('problem with optimalization method...please inspect!')
end








