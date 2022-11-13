function theta_opt = fit_aperiodic(x, y, varargin)
% Fit aperiodic 1/f component to data. 
% 
% 
% Parameters
% ----------
% x : array_like
%     1-D array of frequencies (linear!). 
% y : array_like
%     1-D array of log-power. 
% fit_knee : bool, default=false
%     Whether to use the knee parameter when fitting 1/f. 
% method : {'fminunc', 'lsq'}, default='fminunc'
%     Which matlab fitting method implementation to use. 
% init : array_like
%     [offset, knee, exponent] vector of initialization parameters if fitting
%     with knee, else [offset, exponent]
%
% Returns
% -------
% theta_opt : array_like
%     Array of 1/f parametersestimated by optimizing least-squares cost. 
%     
% 

parser = inputParser; 

addParameter(parser, 'fit_knee', false)
addParameter(parser, 'method', 'fminunc')
addParameter(parser, 'init', [])
addParameter(parser, 'robust', false)

parse(parser, varargin{:})

fit_knee = parser.Results.fit_knee; 
method = parser.Results.method; 
theta_init = parser.Results.init; 
robust = parser.Results.robust; 
          
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
                           'MaxFunctionEvaluations', 5000, ...
                           'Display', 'off'); 
    % first fit 
    if fit_knee
        lsq_fun = @(theta) ...
            (theta(1) - log10(theta(2) + x .^ theta(3))) - y;
    else
        lsq_fun = @(theta) ...
            (theta(1) - log10( x .^ theta(end))) - y;
    end
    theta_opt = lsqnonlin(lsq_fun, theta_init, [], [], options); 
    
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
        theta_opt = lsqnonlin(lsq_fun, theta_opt, [], [], options); 
    end
    
elseif strcmp(method, 'fminunc')
    
    options = optimoptions('fminunc', ...
                           'MaxFunctionEvaluations', 5000, ...
                           'Display', 'off'); 
    try
        theta_opt = fminunc(@(theta) sum((y - aperiodic(theta, x)).^2), ...
                            theta_init, options); 
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
            theta_opt = fminunc(fun, theta_opt, options); 
        end
        
    catch ME
        
        if strcmp(ME.identifier, 'optim:fminusub:UsrObjUndefAtX0')
            warning('Seems like objective function is undefined at initial point. This can be because you are trying to fit clean signal without any noise?')
            theta_opt = zeros(size(theta_init)); 
        end
        
    end
    
else
    error('method "%s" not implemented', method)
end










