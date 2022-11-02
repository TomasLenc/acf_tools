function y = aperiodic(theta, x)
% Calculate aperiodic 1/f component based on array of frequencies. 
%
% Parameters
% ----------
% theta : array_like 
%     Parameters of the 1/f (length 2: no knee, or length 3: with knee). 
% x : array_like
%     1-D array of frequencies for which the 1/f will be calculated. 
%     
% Returns
% -------
% y : array_like
%     Aperiodic component (same shape as x). 
if length(theta) == 3
    y = theta(1) - log10(theta(2) + x .^ theta(3)); 
elseif length(theta) == 2
    y = theta(1) - log10(x .^ theta(2)); 
else
    error('theta must have 2 or 3 parameters, got %d instead', length(theta))
end
