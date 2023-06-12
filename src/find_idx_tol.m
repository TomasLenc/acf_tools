function idx = find_idx_tol(arr, to_find, varargin)
% Find indices for closest elemets in a target array. This finction takes a
% target array `arr` and for each element of `to_find` array, it finds the
% index of the closest element in `arr`. It also checks for repetition in
% `to_find` array and issues a warning if the closest found element is too far
% away. 
%
% Parameters
% ----------
% arr : array_like 
%     Array of elements, e.g. lags, frequencies, time, etc. 
% to_find : array_like
%     Array of elements that we want to search for in `arr`. 
% tol_to_warn : float, optional, default=0.1
%     If the closest found element in `arr` is further than `tol_to_warn`, a
%     warning is issued. This helps to catch bugs, e.g. requesting to find
%     indices of frequencies above the nyqist etc. 
%     
% Returns
% -------
% idx : array_like
%     Array of indices, for each element of `to_find`, we have a corresponding
%     value indexing the closesrt element in `arr`. 

parser = inputParser; 

addParameter(parser, 'tol_to_warn', 0.1); 

parse(parser, varargin{:});

tol_to_warn = parser.Results.tol_to_warn;

% check that we got 1D array inputs
if ~ ( iscolumn(arr) || isrow(arr))
    error('find_idx_tol:badShape', ...
        '`arr` parameter must be a column or row vector!'); 
end
if ~ ( iscolumn(to_find) || isrow(to_find))
    error('find_idx_tol:badShape', ...
        '`to_find` parameter must be a column or row vector!'); 
end

% check if there are repeated values in the `to_find` array
if length(uniquetol(to_find, 1e4*eps(min(to_find)))) < length(to_find)
    warning('Repeated values in `to_find` detected...are you sure?'); 
end

idx = nan(1, length(to_find)); 

for i=1:length(to_find)
    
    [min_diff, idx(i)] = min(abs(arr - to_find(i))); 
        
    if min_diff > tol_to_warn
        warning('find_idx_tol:elementTooFar', ...
            'For %g, the closest element we found is %g. Hope you know what you are doing!', ...
                 to_find(i), arr(idx(i))); 
    end
    
end













