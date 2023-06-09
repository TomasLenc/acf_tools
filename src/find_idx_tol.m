function idx_all = find_idx_tol(arr, to_find, varargin)

parser = inputParser; 

addParameter(parser, 'tol', 1e-9); 

parse(parser, varargin{:});

tol = parser.Results.tol;

if ~ ( iscolumn(arr) || isrow(arr))
    error('find_idx_tol:badShape', ...
        '`arr` parameter must be a column or row vector!'); 
end
if ~ ( iscolumn(to_find) || isrow(to_find))
    error('find_idx_tol:badShape', ...
        '`to_find` parameter must be a column or row vector!'); 
end


idx_all_tmp = nan(1, length(to_find)); 

if length(uniquetol(to_find, tol)) < length(to_find)
    warning('Repeated values detected...are you sure?'); 
end

c = 0; 
for i=1:length(to_find)
    [min_diff, idx] = min(abs(arr - to_find(i))); 
    if min_diff < tol
        c = c+1; 
        idx_all_tmp(c) = idx; 
    else
        warning('No value closer than %g to %g. Skipping!', ...
                 tol, to_find(i)); 
    end
end

idx_all = idx_all_tmp(1:c); 

if isempty(idx_all)
    error('find_idx_tol:NoIdxFound', 'No indices found! Aborting'); 
end
    












