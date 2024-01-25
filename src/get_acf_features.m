function feat = get_acf_features(acf, lags, ...
                                 lags_meter_rel, lags_meter_unrel, ...
                                 varargin)
% Extract features from acf. 
% 
% Parameters
% ----------
% acf : array_like
%     Time lags must be last dimension. 
% lags : array_like
%     1-D array of time lags in seconds
% lags_meter_rel : array_like
%     1-D array of time lags that are meter related. 
% lags_meter_unrel : array_like
%     1-D array of time lags that are meter unrelated. 
% normalize_acf : bool, default=true
%     Whether to normalize the acf between 0 and 1 before calculating features.      
% idx_find_tol : float, optional, default=0.1
%     Tolarance for finding lags of interest in the array of lags. If there is
%     no lag closer than `tol`, a warning will be issued. 
%
% Returns
% -------
% feat : struct
%     Structure with calculated features. 
%     
% Notes
% -----
% 
% 

parser = inputParser; 

addParameter(parser, 'normalize_acf_vals', false)
addParameter(parser, 'idx_find_tol', 0.1)

parse(parser, varargin{:})

normalize_acf_vals = parser.Results.normalize_acf_vals; 
idx_find_tol = parser.Results.idx_find_tol; 

% allocate
z_meter_rel = []; 
ratio_meter_rel = []; 
contrast_meter_rel = []; 

% get indices for lags of interest
lags_meter_rel_idx = ensure_row(...
    find_idx_tol(lags, lags_meter_rel, 'tol_to_warn', idx_find_tol)); 
lags_meter_unrel_idx = ensure_row(...
    find_idx_tol(lags, lags_meter_unrel, 'tol_to_warn', idx_find_tol)); 

% extract all lags of interest 
index = repmat({':'}, 1, ndims(acf)); 
index{end} = [lags_meter_rel_idx, lags_meter_unrel_idx]; 
acf_vals = acf(index{:}); 

% normalize acf between 0 and 1
if normalize_acf_vals
    warning('extracting features: values first normalized between 0 and 1'); 

%     acf_min = min(acf, [], ndims(acf)); 
%     acf_max = max(acf, [], ndims(acf)); 
%     acf_range = acf_max - acf_min; 
%     acf = (acf - acf_min) ./ acf_range; 
    
    % get all lags of interest and normalize by the lowest & highest one
    vals_min = min(acf_vals, [], ndims(acf)); 
    vals_max = max(acf_vals, [], ndims(acf)); 
    acf = (acf - vals_min) ./ (vals_max - vals_min); 
    
%     % just make the acf positive
%     norm_min = min(acf, [], ndims(acf)); 
%     acf = acf + norm_min; 
    
end

% calculate mean acf value for lags of interest
subs_cmd = []; 
subs_cmd.subs = repmat({':'}, 1, ndims(acf)); 
subs_cmd.type = '()'; 

subs_cmd.subs{end} = lags_meter_rel_idx;
acf_mean_meter_rel = mean(subsref(acf, subs_cmd), ndims(acf)); 

subs_cmd.subs{end} = lags_meter_unrel_idx;
acf_mean_meter_unrel = mean(subsref(acf, subs_cmd), ndims(acf)); 

% z-score
subs_cmd.subs{end} = [lags_meter_rel_idx, lags_meter_unrel_idx];
z = zscore(subsref(acf, subs_cmd), [], ndims(acf)); 

subs_cmd.subs{end} = [1 : length(lags_meter_rel_idx)];
z_meter_rel = mean(subsref(z, subs_cmd), ndims(acf)); 

% acf ratio 
ratio_meter_rel = acf_mean_meter_rel ./ acf_mean_meter_unrel;

% acf contrast 
contrast_meter_rel = (acf_mean_meter_rel-acf_mean_meter_unrel) ./ ...
                     (acf_mean_meter_rel+acf_mean_meter_unrel); 

feat = []; 
feat.z_meter_rel = z_meter_rel; 
feat.ratio_meter_rel = ratio_meter_rel; 
feat.contrast_meter_rel = contrast_meter_rel; 
feat.vals = acf_vals; 


