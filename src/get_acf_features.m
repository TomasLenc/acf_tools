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
% lags_meter_unrel_left : array_like, optional 
%     1-D array of time lags that are meter unrelated and on the "left" (i.e. 
%     smaller) than meter-related lags. 
% lags_meter_unrel_right : array_like, optional 
%     1-D array of time lags that are meter unrelated and on the "right" (i.e. 
%     larger) than meter-related lags. 
% normalize_acf : bool, default=true
%     Whether to normalize the acf between 0 and 1 before calculating features.      
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
addParameter(parser, 'lags_meter_unrel_left', [])
addParameter(parser, 'lags_meter_unrel_right', [])

parse(parser, varargin{:})

normalize_acf_vals = parser.Results.normalize_acf_vals; 
lags_meter_unrel_left = parser.Results.lags_meter_unrel_left; 
lags_meter_unrel_right = parser.Results.lags_meter_unrel_right; 

% allocate
z_meter_rel = []; 
ratio_meter_rel = []; 
ratio_meter_rel_left = []; 
ratio_meter_rel_right = []; 
contrast_meter_rel = []; 

% get indices for lags of interest
lags_meter_rel_idx = ensure_row(find_idx_tol(lags, lags_meter_rel)); 
lags_meter_unrel_idx = ensure_row(find_idx_tol(lags, lags_meter_unrel)); 

if ~isempty(lags_meter_unrel_left)
    lags_meter_unrel_left_idx = ensure_row(...
        find_idx_tol(lags, lags_meter_unrel_left)); 
end
if ~isempty(lags_meter_unrel_right)
    lags_meter_unrel_right_idx = ensure_row(...
        find_idx_tol(lags, lags_meter_unrel_right)); 
end

% normalize acf between 0 and 1
if normalize_acf_vals
    warning('extracting features: values first normalized between 0 and 1'); 

%     acf_min = min(acf, [], ndims(acf)); 
%     acf_max = max(acf, [], ndims(acf)); 
%     acf_range = acf_max - acf_min; 
%     acf = (acf - acf_min) ./ acf_range; 
    
    % first get all lags of interest and normalize by the lowest & highest one
    index = repmat({':'}, 1, ndims(acf)); 
    index{end} = [lags_meter_rel_idx, lags_meter_unrel_idx]; 
    acf_vals = acf(index{:}); 
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

if ~isempty(lags_meter_unrel_left)
    subs_cmd.subs{end} = lags_meter_unrel_left_idx;
    acf_mean_meter_unrel_left = mean(subsref(acf, subs_cmd), ndims(acf)); 
end
if ~isempty(lags_meter_unrel_right)
    subs_cmd.subs{end} = lags_meter_unrel_right_idx;
    acf_mean_meter_unrel_right = mean(subsref(acf, subs_cmd), ndims(acf)); 
end

% z-score
subs_cmd.subs{end} = [lags_meter_rel_idx, lags_meter_unrel_idx];
z = zscore(subsref(acf, subs_cmd), [], ndims(acf)); 

subs_cmd.subs{end} = [1 : length(lags_meter_rel_idx)];
z_meter_rel = mean(subsref(z, subs_cmd), ndims(acf)); 

% acf ratio 
ratio_meter_rel = acf_mean_meter_rel ./ acf_mean_meter_unrel;

if ~isempty(lags_meter_unrel_left)
    ratio_meter_rel_left = acf_mean_meter_rel ./ acf_mean_meter_unrel_left;
end
if ~isempty(lags_meter_unrel_right)
    ratio_meter_rel_right = acf_mean_meter_rel ./ acf_mean_meter_unrel_right;
end

% acf contrast 
contrast_meter_rel = (acf_mean_meter_rel-acf_mean_meter_unrel) ./ ...
                     (acf_mean_meter_rel+acf_mean_meter_unrel); 

feat = []; 
feat.z_meter_rel = z_meter_rel; 
feat.ratio_meter_rel = ratio_meter_rel; 
feat.ratio_meter_rel_left = ratio_meter_rel_left; 
feat.ratio_meter_rel_right = ratio_meter_rel_right; 
feat.contrast_meter_rel = contrast_meter_rel; 



