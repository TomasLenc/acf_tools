function [lags_meter_rel, lags_meter_unrel] = get_meter_lags(...
                    max_lag, ...
                    lag_base_incl_meter_rel, lag_base_excl_meter_rel, ...
                    lag_base_incl_meter_unrel, lag_base_excl_meter_unrel, ...
                    varargin)

parser = inputParser; 

addParameter(parser, 'min_lag', 0); 

parse(parser, varargin{:}); 

min_lag = parser.Results.min_lag; 

% meter-related lags 
% ------------------

% Make sure there's no overlap with muiltiples of meter-unrelated lags, and 
% also the pattern repetition period. 
lags_meter_rel = get_lag_harmonics(...
                            lag_base_incl_meter_rel, ...
                            max_lag,...
                            'lag_harm_to_exclude', lag_base_excl_meter_rel ...
                            ); 

% meter-unrelated lags 
% --------------------

% Make sure there's no overlap with muiltiples of meter-related lags 
% (even 0.4 seconds!), and also the pattern repetition period. 

lags_meter_unrel = get_lag_harmonics(...
                            lag_base_incl_meter_unrel, ...
                            max_lag,...
                            'lag_harm_to_exclude', lag_base_excl_meter_unrel ...
                            ); 


lags_meter_rel = lags_meter_rel(lags_meter_rel > min_lag); 
lags_meter_unrel = lags_meter_unrel(lags_meter_unrel > min_lag); 
                        
% make sure one more time that there's no overlap between meter-rel and -unrel !!! 
assert(~any( min(abs(bsxfun(@minus, lags_meter_rel', lags_meter_unrel))) < 1e-9 ))

