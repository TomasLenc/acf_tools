function lags = get_lag_harmonics(lag_bases, max_lag, varargin)
% This function returns all higher "harmonics" (i.e. multiples") of a lag, up
% to the requested max_lag. Also, it will make sure that all multiples of any
% lag pass in the array lag_harm_to_exclude will not be included. 
% 
% Parameters
% ----------
% lag_bases : array of floats
%     Lag(s) to get harmonics from.
% max_lag : float
%     Value up to which to get harmonics.
% lag_harm_to_exclude : array of floats
%     Lag values for which, all overlapping harmoincs with the requested lag will
%     be excluded.
% 
% Returns
% -------
% lags : array of floats
%     All harmonics of the lag, that don't overlap with the multiples of 
%     lag_harm_to_exlclude
%

parser = inputParser; 

addParameter(parser, 'lag_harm_to_exclude', []); 

parse(parser, varargin{:}); 

lag_bases_to_exclude = parser.Results.lag_harm_to_exclude; 


% generate all nonoveralping multiples of the requested lags
lags = []; 
for i_lag=1:length(lag_bases)
    lags = [lags, lag_bases(i_lag) : lag_bases(i_lag) : max_lag]; 
end
lags = uniquetol(lags, 1e-8); 

% remove lags overapping with multiples of lag_bases_to_exclude
mask_to_rm = zeros(1, length(lags), 'like', false); 

for i_excl=1:length(lag_bases_to_exclude)
    
    lag_to_exclude_harm = [...
        lag_bases_to_exclude(i_excl) : lag_bases_to_exclude(i_excl) : max_lag...
        ]; 
        
    for li=1:length(lags)
       
        if any(abs(lags(li) - lag_to_exclude_harm) < 1e-8)
            mask_to_rm(li) = true; 
        end
        
    end
    
end

lags = lags(~mask_to_rm); 

