function lags = get_lag_harmonics(lag, max_lag, varargin)
% This function returns all higher "harmonics" (i.e. multiples") of a lag, up
% to the requested max_lag. Also, it will make sure that all multiples of any
% lag pass in the array lag_harm_to_exclude will not be included. 
% 
% Parameters
% ----------
% lag : float
%     Lag to get harmonics from.
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

lag_harm_to_exclude = parser.Results.lag_harm_to_exclude; 


lags = [lag : lag : max_lag]; 

mask_to_rm = zeros(1, length(lags), 'like', false); 

for i_excl=1:length(lag_harm_to_exclude)
    
    lag_to_exclude_harm = [...
        lag_harm_to_exclude(i_excl) : lag_harm_to_exclude(i_excl) : max_lag...
        ]; 
        
    for li=1:length(lags)
       
        if any(abs(lags(li) - lag_to_exclude_harm) < 1e-8)
            mask_to_rm(li) = true; 
        end
        
    end
    
end

lags = lags(~mask_to_rm); 