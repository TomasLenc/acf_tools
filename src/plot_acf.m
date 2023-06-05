function plot_acf(ax, acf, lags, varargin)
% Generate plot of the ACF. 
% 
% Parameters
% ----------
% ax : matlab axes object
%     The plot will be placed into these axes. 
% acf : array_like
%     Time lags must be last dimension. 
% lags : array_like
%     1-D array of time lags in seconds. 
% lags_meter_rel : array_like, optional
%     1-D array of time lags that are meter related. These will be highlighted
%     in the plot using vertical lines. 
% lags_meter_unrel : array_like, optional
%     1-D array of time lags that are meter unrelated. These will be highlighted
%     in the plot using vertical lines. 
% min_lag : float, optional
%     The smallest lag that will be shown. Default is the smallest lag
%     available in lags array. 
% max_lag : float, optional
%     The largest lag that will be shown. Default is the highest lag
%     available in lags array. 
% prec : int, optional, default=1e2
%     Rounding precision to obtain axis limits and for printing text. The
%     higher the number, the higher the precision (more decimal points are
%     considered). 
% feat : struct, optional 
%     Structure with calculated ACF features. If provided, these will be
%     printed as plot title. 
% col_acf : rgb triplet
%     Color of the plotted ACF. 
% col_meter_rel : rgb triplet
%     Color of the plotted vertical lines marking meter-related lags. 
% col_meter_unrel : rgb triplet
%     Color of the plotted vertical lines marking meter-unrelated lags. 
% linew : float
%     Linewidth of the plotted ACF. 
% linew_lagz : float
%     Linewidth of the vertical lines marking meter-unrelated lags. 
% opacity : float between 0 and 1
%     Alpha value of the plotted ACF. 
% opacity_lagz : float between 0 and 1
%     Alpha value of the vertical lines marking meter-unrelated lags. 
%
% Returns
% -------
% feat : struct
%     Structure with calculated ACF features. 
%     
    
parser = inputParser; 
addRequired(parser, 'ax'); 
addRequired(parser, 'acf', @isnumeric); 
addRequired(parser, 'lags', @(x) isnumeric(x)); 
addParameter(parser, 'min_lag', min(lags), @isnumeric); 
addParameter(parser, 'max_lag', max(lags), @isnumeric); 
addParameter(parser, 'lags_meter_rel', [], @isnumeric); 
addParameter(parser, 'lags_meter_unrel', [], @isnumeric); 
addParameter(parser, 'features', [], @isstruct); 
addParameter(parser, 'prec', 100, @isnumeric); 

addParameter(parser, 'col_acf', [0, 0, 0]); 
addParameter(parser, 'col_meter_rel', [0.8706    0.1765    0.1490]); 
addParameter(parser, 'col_meter_unrel', [0.1922    0.5098    0.7412]); 
addParameter(parser, 'linew', 2, @isnumeric); 
addParameter(parser, 'linew_lagz', 2, @isnumeric); 
addParameter(parser, 'opacity', 1, @isnumeric); 
addParameter(parser, 'opacity_lagz', 1, @isnumeric); 

parse(parser, ax, acf, lags, varargin{:});

min_lag                 = parser.Results.min_lag; 
max_lag                 = parser.Results.max_lag; 
lags_meter_rel          = parser.Results.lags_meter_rel; 
lags_meter_unrel        = parser.Results.lags_meter_unrel; 
features                = parser.Results.features; 
prec                    = parser.Results.prec; 

col_acf         = parser.Results.col_acf; 
col_meter_rel   = parser.Results.col_meter_rel; 
col_meter_unrel = parser.Results.col_meter_unrel; 
linew           = parser.Results.linew; 
linew_lagz      = parser.Results.linew_lagz; 

opacity = parser.Results.opacity; 
opacity_lagz = parser.Results.opacity_lagz; 



%%
if isrow(lags)
    min_lag_idx = dsearchn(lags', min_lag);
    max_lag_idx = dsearchn(lags', max_lag);
else
    min_lag_idx = dsearchn(lags, min_lag);
    max_lag_idx = dsearchn(lags, max_lag);
end
y_lims = [min(floor(acf(min_lag_idx:max_lag_idx)*prec)/prec), ...
          max(ceil(acf(min_lag_idx:max_lag_idx)*prec)/prec)]; 
      
hold(ax, 'on');

if ~isempty(lags_meter_rel)
    lags_meter_rel_idx = dsearchn(lags', lags_meter_rel'); 
    h = plot(ax, [lags(lags_meter_rel_idx); lags(lags_meter_rel_idx)], y_lims,...
            '-', 'linew', linew_lagz, 'color', col_meter_rel); 
    for i=1:length(h)
        h(i).Color(4) = opacity_lagz; 
    end
end

if ~isempty(lags_meter_unrel)
    lags_meter_unrel_idx = dsearchn(lags', lags_meter_unrel'); 
    h = plot(ax, [lags(lags_meter_unrel_idx); lags(lags_meter_unrel_idx)], y_lims,...
            '-', 'linew', linew_lagz, 'color', col_meter_unrel);
    for i=1:length(h)
        h(i).Color(4) = opacity_lagz; 
    end
end

h = plot(ax, lags, acf, 'marker', 'none', 'linew', linew, 'color', col_acf); 
h.Color(4) = opacity; 

x_lims = [min_lag, max_lag]; 
ax.XTick = x_lims; 
ax.XLim = x_lims; 
ax.TickDir = 'out'; 
if y_lims(1) < y_lims(2)
    ax.YTick = y_lims; 
    ax.YLim = [floor(y_lims(1)*prec)/prec, ...
               ceil(y_lims(2)*prec)/prec]; 
end

if ~isempty(features)
    tit = ''; 
    keys = fieldnames(features); 
    for i_key=1:length(keys)
        tit = [tit, sprintf('%s=%.2f   ', keys{i_key}, features.(keys{i_key}))];    
    end
    h_tit = title(ax, tit, 'Interpreter', 'none'); 
    h_tit.HorizontalAlignment = 'left'; 
    h_tit.Units = 'normalized'; 
    h_tit.Position(1) = 0; 
    h_tit.Position(2) = 1.1; 

end












