function plot_acf(ax, acf, lags, varargin)
     
parser = inputParser; 
addRequired(parser, 'ax'); 
addRequired(parser, 'acf', @isnumeric); 
addRequired(parser, 'lags', @(x) isnumeric(x)); 

addParameter(parser, 'ap', [], @isnumeric); 

addParameter(parser, 'min_lag', min(lags), @isnumeric); 
addParameter(parser, 'max_lag', max(lags), @isnumeric); 
addParameter(parser, 'lags_meter_rel', [], @isnumeric); 
addParameter(parser, 'lags_meter_unrel', [], @isnumeric); 
addParameter(parser, 'features', [], @isstruct); 
addParameter(parser, 'prec', 100, @isnumeric); 

addParameter(parser, 'col_acf', [0, 0, 0]); 
addParameter(parser, 'col_meter_rel', [0.8706    0.1765    0.1490]); 
addParameter(parser, 'col_meter_unrel', [0.1922    0.5098    0.7412]); 
addParameter(parser, 'col_ap', [0.8784    0.4588    0.1137]); 
addParameter(parser, 'linew', 2, @isnumeric); 
addParameter(parser, 'linew_ap', 2, @isnumeric); 
addParameter(parser, 'linew_lagz', 2, @isnumeric); 
addParameter(parser, 'opacity', 1, @isnumeric); 
addParameter(parser, 'opacity_lagz', 1, @isnumeric); 

parse(parser, ax, acf, lags, varargin{:});

ap = parser.Results.ap; 

min_lag                 = parser.Results.min_lag; 
max_lag                 = parser.Results.max_lag; 
lags_meter_rel          = parser.Results.lags_meter_rel; 
lags_meter_unrel        = parser.Results.lags_meter_unrel; 
features                = parser.Results.features; 
prec                    = parser.Results.prec; 

col_acf         = parser.Results.col_acf; 
col_meter_rel   = parser.Results.col_meter_rel; 
col_meter_unrel = parser.Results.col_meter_unrel; 
col_ap          = parser.Results.col_ap; 
linew           = parser.Results.linew; 
linew_ap        = parser.Results.linew_ap; 
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

if ~isempty(ap)
    plot(ax, lags, ap, '--', 'marker', 'none', 'linew', linew_ap, 'color', col_ap)
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












