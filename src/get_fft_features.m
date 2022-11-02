function feat = get_fft_features(mX, freq, freq_meter_rel, freq_meter_unrel)

if ~isrow(freq)
    freq = freq';
end
if ~isrow(freq_meter_rel)
    freq_meter_rel = freq_meter_rel';
end
if ~isrow(freq_meter_unrel)
    freq_meter_unrel = freq_meter_unrel';
end

idx_meter_rel = dsearchn(freq', freq_meter_rel')'; 
idx_meter_unrel = dsearchn(freq', freq_meter_unrel')'; 

index = cell(1, ndims(mX));
index(:) = {':'};
index{end} = [idx_meter_rel, idx_meter_unrel];
amps = mX(index{:});

z = zscore(amps, [], ndims(mX)); 

index = cell(1, ndims(mX));
index(:) = {':'};
index{end} = [1 : length(idx_meter_rel)];
z_meter_rel = mean(z(index{:}), ndims(mX)); 

feat = []; 
feat.z_meter_rel = z_meter_rel; 
