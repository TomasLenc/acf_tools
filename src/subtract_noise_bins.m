function mX_subtracted = subtract_noise_bins(mX, from, to)

mX_subtracted = mX; 

N = size(mX, ndims(mX)); 

for i_f=1:length(mX)
    
    idx_1 = max(i_f - to, 1); 
    idx_2 = max(i_f - from, 1); 
    idx_3 = min(i_f + from, N); 
    idx_4 = min(i_f + to, N); 

    index = cell(1, ndims(mX));
    index(:) = {':'};
    index{end} = [idx_1:idx_2, idx_3:idx_4];
    mean_noise = mean(mX(index{:}), ndims(mX)); 

    index = cell(1, ndims(mX));
    index(:) = {':'};
    index{end} = i_f;
    mX_subtracted(index{:}) = mX(index{:}) - mean_noise; 
end