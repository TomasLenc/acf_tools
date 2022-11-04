function r = correlation(x1, x2)
% Pearson correlationfor multi-dimensional arrays. Time must be last dimension!

nd = ndims(x1); 

x1_norm = x1 - mean(x1, nd); 
x2_norm = x2 - mean(x2, nd); 

num = sum(x1_norm .* x2_norm, nd); 

den = sqrt(sum(x1_norm.^2, nd) .* sum(x2_norm.^2, nd)); 

r = num ./ den; 