function tests = test_index_finding
    tests = functiontests(localfunctions);
end


function test_find_idx_tol(test_case)
    
    arr = [0.00001, 1, 1.999, 2.000001, 3, 4, 5, 6.00001, 10]; 
    
    idx = find_idx_tol(arr, [2, 4, 5.000001, 8]);
    assert(isequal(idx, [6]));
    
    idx = find_idx_tol(arr,  [2, 4, 5.000001, 8], 'tol', 1e-3);
    assert(isequal(idx, [4, 6, 7]));
    
end


function test_find_idx_tol_transposed(test_case)
    
    arr = [0.00001, 1, 1.999, 2.000001, 3, 4, 5, 6.00001, 10]; 
    
    idx = find_idx_tol(arr,  [2, 4, 5.000001, 8]', 'tol', 1e-3);
    assert(isequal(idx, [4, 6, 7]));
    
    idx = find_idx_tol(arr', [2, 4, 5.000001, 8], 'tol', 1e-3);
    assert(isequal(idx, [4, 6, 7]));

end


function test_find_idx_tol_nothing_found(test_case)

    verifyError(test_case, ...
        @() find_idx_tol([1,2,3], [5]),...
        'find_idx_tol:NoIdxFound');
    
end


function test_find_idx_tol_bad_shape(test_case)

    verifyError(test_case, ...
        @() find_idx_tol([1,2,3], [5,5;6,6]),...
        'find_idx_tol:badShape');
    
    verifyError(test_case, ...
        @() find_idx_tol([5,5;6,6], [1,2,3]),...
        'find_idx_tol:badShape');
end