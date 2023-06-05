function tests = test_lag_selection
    tests = functiontests(localfunctions);
end


function test_get_lag_harmonics_single(test_case)
    
    max_lag = 11; 
    lags = get_lag_harmonics(2, max_lag); 
    
    assert(isequal(lags, [2,4,6,8,10]))
    
end


function test_get_lag_harmonics_single_excl_single(test_case)
    
    max_lag = 11; 
    lags = get_lag_harmonics(2, max_lag, 'lag_harm_to_exclude', 3); 
    
    assert(isequal(lags, [2,4,8,10]))
    
end


function test_get_lag_harmonics_multi(test_case)
    
    max_lag = 11; 
    lags = get_lag_harmonics([2, 3], max_lag); 
    
    assert(isequal(lags, [2,3,4,6,8,9,10]))
    
end


function test_get_lag_harmonics_multi_excl_multi(test_case)
    
    max_lag = 17; 
    lags = get_lag_harmonics([2, 5], max_lag, 'lag_harm_to_exclude', [3, 7]); 
    
    assert(isequal(lags, [2, 4, 5, 8, 10, 16]))
    
end


function test_get_lag_harmonics_multi_excl_multi_noninteger(test_case)
    
    max_lag = 1.7; 
    lags = get_lag_harmonics([0.2, 0.5], max_lag, ...
                             'lag_harm_to_exclude', [0.3, 0.7]); 
    
    assert(isequal(lags, [.2, .4, .5, .8, 1.0, 1.6]))
    
end


