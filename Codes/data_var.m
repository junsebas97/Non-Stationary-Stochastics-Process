function [var_t, time] = data_var(x_i, win_len, jump, dt)
%{
This function calcules the evolution of the signal's variance through time
    
Inputs:
    x_i:     signal record
    win_len: size of the interval
    jump:    number of no-overlapping samples between evaluations
    dt:      time between measurements
    
Outputs:
    var_t: variance evaluations
    time:  instants when the variance is assessed
%}

    n_xi   = length(x_i);
    n_eval = floor((n_xi - win_len)/jump);    % ammount of variance
                                              % evaluations  
    var_t = nan([1, n_eval]);
    time  = nan([1, n_eval]);
    
    for i = 1:n_eval
        
        % the samples of the temporal signal are extracted
        x_temp = x_i(1, (i - 1)*jump + 1 : (i - 1)*jump + win_len);
        
        % and the temporal variance is computed
        var_t(1, i) = var(x_temp);      
        time(1, i)  = (win_len/2 + (i - 1)*jump)*dt; 
    end
end