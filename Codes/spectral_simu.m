function [ft, t_i] = spectral_simu(psd, freq, time)
%{
This function produces realizations, using an spectral representation, of a
non-stationary process described by the given evolutionary PSD
    
Inputs:
    PSD:  matrix containing the PSD evaluations in time (columns) and
          frequencies (rows)
    freq: a row vector containing the frequencies where PSD was evaluated
    time: a column vector the instants when the PSD was evaluated

outputs:
    ft:  magnitudes of the process realization through time
    t_i: evaluation's instants of the realization
    
%}
    %% time discretization:
    
    % between the PSD evaluation times the process is also evaluated, in
    % any interval the process is modeled as a stationary process 
    ti_be = 10;                    % ammount of stationary time instants
    n_ti  = length(time)*ti_be;    % ammount of process evaluation
    t_i   = linspace(0, max(time), n_ti);
    ft    = nan([n_ti, 1]);
    

    n_eval = size(psd, 2);         % ammount of time evaluation of the PSD
    
    %% representation's parameters:
    
    wn     = 2*pi*freq;            % considered circular frequencies
    dw     = freq(2) - freq(1);
    N      = size(freq, 1);
    phin   = 2*pi*rand(1, N);      % random phases of the process page 3
    
    %% realization:
    
    % for each PSD time evaluation:
    for ee = 0:n_eval
        
        % 1) the average PSD of the adjacent evaluations is computed
        if ee == 0
            spect_aux = (psd(:, ee + 1) + 0)./2;    % here is assumed that
                                                    % in t = 0 PSD = 0
        elseif ee == n_eval
            spect_aux = (psd(:, ee) + 0)./2;        % here is  assumed that
                                                    % in t = T PSD = 0
        else
            spect_aux = (psd(:, ee + 1) + psd(:, ee))./2;    
            
        end
        
        % 2) the magnitude of cosines are computed
        bn_t = sqrt(2*spect_aux*dw);    % eq. 39
        
        % 3) the interval is refined and in each instant the realization is
        %    produced
        for st = 1:ti_be          
            id     = ee*ti_be + st;     % row index
            ft(id) = sqrt(2)*sum(bn_t'.*cos(wn'*t_i(id) + phin));  % eq. 52
        end   
    end        
end