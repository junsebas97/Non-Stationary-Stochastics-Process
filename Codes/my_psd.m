function [psd, freq, time] = my_psd(x_i, window, jump, dt)
%{
This code assesss the power spectral density of a given time series using
the STFT algorithm. It returns the amplitudes, frequencies ans times

Inputs:
   x_i:    array containing the measurements
   window: array 
   jump:   number of no-overlapping samples between evaluations
   dt:     time between measurements

outputs:
   espec: values of the one-sided spectogram
   freq : frequencies that were considered
   time : time-line of the STFT
%}   
    %% Spectrogram: 
    
    % the size of the time series, the size of the window (size of the 
    % temporal evaluation) and the ammount of temporal fft
    n_xi   = length(x_i);
    ln_win = length(window);
    n_stft = floor((n_xi - ln_win)/jump);
    
    psd = cell([1, n_stft]);
    time = zeros([1, n_stft]);
    
    % for each temporal - serie:
    for i = 1:n_stft
        
        % the samples of the temporal signal are extracted
        x_temp = x_i(1, (i - 1)*jump + 1 : (i - 1)*jump + ln_win);
        
        % The temporal PSD is assessed and stored. The frequencies are
        % taken too
        [psd{1, i}, freq] = periodogram(x_temp, window, [], 1/dt);
        
        % and the time of the medium sample is stored
        time(1, i) = (ln_win/2 + (i - 1)*jump)*dt; 
    end
    
    psd = cell2mat(psd);
end