close all, clear all, clc

%% Description:
%{
This codes implements the paper Simulation of Nonstationary Stochastic 
Processes by Spectral Representation, Liang 2007.

Here the Evolutionary-PSD that describes the N-S component of El Centro 
earthquake is assessed with 4 methods, then the spectral representation is
formulated and sample functions are generated, with the PSD estimated with
STFT, in order to compute the mean acceleration Spectra of SDOF.

Bibliography:

-    Simulation of Nonstationary Stochastic Processes by Spectral
     Representation, Liang 2007

Made by: Juan Sebastián Delgado Trujillo
%}

%% Time series:

% the N-S component of El Centro earthquake is loaded
eq   = load('elcentro_NS.dat');
time = eq(:, 1);        
data = eq(:, 2);

ts      = time(2) - time(1);    % sampling period
dat_ene = trapz(ts, data.^2);   % earthquake energy eq.92 (Parseval's
                                % theorem

%% PSD estimation: 
%{
The PSD is estimated using the STFT by 2 means, matlab's function
SPECTROGRAM and my_PSD. Also, in the Wavelets transform the PSD will be
assumed equal to the scalogram and discrete wavelet spectre; in the
Hilbert - Huang transform the produced spectre is assumed as the PSD; and
in Wigner-Vile the distibution is assumed as the PSD
%}

% this window is created (following the paper a Gaussian window with an
% standard deviation of 0.25 is used) for the STFT method
win_n = 128;
s_dev = 0.25;
my_gw = @(sd, L)  (1/(sqrt(2*pi)*sd))*exp(-0.5*(((1:L) - L/2)/(sd*L/2)).^2);    % eq. 94
win   = my_gw(s_dev, win_n);
Knorm = trapz(ts, win.^2);     % energy of the window
win   = win/sqrt(Knorm);       % normalization to satisfy eq. 93
jump  = 25;                    % ammount of non-overlapping samples

% Matlab's Spectogram:
[aux, freq_1, t_1, psd_1] = spectrogram(data, win, win_n - jump, [], 1/ts);

% my Spectogram:
[psd_2, freq_2, t_2]      = my_psd(data', win, jump, ts);

% wavelets:
level = 8;            wpt   = wpdec(data, level, 'sym6');
[psd_3a, t_3a, freq_3a] = wpspectrum(wpt, 1/ts);      % DW spectre
[psd_3b, freq_3b]       = cwt(data, 'amor', 1/ts);    % scalogram      
t_3b  = time';        psd_3b = abs(psd_3b).^2;

% Hilbert -Huang transform:
imf                  = emd(data);
[psd_4, freq_4, t_4] = hht(imf, 1/ts);

% Wigner-Vile distribution:
[psd_5, freq_5, t_5] = wvd(data, 1/ts, 'smoothedPseudo');

figure
%subplot(3, 2, 1)
mesh(t_1, freq_1, psd_1)
axis tight
title 'PSD (Matlab spectrogram)'
xlim([0, 14])
ylim([0, 10])
xlabel 'time (s)'
ylabel 'Frequencies (Hz)'
colorbar

subplot(3, 2, 2)
contourf(t_2, freq_2, psd_2)
axis tight
title 'PSD (my PSD)'
xlim([0, 14])
ylim([0, 10])
xlabel 'time (s)'
ylabel 'Frequencies (Hz)'
colorbar

subplot(3, 2, 3)
contourf(t_3a, freq_3a, psd_3a)
axis tight
title ' PSD (discrete wavelet spectre)'
xlim([0, 14])
ylim([0, 10])
xlabel 'time (s)'
ylabel 'Frequencies (Hz)'
colorbar

subplot(3, 2, 4)
contourf(t_3b, freq_3b, psd_3b)
axis tight
xlim([0, 14])
ylim([0, 10])
title 'PSD (scalogram)'
xlabel 'time (s)'
ylabel 'Frequencies (Hz)'
colorbar

subplot(3, 2, 5)
contourf(t_4, freq_4, psd_4)
axis tight
title 'PSD (HHT spectre)'
xlim([0, 14])
ylim([0, 10])
xlabel 'time (s)'
ylabel 'Frequencies (Hz)'
colorbar

subplot(3, 2, 6)
contourf(t_5, freq_5, psd_5)
axis tight
title 'PSD (WV distribution)'
xlim([0, 14])
ylim([0, 10])
xlabel 'time (s)'
ylabel 'Frequencies (Hz)'
colorbar

%% spectral representation:

% earthqueakes are simulated using the PSD, assessed with my_PSD and
% spectograms function, with the spectral representation, i. e. multiples
% realizations are performed
[ft_s1, t_s1] = spectral_simu(psd_1, freq_1, [t_1, max(time)]);
[ft_s2, t_s2] = spectral_simu(psd_2, freq_2, [t_2, max(time)]);

y_lim = max([max(abs(data)), max(abs(ft_s1)), max(abs(ft_s2))]);

figure

subplot(3, 1, 1)
plot(time, data)
axis tight
grid on
title 'El Centro N-S component'
xlabel 'time (s)'
ylabel 'Aceleration (m/s²)'
ylim([-y_lim, y_lim])

subplot(3, 1, 2)
plot(t_s1, ft_s1, 'r')
axis tight
grid on
title 'Simulation with Matlabs spectogram'
xlabel 'time (s)'
ylabel 'Aceleration (m/s²)'
ylim([-y_lim, y_lim])

subplot(3, 1, 3)
plot(t_s2, ft_s2, 'r')
axis tight
grid on
title 'Simulation with my PSD'
xlabel 'time (s)'
ylabel 'Aceleration (m/s²)'
ylim([-y_lim, y_lim])

%% spectral energy: 
% in the PSDs computed with STFT, the total energy of the time - frequency
% domain is calculated using the Parseval's theorem (eq. 92)

% Matlab's spectogram:
dw_1   = freq_1(2) - freq_1(1);
dt_1   = t_1(2) - t_1(1);
sp_en1 = trapz(dt_1, trapz(dw_1, psd_1, 1));

% my_PSD: 
dw_2   = freq_2(2) - freq_2(1);
dt_2   = t_2(2) - t_2(1);
sp_en2 = trapz(dt_2, trapz(dw_2, psd_1, 1));

fprintf(['The seismic energy is %.3f ' ...
         '\nSTFT: The energy of the PSD is %.3f , %.3f of the energy'...
         ' measured by the data'...
         '\nmy STFT: The energy of the PSD is %.3f, %.3f of the energy'...
         ' measured by the data'], dat_ene, sp_en1, ...
         100*sp_en1/dat_ene, sp_en2, 100*sp_en2/dat_ene)

%% variance of the process:

% a PSD of the process is performed using a unit rectangular window with a
% duration of 1 second
jump_v   = 1;
recw_n   = 1/ts;
[psd_var, freq_var, t_var] = my_psd(data', rectwin(recw_n), jump_v, ts);

% to assess the variance through time, in each instant a integral in the
% whole frequency domain is calculated 
dw    = freq_var(2) - freq_var(1);
var_t = trapz(dw, psd_var, 1);

% the variance in the same times are calculated too in the earthquake data
[var_2, tvar_2] = data_var(data', recw_n, jump_v, ts);

figure
hold on
plot(t_var, var_t)
plot(tvar_2, var_2, 'r--')
legend ('VAR[x] in PSD', 'VAR[x] in  data record')
axis tight
grid on
title 'Variance evolution through time'
xlabel 'time (s)'
ylabel 'Variance (m²/s⁴)'
     
%% SDOF's response:
%{
in order to estimate the accuracy of the spectral representation, in the
structure's response, multiples realizations of the quake are performed and
the mean acceleration spectre is computed and compared with the El CENTRO's
spectre
%}

% period interval
p_min = 0.2;        p_max = 5;

% acceleration spectre of the el centro earthquake
[T_ec, maxa_ec] = acc_spec(data, p_min, p_max, ts);


n_simu = 25;      % ammount of quake simulations
simu_e = cell([1, n_simu]);

% for each simulation:
for i = 1:25
    % 1) a realization of the process is performed
    [ft_ac, aux] = spectral_simu(psd_2, freq_2, [t_2, max(time)]);  
    % 2) the respective acceleration spectre is generated
    [T_sm, simu_e{1, i}] = acc_spec(ft_ac, p_min, p_max, aux(2) - aux(1));
end

% in each instant the mean of all spectre is taken
simu_e    = cell2mat(simu_e);
simu_spec = mean(simu_e, 2);

figure
hold on
plot(T_ec, maxa_ec)
plot(T_sm, simu_spec, 'r --')
legend('El Centro earthquake Spectre', 'Mean spectre of simulations')
hold on
axis tight
grid on
title 'Accelartion Spectre'
xlabel 'Period (s)'
ylabel 'Acceleration (g)'