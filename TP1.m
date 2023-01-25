clc; clear; close all

N = 2048;
nu = 0.2;
sigma = 1;
bruit = sigma*randn(N, 1);


%% correlogramme
R_est = autocorrel(bruit, N);
S_est = fft(R_est); %correlogramme

figure(1)
stem(bruit, '.')
grid()
title("Bruit")

n_aff = (0:N-1)';
figure(2)
stem(n_aff, R_est, '.')
grid()
title("Auto correlation estimee du bruit")

figure(3)
plot(abs(S_est))
grid()
title("Correlogramme")

%faire correlog, peridogramme, periodogramme moyenne, puis moyenn et
%fenetre

%% periodogramme bruit blanc

p_bruitBlanc = periodogramme(bruit, N);
figure(4)
plot(p_bruitBlanc)
grid()
title("Periodogramme du bruit")




%% periodogramme bruit blanc + cos

a1 = 1;
a2 = 0.2;
f1 = 0.05; %Hz
f2 = 0.06;

s1 = a1 * cos(2 * pi * f1 * n_aff);
s2 = a2 * cos(2 * pi * f2 * n_aff);

signal = bruit + s1 + s2;

p_signal = periodogramme(signal, N);
p_signalFenetre = periodogramme(signal .* hann(N), N);

figure(5)
subplot(2, 1, 1)
plot(p_signal)
grid()
title("Periodogramme du bruit blanc + 2 cos")
subplot(2,1,2)
plot(p_signalFenetre)
grid()
title("Periodogramme fenetre")
%% debug avec fonctions matlab
% figure(4)
% [rho, lags] = xcorr(bruit);
% stem(lags, rho)
% title("fonction xcorr")
% 
% figure(5)
% S = fft(rho(N:end));
% plot(abs(S))
