clc; clear; close all

N = 2048;
nu = 0.2;
sigma = 1;
bruit = sigma*randn(N, 1);


%% correlogramme
for k = 0:N-1
    R_est(k+1) = sum(bruit(1+k:end) .* conj(bruit(1:end-k))) / N;
end

S_est = fft(R_est);

figure(1)
stem(bruit)
grid()
title("Bruit")

n_aff = (0:N-1)';
figure(2)
stem(n_aff, R_est)
grid()
title("Auto correlation estimee du bruit")

figure(3)
plot(abs(S_est))
grid()
title("Correlogramme")

%faire correlog, peridogramme, periodogramme moyenne, puis moyenn et
%fenetre

%% debug avec fonctions matlab
% figure(4)
% [rho, lags] = xcorr(bruit);
% stem(lags, rho)
% title("fonction xcorr")
% 
% figure(5)
% S = fft(rho(N:end));
% plot(abs(S))