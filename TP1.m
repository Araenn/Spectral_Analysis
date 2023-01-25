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

%% peridogramme

perido = (1/N) * (abs(fft(bruit)).^2);
figure(4)
plot(perido)
grid()
title("Periodogramme du bruit")
%% debug avec fonctions matlab
% figure(4)
% [rho, lags] = xcorr(bruit);
% stem(lags, rho)
% title("fonction xcorr")
% 
% figure(5)
% S = fft(rho(N:end));
% plot(abs(S))

%% eqm
figure(5)
hold on
for i = 1:10
    newBruit = randn(N, 1) * sigma;
    newR_est = autocorrel(newBruit, N);
    newS_est = fft(newR_est);
    newPeriod = periodogramme(newBruit, N);
    
    %faire eqm
    e = newPeriod - newS_est';
    eq = abs(e).^2;
    eqm = eq/length(eq);
    plot(eqm)
end
grid()