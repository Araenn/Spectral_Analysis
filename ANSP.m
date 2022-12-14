clc; clear; close all

N = 2048;
nu = 0.2;
sigma = 1;
bruit = sigma*randn(N, 1);


%% correlogramme
% k = (1:N)';
%R_est = autocorr(bruit);

for k = 2:N %on ne calcule pas le 0
    for n = k + 1:N
        R_est(k) = (1/(N-k)) * (bruit(n)* conj(bruit(n - k)));
    end
end

R_est(1) = (1/(N)) * conv(bruit(1), conj(bruit(1))); %calcul du 0
R_est2 = [fliplr(R_est)';R_est'];
R_est = [R_est';fliplr(R_est)'];

S_est = fft(R_est2);

figure(1)
stem(bruit)
grid()
title("Bruit")

n_aff = (-N+1:N-2)';
figure(2)
stem(n_aff, R_est)
grid()
title("Auto correlation estimee du bruit")



figure(3)

plot(S_est)
grid()
title("Correlogramme")

%faire correlog, peridogramme, periodogramme moyenne, puis moyenn et
%fenetre