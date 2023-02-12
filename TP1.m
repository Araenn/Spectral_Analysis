clc; clear; close all

N = 2048;
f = -0.5:1/N:0.5 - 1/N; %axe des frequences
n_aff = (0:N-1)'; %axe des echantillons (temps numerique)

sigma = 1;
bruit = sigma*randn(N, 1);


%% correlogramme
R_est = autocorrel(bruit, N);
S_est = fft(R_est); %correlogramme

figure(1)
stem(bruit, '.')
grid()
title("Bruit")

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
a2 = 0.5;
f1 = 0.05; %Hz
f2 = 0.06;

s1 = a1 * cos(2 * pi * f1 * n_aff);
s2 = a2 * cos(2 * pi * f2 * n_aff);

signal = bruit + s1 + s2;

p_signal = periodogramme(signal, N);
p_signalFenetre = periodogramme(signal .* hann(N), N);

T = 100;
fenetreRect = zeros(N, 1);
fenetreRect(N/2-T:N/2+T) = 1;

p_signalFenRect = periodogramme(signal .* fenetreRect, N);

nb_blocs = 2; %nombre de blocs, plus le nombre est petit, plus le perio est precis
somme = 0;
somme2 = 0;
for i = 1:nb_blocs
    p_signalMoyenne = periodogramme(signal(N/nb_blocs * (i - 1) + 1:N/nb_blocs * i), N/nb_blocs);
    p_signalMoyenneFenetre = periodogramme(signal(N/nb_blocs * (i - 1) + 1:N/nb_blocs * i) .* hann(floor(N/nb_blocs)), N/nb_blocs);
    somme = somme + p_signalMoyenne;
    somme2 = somme2 + p_signalMoyenneFenetre;
end
p_signalMoyenne = somme;
p_signalMoyenneFenetre = somme2;


figure(5)
subplot(4, 1, 1)
plot(f, p_signal)
grid()
title("Periodogramme du bruit blanc + 2 cos")
subplot(4,1,2)
plot(f, p_signalFenetre)
hold on
plot(f, p_signalFenRect)
grid()
title("Periodogramme fenetre Hanning et rectangle")
legend("Hanning", "rectangle")
subplot(4,1,3)
fe = -0.5: 1/length(p_signalMoyenne):0.5-1/length(p_signalMoyenne);
plot(fe, p_signalMoyenne)
grid()
title("Periodogramme moyenne")
subplot(4,1,4)
fe = -0.5: 1/length(p_signalMoyenneFenetre):0.5-1/length(p_signalMoyenneFenetre);
plot(fe, p_signalMoyenneFenetre)
grid()
title("Periodogramme moyenne fenetre")

%tests des eqm
nb_tests = 200;
eqm_simple = zeros(nb_tests, 1);
eqm_fenetre = zeros(nb_tests, 1);
eqm_moyenne = zeros(nb_tests, 1);
eqm_moyenneFenetre = zeros(nb_tests, 1);
for k = 1:nb_tests
    bruit = sigma*randn(N, 1);
    signal = bruit + s1 + s2;
    p_signalFenRect = periodogramme(signal .* fenetreRect, N);
    p_signal = periodogramme(signal, N);
    p_signalFenetre = periodogramme(signal .* hann(N), N);
    nb_blocs = 2; %nombre de blocs, plus le nombre est petit, plus le perio est precis
    somme = 0;
    somme2 = 0;
    for j = 1:nb_blocs
        p_signalMoyenne = periodogramme(signal(N/nb_blocs * (i - 1) + 1:N/nb_blocs * i), N/nb_blocs);
        p_signalMoyenneFenetre = periodogramme(signal(N/nb_blocs * (i - 1) + 1:N/nb_blocs * i) .* hann(floor(N/nb_blocs)), N/nb_blocs);
        somme = somme + p_signalMoyenne;
        somme2 = somme2 + p_signalMoyenneFenetre;
    end
    p_signalMoyenne = somme;
    p_signalMoyenneFenetre = somme2;

    eqm_simple(k) = EQM(signal, p_signal, N);
    eqm_fenetre(k) = EQM(signal, p_signalFenetre, N);
    eqm_moyenne(k) = EQM(signal, p_signalMoyenne, N);
    eqm_moyenneFenetre(k) = EQM(signal, p_signalMoyenneFenetre, N);
end

figure(6)
plot(eqm_simple)
hold on
plot(eqm_fenetre)
plot(eqm_moyenne)
plot(eqm_moyenneFenetre)
grid()
title("Erreurs quadratiques moyennes sur 200 réalisations")
legend("Simple", "Fenetre", "Moyenne", "Moyenne fenetre")
%% Moindres Carrés
figure(7)

M = [sin(2*pi*f1*n_aff) cos(2*pi*f1*n_aff) sin(2*pi*f2*n_aff) cos(2*pi*f2*n_aff)];
theta_c = [inv(M'*M)*M'*signal]; 
y1 = M*theta_c;
a_c1 = sqrt(theta_c(1)^2 + theta_c(2)^2)
a_c2 = sqrt(theta_c(3)^2 + theta_c(4)^2)

plot(n_aff, y1)
hold on
plot(n_aff, signal-bruit)
grid()
title("Signaux retrouvé par MC sans bruit")
legend("Signal MC", "original sans bruit")

%% Spectre rationnel, filtre passe-bas
fc = 300; %frequence de coupure
fs = 1000; %frequence d'echantillonnage
gain = 3;
[B, A] = butter(6, fc/fs);

y = filter(gain * B, A, bruit);

figure(8)
subplot(4, 1, 1)
p_butter = periodogramme(y, N);
plot(f, p_butter)
grid()
title("Periodogramme du bruit filtré par un filtre de Butterworth passe-bas")

subplot(4, 1, 2)
p_butterFenetreRect = periodogramme(y .* fenetreRect, N);
p_butterFenetreHann = periodogramme(y .* hann(N), N);
plot(f, p_butterFenetreRect)
hold on
plot(f, p_butterFenetreHann)
title("Periodogramme butter avec fenetres rect et hanning")
grid()
legend("rect", "Hanning")

subplot(4, 1, 3)
nb_blocs = 2; %nombre de blocs, plus le nombre est petit, plus le perio est precis
somme = 0;
somme2 = 0;
for i = 1:nb_blocs
    p_butterMoyenne = periodogramme(y(N/nb_blocs * (i - 1) + 1:N/nb_blocs * i), N/nb_blocs);
    p_butterMoyenneHann = periodogramme(y(N/nb_blocs * (i - 1) + 1:N/nb_blocs * i) .* hann(floor(N/nb_blocs)), N/nb_blocs);
    somme = somme + p_butterMoyenne;
    somme2 = somme2 + p_butterMoyenneHann;
end
p_butterMoyenne = somme;
p_butterMoyenneHann = somme2;

fe = -0.5: 1/length(p_butterMoyenneHann):0.5-1/length(p_butterMoyenneHann);
plot(fe, p_butterMoyenne)
grid()
title("Butter moyenne")

subplot(4, 1, 4)
plot(fe, p_butterMoyenneHann)
grid()
title("Butter moyenne fenetre")

%tests des eqm
figure(9)
nb_tests = 200;

eqm_simple = zeros(nb_tests, 1);
eqm_fenetreRect = zeros(nb_tests, 1);
eqm_fenetreHann = zeros(nb_tests, 1);
eqm_moyenne = zeros(nb_tests, 1);
eqm_moyenneFenetreHann = zeros(nb_tests, 1);
eqm_moyenneFenetreRect = zeros(nb_tests, 1);
for k = 1:nb_tests
    bruit = sigma*randn(N, 1);
    y = filter(gain * B, A, bruit);
    p_butter = periodogramme(y, N);
    p_butterFenetreRect = periodogramme(y .* fenetreRect, N);
    p_butterFenetreHann = periodogramme(y .* hann(N), N);

    nb_blocs = 2; %nombre de blocs, plus le nombre est petit, plus le perio est precis
    somme = 0;
    somme2 = 0;
    somme3 = 0;
    for i = 1:nb_blocs
        p_butterMoyenne = periodogramme(y(N/nb_blocs * (i - 1) + 1:N/nb_blocs * i), N/nb_blocs);
        p_butterMoyenneHann = periodogramme(y(N/nb_blocs * (i - 1) + 1:N/nb_blocs * i) .* hann(floor(N/nb_blocs)), N/nb_blocs);
        p_butterMoyenneRect = periodogramme(y(N/nb_blocs * (i - 1) + 1:N/nb_blocs * i) .* fenetreRect(nb_tests), N/nb_blocs);
        somme = somme + p_butterMoyenne;
        somme2 = somme2 + p_butterMoyenneHann;
        somme3 = somme3 + p_butterMoyenneRect;
    end
    p_butterMoyenne = somme;
    p_butterMoyenneHann = somme2;
    p_butterMoyenneRect = somme3;

    eqm_simple(k) = EQM(signal, p_butter, N);
    eqm_fenetreRect(k) = EQM(signal, p_butterFenetreRect, N);
    eqm_fenetreHann(k) = EQM(signal, p_butterFenetreHann, N);
    eqm_moyenne(k) = EQM(signal, p_butterMoyenne, N);
    eqm_moyenneFenetreHann(k) = EQM(signal, p_butterMoyenneHann, N);
    eqm_moyenneFenetreRect(k) = EQM(signal, p_butterMoyenneRect, N);
end

plot(eqm_simple)
hold on
plot(eqm_fenetreRect)
plot(eqm_fenetreHann)
plot(eqm_moyenne)
plot(eqm_moyenneFenetreHann)
plot(eqm_moyenneFenetreRect)
grid()
title("Erreurs quadratiques moyennes sur 200 réalisations")
legend("Simple", "Fenetre rect","Fenetre Hann", "Moyenne", "Moyenne fenetre rect", "Moyenne fenetre hann")

%% Methode Yule-Walker
theta = [1 -1.5 0.7];
bruit = sigma*randn(N, 1);
s = filter(1, theta, bruit);
r = xcorr(s);
r = r';
R = r(N:2*N - 1);

rxm = [R(2); R(3)];
Rxm = [R(1) R(2); R(2) R(1)];

theta_est = inv(-Rxm) * rxm