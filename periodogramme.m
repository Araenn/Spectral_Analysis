function R = periodogramme(bruit, N)
    R = (1/N) * (abs(fft(bruit)).^2);
end