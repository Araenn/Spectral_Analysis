function e = EQM(signal, estimation, N)    
    if (length(estimation) ~= length(signal))
        estimation = [estimation; zeros(length(signal)-length(estimation), 1)];
    end
    e = esperance(abs(estimation - signal).^2, N);
end