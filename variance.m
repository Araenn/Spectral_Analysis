function v = variance(signal, N)
    esp = esperance(signal, N);
    v = sum((signal - esp).^2)/(N-1);
end