function R = autocorrel(bruit, N)
    for k = 0:N-1
        R(k+1) = sum( bruit(1+k:end) .* conj(bruit(1:end-k)) ) / (N); %autocorrelation
    end
end