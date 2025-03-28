function noise = power_law_noise(beta, v0, dims)
   
    n=dims(end);
    
    f = (0:floor(n/2))/n;
    if ~(0 <= v0 && v0 <= 0.5)
        error('v0 must be in the range [0, 0.5]');
    end
    v0 = max(v0, 1/n);
    
    
    S = (max(f, v0)).^(-beta/2);
    
    
    X = randn([dims(1:end-1), length(f)]) ;
    Y = randn([dims(1:end-1), length(f)]) ;
    
   
    Y(1,:) = 0;                   
    if mod(n,2) == 0
        Y(:,end) = 0;            
    end
    
     F = S .* (X + 1i * Y);
    
    
    if mod(n,2) == 0
        F_full = [F, conj(F(end-1:-1:2))];
    else
        F_full = [F, conj(F(end:-1:2))];
    end
    

    noise = real(ifft(F_full, n, numel(dims)));
    noise = noise/std(noise);  
end