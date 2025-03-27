function noise = power_law_noise(beta, v0, dims)
   
    n=dims(end);
    
    f = (0:floor(n/2))/n;
    f(1) = 1/n; 
    
    
    S = (max(f, v0)).^(-beta/2);
    
    
    X = randn(1,length(f)) .* S;
    Y = randn(1,length(f)) .* S;
    
   
    Y(1) = 0;                   
    if mod(n,2) == 0
        Y(end) = 0;            
    end
    
    
    F = X + 1i*Y;
    
    
    if mod(n,2) == 0
        F_full = [F, conj(F(end-1:-1:2))];
    else
        F_full = [F, conj(F(end:-1:2))];
    end
    

    noise = real(ifft(F_full));
    noise = noise/std(noise);  
end