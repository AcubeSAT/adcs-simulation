function noise = power_law_noise(beta, v0, dims)
  
    n = dims(end);
    
    v = (0:floor(n/2)) / n;
       
  
    if ~(0 <= v0 && v0 <= 0.5)
        error('v0 must be in the range [0, 0.5]');
    end
    v0 = max(v0, 1/n);
    
    S = (max(v, v0)).^(-beta / 2);
    
 
    X = randn([dims(1:end-1), length(v)]);
    Y = randn([dims(1:end-1), length(v)]);
    
   
    Y(:, 1) = 0; 
    if mod(n, 2) == 0
        Y(:, end) = 0; 
    end
    
    
    F = S .* (X + 1i * Y);
    
    noise = real(ifft(F, n, numel(dims)));
   
end


