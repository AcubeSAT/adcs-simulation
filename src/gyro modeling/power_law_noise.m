% This function generates power law noise with a specified exponent (beta),
% and scaling parameter(v0) , in a specified dimensional space(dims). The
% power law scaling is applied in the frequency domain , followed by an
% inverse FFT to generate the time-domain noise signal . The noise is then
% normalized based on the scaling parameters . 
%
%
% INPUTS:
% - beta: The exponent of the power law .
% - v0: The minimum frequency  in the range [0, 0.5].
% - dims: A vector defining the dimensions of the generated noise, where 
%         dims(end) specifies the length of the time-domain signal.
%
% OUTPUTS:
% - noise: The generated power law noise signal in the time domain.
%
% AUTHOR:
% Ujjwal Panda (ujjwalpanda97@gmail.com)
%
% LICENSE:
% MIT License
%
% Copyright (c) 2022-2023 Ujjwal Panda <ujjwalpanda97@gmail.com>
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%


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