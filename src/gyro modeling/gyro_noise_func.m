% Noise model for the MEMS gyroscope

% Inputs: current bias,bias instability,gaussian noise std deviation 
% Outputs:[noise,updated bias]

% sigma_u is bias instability std deviation and sigma_v is gaussian noise
% theoretical  sigma_u = 7.7570e-05 sigma_v = 0.0026
function [total_noise,new_bias] = gyro_noise_func(old_bias,dt,sigma_u,sigma_v)

new_bias = old_bias + sigma_u * sqrt(dt) * randn(3,1);

total_noise =  0.5 * (old_bias + new_bias) + ...   % calculate noise due to bias and wgn according to markley p.147
    sqrt((sigma_v^2 / dt + sigma_u^2 * dt / 12)) * randn(3,1); 

end
