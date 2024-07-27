%% Script for modeling noise for the MEMS gyroscope
%
%  Inputs        :
%  old_bias      - current bias
%  dt            - sampling interval
%  sigma_u       - bias instability std deviation
%  sigma_v       - gaussian noise std deviation
%
%  Outputs       :
%  total_noise   - noise
%  new_bias      - updated bias
%
%  Theoratical values :
%  sigma_u = 7.7570e-05
%  sigma_v = 0.0026
%
%
%  Based on Markley-Crassidis p.147

function [total_noise, new_bias] = gyro_noise_func(old_bias, dt, sigma_u, sigma_v)

    new_bias = old_bias + sigma_u * sqrt(dt) * randn(3, 1);

    total_noise = 0.5 * (old_bias + new_bias) + sqrt((sigma_v^2 / dt + sigma_u^2 * dt / 12)) * randn(3, 1);

end
