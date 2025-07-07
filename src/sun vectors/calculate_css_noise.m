% This script estimates the total noise current of a coarse sun sensor.
% The total noise current includes:
%   - Thermal noise from the feedback resistor
%   - Shot noise from photocurrent and dark current
%   - Quantization noise due to ADC resolution
%
% Reference: https://essay.utwente.nl/72092/1/Hollander_BA_EWI.pdf



k = 1.380649e-23;  % Boltzmann constant [J/K]
T = 298;           % Temperature [K]
B = 202;           % Bandwidth [Hz]
R = 3970;          % Feedback resistor [Ohm]
q = 1.602e-19;     % Electron charge [C]
Ip = 170e-6;       % Photocurrent [A]
Id = 1.7e-6;       % Dark current [A]
LSB = 0.7e-3;      % ADC LSB [V]

Thermal_noise= sqrt(2*k*T/R); % [A]
Photocurrent_noise=sqrt(2*q*Ip*B); %[A]
Dark_noise=sqrt(2*q*Id*B); %[A] 
Quantization_noise=(1/sqrt(12))*(LSB/R); %[A] 

total_noise=sqrt(Thermal_noise^2+Photocurrent_noise^2+Dark_noise^2+Quantization_noise^2); %[A]

disp(total_noise);