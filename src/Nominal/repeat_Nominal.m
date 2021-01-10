%% ==========Repeats of Nominal Simulation============= 
close all;
clc;
clear variables;

Kp_gain = 1e-05*diag([0.580635889546100,91.5110525189067,65.7001327677557]); % ADCS Sim
Kd_gain= 1e-03*diag([2.35107597896066,0.910183539524328,26.3068606513231]);

%Kp_gain= 1e-05*diag([39.1103350435605,93.6662854162826,54.5051339290374]); % 11-PM
%Kd_gain= 1e-03*diag([81.6378389258027,61.5756601133430,96.7727607126781]);
[instant_error_perform, Time] = Nominal_Simulation(Kp_gain, Kd_gain);
