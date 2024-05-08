%% ======================================================================== %%
%   Function for testing the various gains.
%   Also used as the main function of the simulations.
%   Bdot-Simulation as well as Nominal-Simulation are called either
%       seperately or together in series.
%   The selected gains are used as inputs for the Nominal Simulation.
%   The Bdot Simulation requires no inputs
% ========================================================================= %

close all;
clc;
clear variables;

%% ================ Adaptive 6PM =========================

%Kp_gain = 1e-05*diag([98.2626, 115.7470, 95.0197]);%[68] %!!! Main one !!!%
%Kd_gain = 1e-03*diag([37.4058, 56.2206, 52.2378]);%[69]

%Kp_gain = 1e-05*diag([79.6779166239878,79.4583028666047,9.62090934528343]);%[70]
%Kd_gain = 1e-03*diag([1.99283330685127,42.2335876399571,57.1276981902645]);%[71]

%Kp_gain = 1e-05*diag([34.5589711718382,63.5290170868800,15.7857122971642]);%[72]
%Kd_gain = 1e-03*diag([68.4378042020425,30.6513593148011,122.622680062188]);%[73]

%% =============== IdealQ 6PM =============================
% 
%Kp_gain = 2e-05*diag([62.1981870797147,91.1805018033854,29.7882136673452]);%[74] %!!! Not so main one !!!%
%Kd_gain = 2e-03*diag([3.59099528274315,20.3838935931316,28.6766211505791]);%[75]

% Kp_gain = 5e-03*diag([3,3,1]);%[76] %!!! Main one !!!%
% Kd_gain = 1e-01*diag([5,5,1]);%[77]
% 
% Kp_gain= 1e-05*diag([11.3 38.5 29.7]);
% Kd_gain= 1e-04*diag([9.6 47.6 89.6]);

% 
% Kp_gain= 1e-05*diag([30 50 60]);
% Kd_gain= 1e-04*diag([9.6 47.6 70]);

%% =============== 11PM =============================

% Kp_gain= 1e-05*diag([20 120 120]);
% Kd_gain= 1e-04*diag([75 65 75]);

%% ============ Worst-Case Inertia gains ============
Kp_gain= 1e-05*diag([20 150 120]); % 500km 
Kd_gain= 1e-04*diag([75 100 75]);

% Kp_gain= 1e-05*diag([20 100 90]); % 600km 
% Kd_gain= 1e-04*diag([90 90 90]);


%%
%Kp_gain = 1e-05*diag([79.6779166239878,79.4583028666047,9.62090934528343]);%[78]
%Kd_gain = 1e-03*diag([2.99283330685127,43.2335876399571,57.1276981902645]);%[79]

%Kp_gain = 1e-05*diag([191.501367086860,186.798649551510,138.965724595163]);%[80]
%Kd_gain = 1e-03*diag([0.612173320368121,56.4550153997157,52.6266508247473]);%[81]

[instant_error_perform, Time, eclipse] = Nadir_Pointing_function(Kp_gain, Kd_gain);
