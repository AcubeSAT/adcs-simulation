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

%Kp_gain = 1e-05*diag([0.580635889546100,91.5110525189067,65.7001327677557]);%[1]
%Kd_gain=1e-03*diag([2.35107597896066,0.910183539524328,26.3068606513231]);%[2]

%Kp_gain = 1e-05*diag([76.6334486792159,72.5325697578856,55.1500663058936]);%[3]
%Kd_gain= 1e-03*diag([10.0479528568468,25.9061817478587,46.4449268849864]);%[4]

%Kp_gain= 1e-05*diag([39.1103350435605,93.6662854162826,54.5051339290374]);%[5]
%Kd_gain=1e-03*diag([81.6378389258027,61.5756601133430,96.7727607126781]);%[6]

%Kp_gain = 1e-05*diag([79.6779166239878,79.4583028666047,9.62090934528343]);%[7]
%Kd_gain=1e-03*diag([1.99283330685127,42.2335876399571,57.1276981902645]);%[8]

%Kp_gain = 1e-05*diag([191.501367086860,186.798649551510,138.965724595163]);%[9]
%Kd_gain= 1e-03*diag([0.612173320368121,55.4550153997157,52.6266508247473]);%[10]

%Kp_gain =1e-05*diag([482.444267599638,396.103664779777,196.238509767084]);%[11]
%Kd_gain=1e-03*diag([347.508061487909,244.757197894116,59.4988407791883]);%[12]
 
%Kp_gain = 1e-05*diag([365.701788816455,475.667465612545,34.0668096044233]);%[13]
%Kd_gain=1e-03*diag([337.574468207506,77.5843516315163,103.406020943830]);%[14]

%Kp_gain =1e-05*diag([20.2587831878118,286.255344211728,100.823017111627]);%[15]
%Kd_gain=1e-03*diag([392.756598662821,389.122258132101,322.639766559320]);%[16]

%Kp_gain =1e-05*diag([291.727298024707,271.447126508074,185.896062600910]);%[17]
%Kd_gain=1e-03*diag([395.932955671478,67.3307615687087,291.407470833244]);%[18]

%Kp_gain =1e-05*diag([34.5589711718382,63.5290170868800,25.7857122971642]);%[19]
%Kd_gain=1e-03*diag([68.4378042020425,30.6513593148011,122.622680062188]);%[20]

%Kp_gain =1e-05*diag([55.6996437734097,191.898485278581,169.825861173755]);%[23]
%Kd_gain=1e-03*diag([34.2373375623124,6.88921610058175,159.164980227413]);%[24]

%Kp_gain =1e-05*diag([148.88,200.79,103.3]);%[25]
%Kp_gain =1e-05*diag([62.24,148.55,93.3]);%[26]
%Kd_gain=1e-03*diag([152.1,91.1,121.63]);%[27]
%Kd_gain=2e-03*diag([28.5,40.94,80.43]);%[28]

%Kp_gain =1e-05*diag([95.4475076777000,97.3909244779433,68.2295286594339]);%[29]
%Kd_gain=1e-03*diag([84.9657965004803,23.3099746930591,40.43]);%[30]

%Kp_gain =1e-05*diag([98.2626, 135.7470, 165.0197]);%[31]
%Kp_gain =1e-05*diag([98.2626, 115.7470, 95.0197]);%[32]
%Kd_gain=1e-03*diag([37.4058, 55.2206, 51.2378]);%[33]

%Kp_gain =2e-05*diag([62.1981870797147,91.1805018033854,29.7882136673452]);%[34]
%Kd_gain=2e-03*diag([2.59099528274315,18.3838935931316,27.6766211505791]);%[35]

%Kp_gain =2e-05*diag([85.1629305868777,70.6671088019609,34.1635726666133]);%[36]
%Kd_gain=2e-03*diag([19.4095250431208,13.1156208473730,99.6134716626885]);%[37]

%Kp_gain =2e-05*diag([63.3296890801963,45.2937001396341,8.68301662731753]);%[38]
%Kd_gain=2e-03*diag([8.88492529288979,12.5194274408755,79.5951730863370]);%[39]

%Kp_gain =2e-05*diag([70.24,69.06,24.23]);%!!!%[40]
%Kd_gain=2e-03*diag([10.3,14.7,68.97]);%!!!%[41]

%Kp_gain =2e-05*diag([62.75,80.95,32]);%!!!%[42]
%Kd_gain=2e-03*diag([10.3,14.7,53.65]);%!!!%[43]

%Kp_gain =2e-05*diag([62.75,82.95,32]);%!!!%[44]
%Kd_gain=2e-03*diag([10.3,14.7,48.65]);%!!!%[45]

%Kp_gain =2e-05*diag([66,85.9,33.55]);%!!!%[48]
%Kd_gain=2e-03*diag([10.4,10.55,46.6]);%!!!%[49]

%Kp_gain =2e-05*diag([69.2,88.9,35.1]); %[50]
%Kd_gain= 2e-03*diag([10,6.4,44.5]); %[51]

%Kp_gain =2e-05*diag([69.1820271090564,88.8470310647548,35.1309716222657]);%[52]
%Kd_gain= 2e-03*diag([10.0029392633184,6.43624923102076,44.4898163655996]);%[53]

%Kp_gain =2e-05*diag([46.96,86.16,43.58]);%[54]
%Kd_gain=2e-03*diag([2.81,11.41,51.46]);%[55]

%Kp_gain=2e-05*diag([63.3296890801963,2*43.8207019007756,4*7.68301662731753]);%[56]
%Kd_gain=2e-03*diag([8.88492529288979,9.57499487347171,74.3018928206815/2]);%[57]

%Kp_gain=2e-05*diag([87.2734236968545,84.6175776618871,18.9989196044459]);%[58]
%Kd_gain=2e-03*diag([15.0088177899553,1.82129859851249,47.9592654623985]);%[59]

%Kp_gain =1e-05*diag([55.6996437734097,191.898485278581,199.825861173755]);%[60]
%Kd_gain=1e-03*diag([34.2373375623124,6.88921610058175,189.164980227413]);%[61]

%Kp_gain =1e-05*diag([55.7,191.9,199.8]);%[62]
%Kd_gain=1e-03*diag([34.2,6.9,189.15]);%[63]

%Kp_gain =1e-05*diag([36.5589711718382,63.5290170868800,15.7857122971642]);%[64]
%Kd_gain=1e-03*diag([69.4378042020425,30.6513593148011,122.622680062188]);%[65]

%Kp_gain =2e-05*diag([72.1981870797147,101.1805018033854,39.7882136673452]);%[66]
%Kd_gain=2e-03*diag([12.59099528274315,28.3838935931316,37.6766211505791]);%[67]

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

Kp_gain = 5e-03*diag([3,3,1]);%[76] %!!! Main one !!!%
Kd_gain = 1e-01*diag([5,5,1]);%[77]
% 
% Kp_gain= 1e-05*diag([11.3 38.5 29.7]);
% Kd_gain= 1e-04*diag([9.6 47.6 89.6]);

% 
% Kp_gain= 1e-05*diag([30 38.5 40]);
% Kd_gain= 1e-04*diag([9.6 47.6 70]);

%Kp_gain = 1e-05*diag([79.6779166239878,79.4583028666047,9.62090934528343]);%[78]
%Kd_gain = 1e-03*diag([2.99283330685127,43.2335876399571,57.1276981902645]);%[79]

%Kp_gain = 1e-05*diag([191.501367086860,186.798649551510,138.965724595163]);%[80]
%Kd_gain = 1e-03*diag([0.612173320368121,56.4550153997157,52.6266508247473]);%[81]

[instant_error_perform, Time] = Nominal_Simulation(Kp_gain, Kd_gain);
