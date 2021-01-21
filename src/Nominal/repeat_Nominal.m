%% ==========Repeats of Nominal Simulation============= 
close all;
clc;
clear variables;

%Kp_gain = 1e-05*diag([0.580635889546100,91.5110525189067,65.7001327677557]); [1] % ADCS Sim
%Kd_gain=1e-03*diag([2.35107597896066,0.910183539524328,26.3068606513231]); [2]

%Kp_gain = 1e-05*diag([76.6334486792159,72.5325697578856,55.1500663058936]);% [3]  11-PM
%Kd_gain= 1e-03*diag([10.0479528568468,25.9061817478587,46.4449268849864]); %[4]

%Kp_gain= 1e-05*diag([39.1103350435605,93.6662854162826,54.5051339290374]); [5]% 11-PM, ^2 
%Kd_gain=1e-03*diag([81.6378389258027,61.5756601133430,96.7727607126781]);[6]

%Kp_gain = 1e-05*diag([79.6779166239878,79.4583028666047,9.62090934528343]);%[7] % 11-PM, abs()
%Kd_gain=1e-03*diag([1.99283330685127,42.2335876399571,57.1276981902645]);%[8]

%Kp_gain = 1e-05*diag([191.501367086860,186.798649551510,138.965724595163]);%[9] % 6-PM, looks a bit better
%Kd_gain= 1e-03*diag([0.612173320368121,55.4550153997157,52.6266508247473]);%[10]

%Kp_gain =1e-05*diag([482.444267599638,396.103664779777,196.238509767084]);%[11]
%Kd_gain=1e-03*diag([347.508061487909,244.757197894116,59.4988407791883]);%[12]

%Kp_gain = 1e-05*diag([365.701788816455,475.667465612545,34.0668096044233]);%[13] no eclipse, error
%Kd_gain=1e-03*diag([337.574468207506,77.5843516315163,103.406020943830]);%[14]

%Kp_gain =1e-05*diag([20.2587831878118,286.255344211728,100.823017111627]);%[15]
%Kd_gain=1e-03*diag([392.756598662821,389.122258132101,322.639766559320]);%[16]

%Kp_gain =1e-05*diag([291.727298024707,271.447126508074,185.896062600910]);%[17]
%Kd_gain=1e-03*diag([395.932955671478,67.3307615687087,291.407470833244]);%[18]

%Kp_gain =1e-05*diag([34.5589711718382,63.5290170868800,15.7857122971642]);%[19]
%Kd_gain=1e-03*diag([68.4378042020425,30.6513593148011,122.622680062188]);%[20]

%Kp_gain =1e-05*diag([55.6996437734097,191.898485278581,169.825861173755]);%[23]
%Kd_gain=1e-03*diag([34.2373375623124,6.88921610058175,159.164980227413]);%[24]

%Kp_gain =1e-05*diag([148.88,200.79,93.3]);%[25]
%Kp_gain =1e-05*diag([62.24,148.55,93.3]);%[25]
%Kd_gain=1e-03*diag([152.1,91.1,121.63]);%[26]
%Kd_gain=2e-03*diag([28.5,40.94,80.43]);%[26]

%Kp_gain =1e-05*diag([95.4475076777000,97.3909244779433,68.2295286594339]);
%Kd_gain=1e-03*diag([84.9657965004803,23.3099746930591,40.43]);

%Kp_gain =1e-05*diag([98.2626, 135.7470, 165.0197]);
%Kp_gain =1e-05*diag([98.2626, 115.7470, 95.0197]);
%Kd_gain=1e-03*diag([37.4058, 55.2206, 51.2378]);

%Kp_gain =2e-05*diag([62.1981870797147,91.1805018033854,29.7882136673452]);
%Kd_gain=2e-03*diag([2.59099528274315,18.3838935931316,27.6766211505791]);

Kp_gain =2e-05*diag([85.1629305868777,70.6671088019609,34.1635726666133]);
Kd_gain=2e-03*diag([19.4095250431208,13.1156208473730,99.6134716626885]);


[instant_error_perform, Time] = Nominal_Simulation(Kp_gain, Kd_gain);
