% This function calculates the total gyro noise and the gyro noise components: White Noise,
% Pink Noise, and Red Noise, based on the given noise parameters (ARW, RRW, and BI). 
% The function uses the `power_law_noise` function to generate each noise type 
% with specific exponents (0 for White, 1 for Pink, and 2 for Red).
%
% INPUTS:
% - ARW: The Angular Random Walk (ARW) parameter, used to scale White Noise((rad/sec)/sqrt(Hz).
% - RRW: The Random Rate Walk (RRW) parameter, used to scale Red Noise(rad/sec)sqrt(Hz).
% - BI: The Bias Instability (BI) parameter, used to scale Pink Noise (rad/sec).
% - N: The dimensions of the output noise signal.
%
% OUTPUTS:
% - White_Noise: The generated White Noise signal, scaled by ARW.
% - Pink_Noise: The generated Pink Noise signal, scaled by BI.
% - Red_Noise: The generated Red Noise signal, scaled by RRW.
% - total_noise: The sum of White, Pink, and Red noise components.

function [White_Noise,Pink_Noise,Red_Noise,total_noise] = calculate_gyro_noise(ARW,RRW,BI,N)

v0 = 0.1;

White_Noise = ARW * power_law_noise(0,v0,N);
Pink_Noise = BI* power_law_noise(1,v0,N);
Red_Noise = RRW * power_law_noise(2,v0,N);

total_noise= White_Noise + Pink_Noise + Red_Noise;

end









