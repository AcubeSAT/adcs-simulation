function [White_Noise,Pink_Noise,Red_Noise,total_noise] = calculate_gyro_noise(ARW,RRW,BI,N)

v0 = 0.1;

White_Noise = ARW * power_law_noise(0,v0,N);
Pink_Noise = BI* power_law_noise(1,v0,N);
Red_Noise = RRW * power_law_noise(2,v0,N);

total_noise= White_Noise + Pink_Noise + Red_Noise;

end









