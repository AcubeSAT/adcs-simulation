function Const = constants()

    global orbits;

    m = 4;       % Satellite mass [kg]
    lx = 0.1;    % X-axis length
    ly = 0.1;    % Y-axis length
    lz = 0.3405;    % Z-axis length

    N_coils = [400 400 800];                % Number of turns of the coil
    A_coils = [0.0022 0.0022 0.001566];  % Total area of the coil (m^2)
    R_coils = [110 110 31.32];              % Coil resistance (Ohm)
    mt = [5 5 1];                   % For torquer rods: mt = 1.66 * (length/diameter)^1.5.
                                            % For air core magnetorquers: mt = 1.
    W_up_limit = [0.2273 0.2273 0.8017];     % Upper limit of absolute power of mtq's (W)
    I_up_limit = [0.0455 0.0455 0.1566];     % Upper limit of absolute current of mtq's (A)
    
    K = eye(3) * diag([N_coils(1) * A_coils(1) / R_coils(1), N_coils(2) * A_coils(2) / R_coils(2), N_coils(3) * A_coils(3) / R_coils(3)]);
    K_inv = eye(3) / K;

    Km = 20e-5;       % Motor torque constant (RW)
    Rb = 10;           % Armature resistance (of RW)
    Ai = 1;             % The influence of the RW on the angular acceleration
                            % of the satellite (on the axis the RW is placed)
    Kv = 1/Km;        % Motor velocity constant (back-EMF constant)
    T_friction = 0.002; % The torque due to Friction Force
    b_friction = 9.5e-9;     % Viscuous friction
    c_friction = 1.9e-7;      % Coulomb friction
    Jw = 1.9e-6;         % Inertia of the RW
    A = 0.12;            % Used when desaturation of the RW
    lim_dz = 300;         % The absolute value of the limits of the deadzone (in rpm)
  
    p_400 = 7.55e-12;   % Atmospheric density [kg/m^3]
    p_500 = 1.80e-12;
    
    Ix = (m / 12) * (ly^2 + lz^2); % Xaxis inertia
    Iy = (m / 12) * (lx^2 + lz^2); % Yaxis inertia
    Iz = (m / 12) * (lx^2 + ly^2); % Zaxis inertia

%% Current Inertia
%     PMI = diag([0.03868845951 0.03899129965 0.00696263029]); %Principal Moments of Inertia
%     PAI = [-0.89 0.46 0; 0.46 0.89 -0.01; 0 0.01 1]; %Principal Axes of Inertia
%     Cm =  [0.03228 -0.02314 0.08244]'; % Center of mass
%% Modified Inertia
%     PMI = diag([0.03928052501	0.03948290041	0.00720041142]); %Principal Moments of Inertia
%     PAI = [-1 0.07 -0.01; 0.07 1 -0.02; -0.01 0.02 1]; %Principal Axes of Inertia
%     Cm = [0.03111 -0.02099 0.08135]';  % Center of mass

%% THIS INERTIA MAKES ME HAPPY
    PMI = diag([0.03552528444 0.03572444349 0.00626757327]); %Principal Moments of Inertia
    PAI = [0.98, 0.2, -0.02; -0.20, 0.98, 0.00 ; 0.02, 0.01, 1.00]; %Principal Axes of Inertia
    Cm =  [0.00121 0.00057 0.00188]'; % Center of mass



    for j=1:3
        PAI(:,j) = PAI(:,j)/norm(PAI(:,j));
    end
    I = PAI*PMI*PAI';
    I_inv = eye(3,3)/I;
       
    Re = 6371.2e3;                        % Earth radius [m]
    Rs = 500e3;                           % Satellite altitude [m]
    Radius = Re + Rs;                     % Distance from earth center to satellite [m]
    G = 6.67428e-11;                      % Earth gravitational constant
    M = 5.972e24;                         % Earth mass
    w_o = sqrt(G * M / Radius^3);         % Satellite angular velocity relative to Earth
    v_satellite = sqrt(G * M / Radius);   % Satellite's velocity in orbit
    w_o_io = [0 -w_o 0]';

    orbitPeriod = (2 * pi) / (w_o);       % Orbit Period
    n = orbitPeriod * orbits;             % Total simulation time
    mtq_max1 = 0.3;
    mtq_max2 = 0.3;
    mtq_max3 = 0.3;
    rw_max_torque = 1e-4;

%%  Passing the values of the parameters in a struct.

    Const.K = K;
    Const.K_inv = K_inv;

    Const.Radius = Radius;
    Const.w_o = w_o;
    Const.v_satellite = v_satellite;
    Const.n = n;
    Const.orbitPeriod = orbitPeriod;
    Const.I = I;
    Const.I_inv = I_inv;
    Const.w_o_io = w_o_io;
    Const.N_coils = N_coils;
    Const.A_coils = A_coils;
    Const.R_coils = R_coils;
    Const.mt = mt;
    Const.I_up_limit = I_up_limit;
    Const.W_up_limit = W_up_limit;
    Const.Jw = Jw;
    Const.A = A;
    Const.Km = Km;
    Const.Kv = Kv;
    Const.Ai = Ai;
    Const.T_friction = T_friction;
    Const.b_friction = b_friction;
    Const.c_friction = c_friction;
    Const.Rb = Rb;
    Const.lim_dz = lim_dz;
    Const.mtq_max1 = mtq_max1;
    Const.mtq_max2 = mtq_max2;
    Const.mtq_max3 = mtq_max3;
    Const.rw_max_torque = rw_max_torque;
    Const.p = p_400;
    Const.Cm = Cm;
end
