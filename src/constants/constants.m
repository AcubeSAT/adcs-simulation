%% ======================================================================= %%
%   In this function all parameters regarding the physical architecture
%   of the satellite and the ADCS components are set.
%   Input: Nothing
%   Output: Struct where all values of the parameters are included
%
%   Parameters:
%       m             - Satellite mass [kg]
%       lx            - Satellite X-axis length [m]
%       ly            - Satellite Y-axis length [m]
%       lz            - Satellite Z-axis length [m]
%       N_coils       - Number of turns of the coils
%       A_coils       - Total area of the coils [m^2]
%       R_coils       - Coils resistance [Ohm]
%       mt            - Proportional parameter for calculating the current on each MTQ
%                       For torquer rods: mt = 1.66 * (length/diameter)^1.5.
%                       For air core magnetorquers: mt = 1.
%       W_up_limit    - Upper limit of absolute power of mtq's [W]
%       I_up_limit    - Upper limit of absolute current of mtq's [A]
%       Km            - Motor torque constant (RW)
%       Rb            - Armature resistance (RW)
%       Kv            - Motor velocity constant (back-EMF constant) -> 1/Km
%       b_friction    - Viscuous friction
%       c_friction    - Coulomb friction
%       Jw            - Inertia of the RW
%       A             - Parameter which indicates the percentance of RW desaturation torque
%       lim_dz        - The absolute value of the limits of the deadzone (in rpm)
%       p_400         - Atmospheric density in an altitude of 400km [kg/m^3]
%       p_500         - Atmospheric density in an altitude of 500km [kg/m^3]
%       Ix            - Xaxis inertia (simplified)
%       Iy            - Yaxis inertia (simplified)
%       Iz            - Zaxis inertia (simplified)
%       PMI           - Principal Moments of Inertia
%       PAI           - Principal Axes of Inertia
%       Cm            - Center of mass
%       I             - Inertia matrix -> PAI*PMI*PAI'
%       I_inv 	      - Inverse inertia matrix -> eye(3,3)/I;
%       Re            - Earth radius [m]
%       Rs            - Satellite altitude [m]
%       Radius        - Distance from earth center to satellite [m]
%       G             - Earth gravitational constant
%       M             - Earth mass
%       w_o           - Satellite angular velocity relative to Earth
%       v_satellite   - Satellite's velocity in orbit
%       known_rm      - Residual magnetic dipole
%       orbitPeriod   - Orbit Period
%       n             - Total simulation time
%       mtq_max       - Maximum magnetic dipole of magnetorquers
%       rw_max_torque - Maximum torque of the Reaction Wheel
%       N2D_threshold - Threshold for switching from nominal to detumbling
%       D2N_threshold - Threshold for switching from detumbling to nominal

% ======================================================================== %%

function Const = constants()

    global orbits;

    m = 3.44; % mass mentioned in STR's CDR
    %m = 4;      % mass mentioned in ADCS' CDR

    lx = 0.1;
    ly = 0.1;
    lz = 0.3405;

    N_coils = [400, 400, 800];
    A_coils = [0.0022, 0.0022, 0.001566];
    R_coils = [110, 110, 31.32];
    mt = [5, 5, 1];
    W_up_limit = [0.2273, 0.2273, 0.8017];
    I_up_limit = [0.0455, 0.0455, 0.1566];

    Km = 20e-5;
    Rb = 10;
    Ai = 1;

    Kv = 1 / Km;
    b_friction = 9.5e-9;
    c_friction = 1.9e-7;
    Jw = 1.9e-6;
    %Jw = 2.07e-6;
    A = 0.12;
    lim_dz = 300;

    p_400 = 7.55e-12;
    p_500 = 1.80e-12;

    Ix = (m / 12) * (ly^2 + lz^2);
    Iy = (m / 12) * (lx^2 + lz^2);
    Iz = (m / 12) * (lx^2 + ly^2);

    N2D_threshold = 0.08; % rad/sec                     %Threshold for switching from nominal to detumbling
    D2N_threshold = 0.035; % rad/sec                    % Threshold for switching from detumbling to nominal

    %% Old Inertia
    % PMI = diag([0.03868845951 0.03899129965 0.00696263029]); %Principal Moments of Inertia
    % PAI = [-0.89 0.46 0; 0.46 0.89 -0.01; 0 0.01 1]; %Principal Axes of Inertia
    % Cm =  [0.03228 -0.02314 0.08244]'; % Center of mass

    %% Modified Inertia
    % PMI = diag([0.03928052501	0.03948290041	0.00720041142]); %Principal Moments of Inertia
    % PAI = [-1 0.07 -0.01; 0.07 1 -0.02; -0.01 0.02 1]; %Principal Axes of Inertia
    % Cm = [0.03111 -0.02099 0.08135]';  % Center of mass

    %% THIS INERTIA MAKES ADCS HAPPY
    % PMI = diag([0.03552528444 0.03572444349 0.00626757327]);
    % PAI = [0.98, 0.2, -0.02; -0.20, 0.98, 0.00 ; 0.02, 0.01, 1.00];
    % Cm =  [0.00121 0.00057 0.00188]';

    %% CDR inertias
    PAI = [-1, 0.01, -0.04; 0.01, 1, 0.00; -0.04, 0.00, 1.00];

    % PMI mentioned in Structural's CDR
    % good case for STR , worst case for ADCS
    % PMI = diag([0.03634454760 0.03658224482 0.00626274895]);

    % PMI mentioned in ADCS' CDR
    % worst case for STR , good case for ADCS
    PMI = diag([0.04127073921, 0.041018997570, 0.00690030456]);

    Cm = [0.00415, 0.00116, 0.0016]';


    for j = 1:3
        PAI(:, j) = PAI(:, j) / norm(PAI(:, j));
    end
    I = PAI * PMI * PAI';
    I_inv = eye(3, 3) / I;

    Re = 6371.2e3;
    Rs = 500e3;
    Radius = Re + Rs;
    G = 6.67428e-11;
    M = 5.972e24;
    w_o = sqrt(G*M/Radius^3);
    v_satellite = sqrt(G*M/Radius);
    w_o_io = [0, w_o, 0]';

    known_rm = [0.048, 0.051, 0.047];

    orbitPeriod = (2 * pi) / (w_o);
    n = orbitPeriod * orbits;
    mtq_max = [0.2, 0.2, 0.2];
    rw_max_torque = 1e-4;

    %%  Passing the values of the parameters in a struct.

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
    Const.b_friction = b_friction;
    Const.c_friction = c_friction;
    Const.Rb = Rb;
    Const.lim_dz = lim_dz;
    Const.mtq_max = mtq_max;
    Const.rw_max_torque = rw_max_torque;
    Const.p = p_500;
    Const.Cm = Cm;
    Const.known_rm = known_rm;
    Const.N2D_threshold= N2D_threshold;
    Const.D2N_threshold= D2N_threshold;
end
