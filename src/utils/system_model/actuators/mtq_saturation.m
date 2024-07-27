% ======================================================================== %
%   Function for desaturating the MTQs in case that the induced magnetic
%   dipole exceeds the available maximum dipole provided by each
%   magnetorquer.
%
%   Inputs:
%     T_magnetic       - Torque to be provided by the MTQs
%     T_rw             - Torque to be provided by the RW
%     T_commanded      - Commanded torque
%     B_body           - Magnetic Field expressed on Body frame
%     M                - Magnetic Dipole
%     mtq_max          - Maximum Dipole provided by each Magnetorquer
%     known_rm         - Estimated constant residual magnetic dipole
%
%
%   Outputs:
%     T_magnetic	   - Torque to be provided by the MTQs
%     T_rw             - Torque to be provided by the RW
% ======================================================================== %

function [T_magnetic, T_rw] = mtq_saturation(T_magnetic, T_rw, T_commanded, B_body, M, mtq_max, known_rm)

    mtq_maxima = zeros(3, 2);
    mtq_maxima(:, 1) = mtq_max + known_rm;
    mtq_maxima(:, 2) = mtq_max - known_rm;

    M2 = (1 / norm(T_magnetic)) * M;
    Tm = T_magnetic / norm(T_magnetic);
    Tw = T_rw / norm(T_rw);
    Kma = (Tm' - (Tw' * Tm) * Tw') * T_commanded / (1 - (Tw' * Tm)^2);
    Kwa = (Tw' - (Tm' * Tw) * Tm') * T_commanded / (1 - (Tm' * Tw)^2);

    Torque_split_maxima = zeros(3, 1);

    if M(1) > mtq_maxima(1, 1) || M(2) > mtq_maxima(2, 1) || M(3) > mtq_maxima(3, 1) || M(1) < -mtq_maxima(1, 2) || M(2) < -mtq_maxima(2, 2) || M(3) < -mtq_maxima(3, 2)
        for j = 1:3
            if M(j) > 0
                Torque_split_maxima(j) = mtq_maxima(j, 1);
            elseif M(j) < 0
                Torque_split_maxima(j) = mtq_maxima(j, 2);
            end
        end
        Kma_s = min(abs(Torque_split_maxima./M2));
        Kwa_s = Kma_s * Kwa / Kma;
        Ms = Kma_s * M2;
        T_magnetic = skew(B_body)' * Ms;
        T_rw = Kwa_s .* Tw;
    else
        T_magnetic = Kma .* Tm;
        T_rw = Kwa .* Tw;
    end

end
