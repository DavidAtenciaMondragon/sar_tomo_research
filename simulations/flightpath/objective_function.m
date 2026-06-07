function cost = objective_function(xy, tg, Tx_pos, Rx_z, radPattern, Gt, Pt, sigma, lambda, n1, n2, ...
    alpha_res, B, B_helix, beta_helix, target_center_2d)
% OBJECTIVE_FUNCTION  Combined bistatic cost: power vs. 3D resolution tradeoff.
%
% Minimizes a dimensionless log-normalized cost that simultaneously maximizes
% received bistatic power (Brewster effect) and minimizes 3D SAR resolution:
%
%   J = (1 - alpha_res) * (-log10(Pr)) + alpha_res * log10(delta_xy * delta_z)
%
% Both terms are dimensionless (log10 of quantities normalized to 1 W and 1 m²),
% making the weighted sum physically meaningful across different scales.
%
% alpha_res = 0  → pure power maximization (original behavior)
% alpha_res = 1  → pure resolution minimization
% alpha_res ∈ (0,1) → tradeoff; recommended starting value: 0.3
%
% INPUTS:
%   xy               - Candidate Rx horizontal position [x; y] (m)
%   tg               - Target grid positions [3 x Ntg] (m)
%   Tx_pos           - Transmitter position [3x1] (m)
%   Rx_z             - Receiver altitude (m)
%   radPattern       - Antenna radiation pattern matrix [181x181]
%   Gt               - Pre-computed Tx gain toward each target [Ntg x 1]
%   Pt               - Transmitted power (W)
%   sigma            - Target RCS (m²)
%   lambda           - Carrier wavelength (m)
%   n1               - Air refractive index (= 1)
%   n2               - Ground refractive index
%   alpha_res        - Resolution weight ∈ [0, 1]
%   B                - Chirp bandwidth (Hz)  [used for delta_z]
%   B_helix          - Helix arc length (m)  [used for delta_z]
%   beta_helix       - Helix tilt angle (rad)[used for delta_z]
%   target_center_2d - Horizontal centroid of target grid [2x1] (m)

Rx_pos      = [xy(1); xy(2); Rx_z];
total_power = 0;

for i = 1:size(tg, 2)
    target_pos = tg(:, i);

    % Tramo 1: Tx → Target (refracción en interfaz, ida)
    intersec_Tx_Tg    = calculateRefractionPointFermat(Tx_pos, target_pos, n1, n2);
    ang_inc_1         = calculateIncidenceAngle(Tx_pos, intersec_Tx_Tg, target_pos);
    [~,~,~,~, T1, ~]  = calculateTMcoef(ang_inc_1, n1, n2);
    R_T = norm(Tx_pos - intersec_Tx_Tg) + norm(intersec_Tx_Tg - target_pos);

    % Tramo 2: Target → Rx (refracción en interfaz, vuelta)
    intersec_Tg_Rx    = calculateRefractionPointFermat(target_pos, Rx_pos, n2, n1);
    ang_inc_2         = calculateIncidenceAngle(target_pos, intersec_Tg_Rx, Rx_pos);
    [~,~,~,~, T2, ~]  = calculateTMcoef(ang_inc_2, n2, n1);
    R_R = norm(target_pos - intersec_Tg_Rx) + norm(intersec_Tg_Rx - Rx_pos);

    % Ganancia del receptor hacia este target
    delta_Tg_Rx = target_pos - Rx_pos;
    [az_Rx, el_Rx, ~] = cart2sph(delta_Tg_Rx(1), delta_Tg_Rx(2), delta_Tg_Rx(3));
    az_Rx = mod(az_Rx + 2*pi, 2*pi);
    Gr    = getAntennaGain(target_pos, Rx_pos, az_Rx, el_Rx, radPattern);

    % Potencia bistática — ecuación de radar con caminos refractados totales
    % Spreading loss sobre la longitud total de cada tramo: R_T = R1+R2, R_R = R3+R4
    Pr_i = (Pt * Gt(i) * Gr * sigma * T1 * T2 * lambda^2) / ...
           ((4*pi)^3 * R_T^2 * R_R^2);

    total_power = total_power + Pr_i;
end

%% Término de potencia en log10 (escala adimensional)
if total_power > 0 && isfinite(total_power)
    log_power = log10(total_power);   % típico: -8 a -12
else
    log_power = -300;                  % Penalización extrema
end

%% Término de resolución en log10 (escala adimensional)
if alpha_res > 0
    [delta_xy, delta_z] = calculateBistaticResolution(Tx_pos, Rx_pos, target_center_2d, ...
        n2, lambda, B, B_helix, beta_helix);

    if isfinite(delta_xy) && isfinite(delta_z) && delta_xy > 0 && delta_z > 0
        log_res = log10(delta_xy * delta_z);  % típico: -2 a -4
    else
        log_res = 10;    % Penalización: geometría degenerada (Δφ→180° o ψ₀→0°)
    end
else
    log_res = 0;         % No se evalúa si alpha=0
end

%% Función de costo combinada (fmincon minimiza J)
%   -log10(Pr) se minimiza al aumentar Pr              ✓
%   log10(δxy·δz) se minimiza al reducir la resolución ✓
cost = -(1 - alpha_res) * log_power + alpha_res * log_res;

end
