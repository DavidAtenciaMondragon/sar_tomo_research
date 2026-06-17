function cost = objective_function_multislab(xy, tg, Tx_pos, Rx_z, radPattern, Gt, Pt, sigma, lambda, ...
    n_layers, d_layers, f, threshold, alpha_res, B, B_helix, beta_helix, target_center_2d, ...
    n_eff, R_T_pre, T1_pre)
% OBJECTIVE_FUNCTION_MULTISLAB  Bistatic cost for an N-layer ground model.
%
% Extends objective_function.m to multi-slab media.  The single-interface
% Fermat refraction point and Fresnel TM coefficient are replaced by:
%   - calculateSlantRangeMultislab : bisection ray solver across N layers
%   - multislabTMcoef              : ABCD-cascade end-to-end TM transmittance
%
% Radar equation (per target i):
%
%   Pr_i = Pt * Gt_i * Gr_i * sigma * T1_i * T2_i * lambda^2
%          -------------------------------------------------------
%          (4*pi)^3 * R_T_i^2 * R_R_i^2
%
%   where R_T_i / R_R_i are the total physical path lengths (sum of all
%   ray segments in air + soil layers) for the Tx→target and target→Rx
%   legs respectively, and T1_i / T2_i are the corresponding end-to-end
%   TM power transmittances from air to the target layer.
%
% Cost function (dimensionless, fmincon minimizes):
%   J = (1-alpha_res)*(-log10(Pr_total)) + alpha_res*log10(delta_xy*delta_z)
%
% INPUTS:
%   xy               - Candidate Rx horizontal position [x; y] (m)
%   tg               - Target grid [3 x Ntg] (m)
%   Tx_pos           - Transmitter position [3x1] (m)
%   Rx_z             - Receiver altitude (m)
%   radPattern       - Antenna radiation pattern [181x181]
%   Gt               - Pre-computed Tx gain toward each target [Ntg x 1]
%   Pt               - Transmitted power (W)
%   sigma            - Target RCS (m^2)
%   lambda           - Carrier wavelength (m)
%   n_layers         - Soil layer refractive indices [1 x M]
%   d_layers         - Soil layer thicknesses [m] [1 x M]
%   f                - Carrier frequency (Hz)
%   threshold        - Bisection tolerance for ray solver (m)
%   alpha_res        - Resolution weight in [0,1]
%   B                - Chirp bandwidth (Hz)
%   B_helix          - Helix arc length (m)
%   beta_helix       - Helix tilt angle (rad)
%   target_center_2d - Horizontal centroid of target grid [2x1] (m)
%   n_eff            - Effective refractive index of target layer (for resolution)
%   R_T_pre          - Pre-computed Tx→Target path lengths [1 x Ntg] (m)
%   T1_pre           - Pre-computed Tx→Target TM transmittance [1 x Ntg]

Rx_pos      = [xy(1); xy(2); Rx_z];
total_power = 0;

for i = 1:size(tg, 2)
    target_pos = tg(:, i);

    % Tx→Target (precomputed, constant for this Tx position)
    R_T = R_T_pre(i);
    T1  = T1_pre(i);

    % Target→Rx: solve refracted ray with the same multislab geometry
    [Rseg_R, ~, theta_inc_Rx, L] = calculateSlantRangeMultislab( ...
        Rx_pos(1), Rx_pos(2), Rx_pos(3), ...
        target_pos(1), target_pos(2), target_pos(3), ...
        n_layers, d_layers, threshold);

    R_R = sum(Rseg_R);

    n_slabs_Rx = n_layers(1:L-1);
    d_slabs_Rx = d_layers(1:L-1);
    [~, ~, ~, T2] = multislabTMcoef(f, theta_inc_Rx, 1.0, n_slabs_Rx, d_slabs_Rx, n_layers(L));

    % Rx antenna gain toward this target
    delta_Tg_Rx = target_pos - Rx_pos;
    [az_Rx, el_Rx, ~] = cart2sph(delta_Tg_Rx(1), delta_Tg_Rx(2), delta_Tg_Rx(3));
    az_Rx = mod(az_Rx + 2*pi, 2*pi);
    Gr    = getAntennaGain(target_pos, Rx_pos, az_Rx, el_Rx, radPattern);

    % Bistatic radar equation
    Pr_i = (Pt * Gt(i) * Gr * sigma * T1 * T2 * lambda^2) / ...
           ((4*pi)^3 * R_T^2 * R_R^2);

    total_power = total_power + Pr_i;
end

%% Power term (log-scale, dimensionless)
if total_power > 0 && isfinite(total_power)
    log_power = log10(total_power);
else
    log_power = -300;
end

%% Resolution term (log-scale, dimensionless)
if alpha_res > 0
    [delta_xy, delta_z] = calculateBistaticResolution(Tx_pos, Rx_pos, target_center_2d, ...
        n_eff, lambda, B, B_helix, beta_helix);

    if isfinite(delta_xy) && isfinite(delta_z) && delta_xy > 0 && delta_z > 0
        log_res = log10(delta_xy * delta_z);
    else
        log_res = 10;
    end
else
    log_res = 0;
end

%% Combined cost (fmincon minimizes J)
cost = -(1 - alpha_res) * log_power + alpha_res * log_res;

end
