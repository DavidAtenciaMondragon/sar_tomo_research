function [delta_xy, delta_z, psi0, DeltaPhi] = calculateBistaticResolution(Tx_pos, Rx_pos, target_center_2d, n2, lambda, B, B_helix, beta_helix)
% CALCULATEBISTATICRESOLUTION  Estimates 3D bistatic SAR resolution for a Tx-Rx pair.
%
% Based on the k-space model derived in explicacion_resolucion.tex (Góes 2022):
%
%   delta_z  = c / (2*Wz)
%   Wz       = n2*B*cos(theta_t0) + c*Bperp*sin(psi0)*cos(psi0)/(lambda*R0*n2*cos(theta_t0))
%   delta_xy = 0.60*lambda / (pi*sin(psi0)*|cos(DeltaPhi/2)|)
%
% where psi0 is the Rx look angle toward the target centroid, DeltaPhi is the
% azimuthal separation between Tx and Rx as seen from the target centroid,
% and Bperp = B_helix*|cos(beta_helix - psi0)| is the effective tomographic baseline.
%
% INPUTS:
%   Tx_pos           - Transmitter position [3x1] (m)
%   Rx_pos           - Receiver position [3x1] (m)
%   target_center_2d - Horizontal centroid of target grid [2x1] (m)
%   n2               - Ground refractive index (dimensionless)
%   lambda           - Carrier wavelength (m)
%   B                - Chirp bandwidth (Hz)
%   B_helix          - Helix arc length sqrt(Dz^2 + Drho^2) (m)
%   beta_helix       - Helix tilt angle arctan(Dz/Drho) (rad)
%
% OUTPUTS:
%   delta_xy  - Horizontal resolution -3dB [m]  (Inf if geometry is degenerate)
%   delta_z   - Vertical resolution -3dB [m]    (Inf if geometry is degenerate)
%   psi0      - Rx look angle toward target centroid [rad]
%   DeltaPhi  - Azimuthal TX-RX separation as seen from target centroid [rad]

c = physconst('lightspeed');

%% Look angle of Rx toward the target centroid
rho_Rx = norm(Rx_pos(1:2) - target_center_2d);
z_Rx   = max(abs(Rx_pos(3)), 1.0);          % Guard against division by zero
psi0   = atan(rho_Rx / z_Rx);               % [rad]

%% Snell: transmission angle in ground medium
sin_theta_t = min(sin(psi0) / n2, 1.0);
cos_theta_t = sqrt(max(1 - sin_theta_t^2, 0));
cos_theta_t = max(cos_theta_t, 1e-9);       % Numerical guard

%% Slant range Rx → target centroid
R0 = sqrt(rho_Rx^2 + z_Rx^2);

%% Effective tomographic baseline (projection of helix perpendicular to LOS)
B_perp = B_helix * abs(cos(beta_helix - psi0));

%% Vertical resolution (invariant to DeltaPhi, depends on psi0 through Wz)
if sin(psi0) < 1e-9
    delta_z = Inf;
else
    Wz = n2 * B * cos_theta_t + ...
         (c * B_perp * sin(psi0) * cos(psi0)) / (lambda * R0 * n2 * cos_theta_t);
    if Wz > 0
        delta_z = c / (2 * Wz);
    else
        delta_z = Inf;
    end
end

%% Azimuthal TX-RX separation as seen from the target centroid
az_Tx    = atan2(Tx_pos(2) - target_center_2d(2), Tx_pos(1) - target_center_2d(1));
az_Rx    = atan2(Rx_pos(2) - target_center_2d(2), Rx_pos(1) - target_center_2d(1));
DeltaPhi = mod(az_Rx - az_Tx + pi, 2*pi) - pi;    % wrapped to (-pi, pi]

%% Horizontal resolution (J0-type PSF, -3dB criterion, factor 0.60/pi from tex §4.2)
cos_term = abs(cos(DeltaPhi / 2));
if psi0 < 1e-9 || cos_term < 1e-9
    delta_xy = Inf;
else
    delta_xy = 0.60 * lambda / (pi * sin(psi0) * cos_term);
end

end
