function [R_T_pre, T1_pre] = precomputeTxDataMultislab(Tx_pos, tg, n_layers, d_layers, f, threshold)
% PRECOMPUTETXDATAMULTISLAB  Pre-computes Tx→Target path data for a multislab medium.
%
% Replaces precomputeTxData for N-layer ground models. For each target the
% refracted ray is solved by bisection over the Snell invariant p (via
% calculateSlantRangeMultislab), and the end-to-end TM power transmittance
% through the full layer stack is obtained from the ABCD cascade (via
% multislabTMcoef).
%
% The Tx→Target leg is constant for a fixed Tx position, so pre-computing it
% once eliminates the majority of expensive ray-solver calls from the inner
% fmincon loop.
%
% INPUTS:
%   Tx_pos    - Transmitter position [3x1] (m).  Tx_pos(3) > 0 (in air).
%   tg        - Target grid [3 x Ntg] (m).  tg(3,:) < 0 (below surface).
%   n_layers  - Refractive indices of soil layers, surface→down [1 x M]
%   d_layers  - Thicknesses of soil layers [m] [1 x M]
%   f         - Carrier frequency [Hz]
%   threshold - Bisection convergence tolerance [m] (e.g. 1e-6)
%
% OUTPUTS:
%   R_T_pre - Total physical Tx→Target path (sum of all ray segments) [1 x Ntg] (m)
%   T1_pre  - End-to-end TM power transmittance from air to target layer [1 x Ntg]

Ntg     = size(tg, 2);
R_T_pre = zeros(1, Ntg);
T1_pre  = zeros(1, Ntg);

for i = 1:Ntg
    target_pos = tg(:, i);

    [Rseg, ~, theta_inc, L] = calculateSlantRangeMultislab( ...
        Tx_pos(1), Tx_pos(2), Tx_pos(3), ...
        target_pos(1), target_pos(2), target_pos(3), ...
        n_layers, d_layers, threshold);

    R_T_pre(i) = sum(Rseg);

    % Intermediate slabs between air and the target layer (empty if L==1,
    % in which case multislabTMcoef reduces to single-interface Fresnel).
    n_slabs = n_layers(1:L-1);
    d_slabs = d_layers(1:L-1);
    [~, ~, ~, T1] = multislabTMcoef(f, theta_inc, 1.0, n_slabs, d_slabs, n_layers(L));
    T1_pre(i) = T1;
end

end
