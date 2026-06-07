function [R_T_pre, T1_pre] = precomputeTxData(Tx_pos, tg, n1, n2)
% PRECOMPUTETXDATA  Pre-computes Tx→Target refraction data for a fixed Tx position.
%
% The Tx→Target leg is independent of the Rx position being optimized by fmincon.
% Pre-computing it once per Tx position avoids recomputing the Fermat refraction
% point for every objective_function evaluation — eliminating ~50% of Fermat calls.
%
% INPUTS:
%   Tx_pos  - Transmitter position [3x1] (m)
%   tg      - Target grid [3 x Ntg] (m)
%   n1      - Air refractive index
%   n2      - Ground refractive index
%
% OUTPUTS:
%   R_T_pre - Total Tx→interface + interface→Target path lengths [1 x Ntg] (m)
%   T1_pre  - TM power transmittance Tx→Target at each refraction point [1 x Ntg]

Ntg    = size(tg, 2);
R_T_pre = zeros(1, Ntg);
T1_pre  = zeros(1, Ntg);

for i = 1:Ntg
    target_pos = tg(:, i);
    intersec   = calculateRefractionPointFermat(Tx_pos, target_pos, n1, n2);
    ang_inc    = calculateIncidenceAngle(Tx_pos, intersec, target_pos);
    [~,~,~,~, T1, ~] = calculateTMcoef(ang_inc, n1, n2);
    R_T_pre(i) = norm(Tx_pos - intersec) + norm(intersec - target_pos);
    T1_pre(i)  = T1;
end

end
