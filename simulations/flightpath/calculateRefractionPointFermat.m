function intersec_piso = calculateRefractionPointFermat(Tx_pos, target_pos, n1, n2)
% CALCULATE_REFRACTION_POINT_FERMAT  Snell refraction point on the z=0 interface.
%
% For a flat horizontal interface, Fermat's principle reduces the 2D minimization
% to a 1D root-finding problem.  The refraction point P lies on the horizontal
% line from the source projection to the target projection (provable by symmetry),
% parameterized by u in [0, D]:
%
%   f(u) = n1*u/d1 - n2*(D-u)/d2 = 0      (Snell condition)
%   d1 = sqrt(u^2 + sz^2),  d2 = sqrt((D-u)^2 + tz^2)
%   f'(u) = n1*sz^2/d1^3 + n2*tz^2/d2^3 > 0  always => unique root
%
% Newton-Raphson converges in 3-5 iterations (vs. 50-150 fminunc iterations).
%
% SYNOPSIS:
%   intersec_piso = calculateRefractionPointFermat(Tx_pos, target_pos, n1, n2)
%
% INPUT:
%   Tx_pos      - Source positions [3 x N] (x, y, z)
%   target_pos  - Target position  [3 x 1] (x, y, z)
%   n1          - Refractive index of source medium
%   n2          - Refractive index of target medium
%
% OUTPUT:
%   intersec_piso - Refraction points on z=0 plane [3 x N]

N             = size(Tx_pos, 2);
intersec_piso = zeros(3, N);

for i = 1:N
    tx  = Tx_pos(:, i);
    tgt = target_pos;

    dxy = tgt(1:2) - tx(1:2);
    D   = norm(dxy);

    sz = abs(tx(3));    % source distance to interface
    tz = abs(tgt(3));   % target distance to interface

    % Degenerate cases
    if D < 1e-12 || sz < 1e-12
        intersec_piso(:, i) = [tx(1); tx(2); 0];
        continue;
    end
    if tz < 1e-12
        intersec_piso(:, i) = [tgt(1); tgt(2); 0];
        continue;
    end

    dir = dxy / D;

    % Initial guess: linear interpolation weighted by optical depth
    u = D * (sz / n1) / (sz/n1 + tz/n2);
    u = max(0.0, min(D, u));

    % Newton-Raphson on the 1D Snell equation
    for iter = 1:12
        v  = D - u;
        d1 = sqrt(u^2 + sz^2);
        d2 = sqrt(v^2  + tz^2);
        f  =  n1 * u / d1 - n2 * v / d2;
        fp =  n1 * sz^2 / d1^3 + n2 * tz^2 / d2^3;
        du = -f / fp;
        u  =  max(0.0, min(D, u + du));
        if abs(du) < 1e-10 * D
            break;
        end
    end

    intersec_piso(:, i) = [tx(1:2) + u * dir; 0];
end

end
