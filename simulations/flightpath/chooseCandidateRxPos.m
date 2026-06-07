function x0 = chooseCandidateRxPos(startIdx, targetCenter, targetsDistance, oppositeDir, perpDir, distanciaBrewster, azimutOpuesto)
% chooseCandidateRxPos Selects one candidate XY start point for Rx optimization.
%   Keeps the original heuristic:
%   1) opposite to targets direction,
%   2) perpendicular direction,
%   3) opposite perpendicular direction,
%   4-5) around Brewster-guided point with random jitter.

if startIdx == 1
    % Candidate 1: opposite to target cluster (bistatic angle ~180 deg)
    x0 = targetCenter + oppositeDir * targetsDistance * 0.8 + randn(1,2);
elseif startIdx == 2
    % Candidate 2: first perpendicular direction (bistatic angle ~90 deg)
    x0 = targetCenter + perpDir * targetsDistance * 0.8 + randn(1,2);
elseif startIdx == 3
    % Candidate 3: opposite perpendicular direction
    x0 = targetCenter - perpDir * targetsDistance * 0.8 + randn(1,2);
else
    % Candidates 4-5: Brewster-based region with uniform random perturbation
    puntoBrewster = targetCenter + distanciaBrewster * [cos(azimutOpuesto), sin(azimutOpuesto)];
    noiseFactor = 0.2 * distanciaBrewster;
    x0 = puntoBrewster + noiseFactor * (2 * rand(1,2) - 1);
end
end
