function [azTxRef, elTxRef, GtCurrent] = calculateTxGainsForTargets(TxPos, tg, radPattern)
% calculateTxGainsForTargets Computes Tx boresight angles and gains to all targets.

refApontRadar = [0,0,0].';
deltaTxRef = refApontRadar - TxPos;
[azTxRef, elTxRef, ~] = cart2sph(deltaTxRef(1), deltaTxRef(2), deltaTxRef(3));
azTxRef = mod(azTxRef + 2*pi, 2*pi);

numTargets = size(tg, 2);
GtCurrent = zeros(numTargets, 1);
for k = 1:numTargets
    GtCurrent(k) = getAntennaGain(tg(:,k), TxPos, azTxRef, elTxRef, radPattern);
end
end
