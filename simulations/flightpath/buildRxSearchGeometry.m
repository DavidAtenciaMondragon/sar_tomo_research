function searchGeo = buildRxSearchGeometry(TxPos, RxZ, tg, anguloBrewster)
% buildRxSearchGeometry Precomputes geometry used for Rx candidate initialization.

searchGeo = struct();

searchGeo.target_center = [mean(tg(1,:)), mean(tg(2,:))];
txToTargets = searchGeo.target_center - [TxPos(1), TxPos(2)];
searchGeo.targets_distance = norm(txToTargets);

if searchGeo.targets_distance > eps
    searchGeo.opposite_dir = txToTargets / searchGeo.targets_distance;
    searchGeo.perp_dir = [-searchGeo.opposite_dir(2), searchGeo.opposite_dir(1)];
else
    searchGeo.opposite_dir = [1, 0];
    searchGeo.perp_dir = [0, 1];
end

txDirection = [TxPos(1), TxPos(2)] - searchGeo.target_center;
searchGeo.azimut_opuesto = atan2(txDirection(2), txDirection(1)) + pi;

targetDepth = abs(mean(tg(3,:)));
alturaDiferencia = RxZ + targetDepth;
searchGeo.distancia_brewster = alturaDiferencia / tan(anguloBrewster);
end
