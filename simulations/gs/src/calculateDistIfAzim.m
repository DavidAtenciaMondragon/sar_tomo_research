function dist = calculateDistIfAzim(posRadar,strTarget, strRadar, strSystem)

dist = [];

switch strRadar.tipoCobertura
    case 'stripmat'
        Yaw = strRadar.Yaw;
    case 'spotlight'
        % Mover el centro de referencia a la posicion del radar
        posRefAux = strSystem.RefGeografico - posRadar;
        % Calcular azimuth y elevacion
        [Yaw, Roll, ~] = cart2sph(posRefAux(1),posRefAux(2),posRefAux(3));
        
    otherwise
        error('Tipo de cobertura não reconhecido')
end

Yaw = Yaw*180/pi;

% Mover target a la posicion del radar

posTargetAux = strTarget.pos - posRadar;

% Calcular azimuth y elevacion

[azTarget, elTarget, ~] = cart2sph(posTargetAux(1),posTargetAux(2),posTargetAux(3));

anguloMinBusca = Yaw - strRadar.AperturaAzimut/2;
anguloMaxBusca = Yaw + strRadar.AperturaAzimut/2;

if angulo_en_rango(azTarget*180/pi,anguloMinBusca,anguloMaxBusca)
    
    % Calcula a distância entre o radar e o alvo
    dist = norm(posRadar(:) - strTarget.pos);
    
end

end