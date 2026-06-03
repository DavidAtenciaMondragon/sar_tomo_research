function [strRadarTx, strRadarRx] = calculateCoverageBistatic(strRadarTx,strRadarRx)

switch strRadarTx.tipoCobertura
    case 'stripmat'
        disp("------------------------------------------------------------------------------------")
        disp("Cobertura em modo STRIPMAT, angulo de apontamento do radar configurado por parametro")
        disp("------------------------------------------------------------------------------------")
    case 'spotlight'
        disp("----------------------------------------------------------------------------------------------------------------------------")
        disp("Cobertura em modo SPOTLIGHT, angulo de apontamento do radar calculado internamente (sempre apontando num ponto de referencia")
        disp("----------------------------------------------------------------------------------------------------------------------------")
    otherwise
        error('Tipo de cobertura não reconhecido')
end

posInicialCartTx = convertCordinatesToCartesian(strRadarTx.posInicial,strRadarTx.tipoTrajetoria);
posInicialCartRx = convertCordinatesToCartesian(strRadarRx.posInicial,strRadarRx.tipoTrajetoria);

% Transmisor
strRadarTx.AperturaElevRad = strRadarTx.AperturaElev*pi/180; % rad
strRadarTx.InclElevRad  = strRadarTx.InclElev*pi/180; % rad

rangeMinTx = tan(strRadarTx.InclElevRad)*posInicialCartTx(3)*cosd(strRadarTx.Yaw);
rangeMaxTx = tan(strRadarTx.AperturaElevRad + strRadarTx.InclElevRad)*posInicialCartTx(3)*cosd(strRadarTx.Yaw);

strRadarTx.AperturaAzimutRad = strRadarTx.AperturaAzimut*pi/180; % rad
pontoMaxTx = tan(strRadarTx.AperturaAzimutRad/2)*rangeMaxTx*cosd(strRadarTx.Yaw);

% Receptor
strRadarRx.AperturaElevRad = strRadarRx.AperturaElev*pi/180; % rad
strRadarRx.InclElevRad     = strRadarRx.InclElev*pi/180; % rad

rangeMinRx = tan(strRadarRx.InclElevRad)*posInicialCartRx(3)*cosd(strRadarRx.Yaw);
rangeMaxRx = tan(strRadarRx.AperturaElevRad + strRadarRx.InclElevRad)*posInicialCartRx(3)*cosd(strRadarRx.Yaw);

strRadarRx.AperturaAzimutRad = strRadarRx.AperturaAzimut*pi/180; % rad
pontoMaxRx = tan(strRadarRx.AperturaAzimutRad/2)*rangeMaxRx*cosd(strRadarRx.Yaw);

% --

figure

subplot(121)
hold on
plot([posInicialCartTx(1), rangeMinTx(1) + posInicialCartTx(1)], [posInicialCartTx(3),0],'--b')
plot([posInicialCartTx(1), rangeMaxTx(1) + posInicialCartTx(1)], [posInicialCartTx(3),0],'--b')
hold off
grid on
xlabel('x')
ylabel('z')
title("Range min e max (front view)")
title('Tx')
grid minor

subplot(122)
hold on
plot([posInicialCartRx(1), rangeMinRx(1) + posInicialCartRx(1)], [posInicialCartRx(3),0],'--b')
plot([posInicialCartRx(1), rangeMaxRx(1) + posInicialCartRx(1)], [posInicialCartRx(3),0],'--b')
hold off
grid on
xlabel('x')
ylabel('z')
title("Range min e max (front view)")
title('Rx')
grid minor


figure

subplot(121)
hold on
plot([0, rangeMaxTx],[0, pontoMaxTx],'--b')
plot([0, rangeMaxTx],[0, -pontoMaxTx],'--b')
xlabel("x")
ylabel("y")
title("Abertura em azimute - Tx")
grid minor
hold off

subplot(122)
hold on
plot([0, rangeMaxRx],[0, pontoMaxRx],'--b')
plot([0, rangeMaxRx],[0, -pontoMaxRx],'--b')
xlabel("x")
ylabel("y")
title("Abertura em azimute - Rx")
grid minor
hold off

figure

hold on
plot([posInicialCartTx(1), rangeMinTx(1) + posInicialCartTx(1)], [posInicialCartTx(3),0],'--b')
plot([posInicialCartTx(1), rangeMaxTx(1) + posInicialCartTx(1)], [posInicialCartTx(3),0],'--b')
plot([posInicialCartRx(1), rangeMinRx(1) + posInicialCartRx(1)], [posInicialCartRx(3),0],'--r')
plot([posInicialCartRx(1), rangeMaxRx(1) + posInicialCartRx(1)], [posInicialCartRx(3),0],'--r')
hold off
grid on
xlabel('x')
ylabel('z')
title("Cobertura (front view)")
grid minor

end