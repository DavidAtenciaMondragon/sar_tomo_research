function strRadar = calculateCoverage(strRadar)

strRadar.AperturaElevRad = strRadar.AperturaElev*pi/180; % rad
strRadar.InclElevRad  = strRadar.InclElev*pi/180; % rad

rangeMin = tan(strRadar.InclElevRad)*strRadar.posInicial(3);
rangeMax = tan(strRadar.AperturaElevRad + strRadar.InclElevRad)*strRadar.posInicial(3);

strRadar.AperturaAzimutRad = strRadar.AperturaAzimut*pi/180; % rad
pontoMax = tan(strRadar.AperturaAzimutRad/2)*rangeMax;

% --

figure
subplot(121)
hold on
plot([strRadar.posInicial(1), rangeMin(1)], [strRadar.posInicial(3),0],'--b')
plot([strRadar.posInicial(1), rangeMax(1)], [strRadar.posInicial(3),0],'--b')
hold off
grid on
xlabel('x')
ylabel('z')
title("Range min e max (top view)")
title('Radar beamwidth')
grid minor

subplot(122)
hold on
plot([0, rangeMax],[0, pontoMax],'--b')
plot([0, rangeMax],[0, -pontoMax],'--b')
xlabel("x")
ylabel("y")
title("Abertura em azimute")
grid minor
hold off

end