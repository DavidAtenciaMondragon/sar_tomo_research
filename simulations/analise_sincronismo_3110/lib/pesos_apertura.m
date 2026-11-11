function w = pesos_apertura(pos)
% PESOS_APERTURA  Ponderacion en azimut a lo largo de la trayectoria (vuelo circular).
%   Arco parcial: ventana Hann sobre el azimut del radar respecto al centro.
%   Vuelta completa: uniforme (la apertura circular no tiene bordes).
th = unwrap(atan2(pos(:,2), pos(:,1)));
span = max(th) - min(th);
if span < eps || span > 0.98*2*pi, w = ones(size(th)); return; end
w = 0.5 - 0.5*cos(2*pi*(th - min(th))/span);
w = w + 1e-3;                                  % evita peso nulo con 1-2 ventanas
end
