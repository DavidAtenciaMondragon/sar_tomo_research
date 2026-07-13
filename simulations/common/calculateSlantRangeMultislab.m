function [Rseg, nSeg, theta_inc, L] = calculateSlantRangeMultislab(Px,Py,Pz,x0,y0,z0,n_layers,d_layers,threshold)
% CALCULATESLANTRANGEMULTISLAB  Geometria de rayo (Fermat/Snell) para un
% sistema de N capas de suelo planas y horizontales, con el target embebido
% dentro de una de ellas.
%
% El angulo de incidencia (medido desde la normal) y los angulos de
% refraccion dentro de cada capa quedan ligados por el invariante de Snell
% generalizado n(i)*sin(theta(i)) = p = cte (ver models/multiple_slab_tm.md,
% Seccion 2). Dado p, la traza horizontal de cada tramo del rayo es
% h(i)*tan(theta(i)), con h(i) la altura/espesor vertical de ese tramo.
% El parametro de rayo p se obtiene por biseccion de modo que la suma de las
% trazas horizontales de todos los tramos sea igual a la distancia horizontal
% entre el radar y el target.
%
% SYNOPSIS:
%   [Rseg, nSeg, theta_inc] = calculateSlantRangeMultislab(Px,Py,Pz,x0,y0,z0,n_layers,d_layers,threshold)
%
% INPUT:
%   Px,Py,Pz  - Posicion del radar (trayectoria) [1 x Np]. Pz > 0 (aire).
%   x0,y0,z0  - Posicion del target [m]. z0 < 0 (profundidad bajo superficie).
%   n_layers  - Indices de refraccion de las capas de suelo, de superficie
%               hacia abajo [1 x M].
%   d_layers  - Espesores de las capas de suelo [m] [1 x M]. Para la capa que
%               contiene al target solo se usa la fraccion de espesor hasta
%               el target; el resto del valor declarado se ignora.
%   threshold - Tolerancia de convergencia de la biseccion sobre la traza
%               horizontal [m].
%
% OUTPUT:
%   Rseg      - Rangos oblicuos de cada tramo [(L+1) x Np]. Fila 1 = tramo en
%               aire (radar -> superficie). Filas 2..L+1 = tramos dentro de
%               las capas de suelo 1..L, donde L es la capa que contiene al
%               target (la fila L+1 es el tramo parcial hasta el target).
%   nSeg      - Indices de refraccion de cada tramo [(L+1) x 1], con
%               nSeg(1) = 1 (aire).
%   theta_inc - Angulo de incidencia en aire, respecto a la normal [rad]
%               [1 x Np].
%   L         - Indice (1-based) de la capa de suelo que contiene al target.

Np = numel(Px);

% --- Localizar la capa que contiene al target ---------------------------
depthTarget = -z0;
cumDepth    = cumsum(d_layers(:)');
L = find(cumDepth >= depthTarget, 1, 'first');
if isempty(L)
    % La ultima capa es semi-infinita: el target puede estar arbitrariamente
    % profundo dentro de ella, independientemente del d declarado en el JSON.
    L = numel(d_layers);
end

if L == 1
    hLast = depthTarget;
else
    hLast = depthTarget - cumDepth(L-1);
end

% --- Alturas/espesores y refractividad de cada tramo ---------------------
% Fila 1: aire; filas 2..L: capas completas por encima del target;
% fila L+1: tramo parcial dentro de la capa del target.
h    = zeros(L+1, Np);
nSeg = zeros(L+1, 1);

h(1,:)  = Pz(:)';
nSeg(1) = 1;

for i = 1:(L-1)
    h(i+1,:) = d_layers(i);
    nSeg(i+1) = n_layers(i);
end

h(L+1,:)  = hLast;
nSeg(L+1) = n_layers(L);

% --- Distancia horizontal radar-target ------------------------------------
u = sqrt((Px(:)'-x0).^2 + (Py(:)'-y0).^2);

% --- Biseccion vectorizada del parametro de rayo p -------------------------
% p = n(i)*sin(theta(i)) = cte, con 0 <= p < min(nSeg) para que todos los
% angulos sean reales (sin reflexion total interna).
lo = zeros(1,Np);
hi = (min(nSeg) - 1e-12) * ones(1,Np);

maxIter = 100;
for iter = 1:maxIter
    p = (lo+hi)/2;

    Up = zeros(1,Np);
    for i = 1:(L+1)
        Up = Up + h(i,:) .* tan(asin(p./nSeg(i)));
    end

    mask = Up < u;
    lo(mask)  = p(mask);
    hi(~mask) = p(~mask);

    if max(hi-lo) < threshold
        break;
    end
end
p = (lo+hi)/2;

% --- Rangos oblicuos y angulo de incidencia --------------------------------
Rseg = zeros(L+1, Np);
for i = 1:(L+1)
    theta_i = asin(p./nSeg(i));
    Rseg(i,:) = h(i,:) ./ cos(theta_i);
end

theta_inc = asin(p./nSeg(1));

end
