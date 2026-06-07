function intersec_piso = calculateRefractionPointFermat(Tx_pos, target_pos, n1, n2)
% CALCULATE_REFRACTION_POINT_FERMAT Calculates the intersection point on a flat surface
% using Fermat's principle
%
% SYNOPSIS:
%   intersec_piso = calculateRefractionPointFermat(Tx_pos, target_pos, n1, n2)
%
% INPUT:
%   Tx_pos      - Transmitter positions [3 x N] (x, y, z)
%   target_pos  - Target position [3 x 1] (x , y, z)  
%   n1          - Refractive index of medium 1 (air)
%   n2          - Refractive index of medium 2 (subsurface)
%
% OUTPUT:
%   intersec_piso - Intersection points on the plane z=0 [3 x N] (x, y, 0)
%
% DESCRIPTION:
%   This function implements Fermat's principle to find the intersection points
%   on a flat surface (z = 0) that minimize the total optical path length between
%   each transmitter position and the target, considering refraction according
%   to Snell's law.

% NNumber of transmitter positions
N = size(Tx_pos, 2);

% Inicializar matriz de salida
intersec_piso = zeros(3, N);

% Configuración para la optimización
options = optimoptions('fminunc', 'Display', 'off', 'Algorithm', 'quasi-newton');

for i = 1:N
    % Posición actual del transmisor
    tx_current = Tx_pos(:, i);
    
    % Punto inicial para la optimización (proyección directa en el plano)
    % Usamos interpolación lineal como estimación inicial
    if tx_current(3) ~= target_pos(3)
        t_intersec = -tx_current(3) / (target_pos(3) - tx_current(3));
        x0_inicial = tx_current(1) + t_intersec * (target_pos(1) - tx_current(1));
        y0_inicial = tx_current(2) + t_intersec * (target_pos(2) - tx_current(2));
    else
        % Si ambos puntos están en el mismo plano z, usar punto medio
        x0_inicial = (tx_current(1) + target_pos(1)) / 2;
        y0_inicial = (tx_current(2) + target_pos(2)) / 2;
    end
    
    xy0 = [x0_inicial; y0_inicial];
    
    % Definir función objetivo: tiempo óptico total
    objective_function = @(xy) tiempo_optico_total(xy, tx_current, target_pos, n1, n2);
    
    % Optimizar para encontrar el punto de Fermat
    [xy_optimal, ~] = fminunc(objective_function, xy0, options);
    
    % Guardar resultado
    intersec_piso(:, i) = [xy_optimal(1); xy_optimal(2); 0];
end

end

function T = tiempo_optico_total(xy, tx_pos, target_pos, n1, n2)
% TIEMPO_OPTICO_TOTAL Calcula el tiempo óptico total para un punto de intersección
%
% INPUT:
%   xy         - Coordenadas [x, y] del punto de intersección en z=0
%   tx_pos     - Posición del transmisor [3x1]
%   target_pos - Posición del target [3x1]
%   n1         - Índice de refracción del medio 1
%   n2         - Índice de refracción del medio 2

% Punto de intersección en la superficie
intersec_point = [xy(1); xy(2); 0];

% Distancia en el medio 1 (desde transmisor hasta intersección)
d1 = norm(intersec_point - tx_pos);

% Distancia en el medio 2 (desde intersección hasta target)
d2 = norm(target_pos - intersec_point);

% Tiempo óptico total (proporcional a n*distancia)
T = n1 * d1 + n2 * d2;

end