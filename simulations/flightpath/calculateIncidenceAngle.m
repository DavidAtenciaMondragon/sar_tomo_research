function angulo_incidencia = calculateIncidenceAngle(tx_pos, intersec_point, target_pos)
% CALCULATE_INCIDENCE_ANGLE Calculates the incidence angle of the ray with the normal
%
% SYNOPSIS:
%   angulo_incidencia = calculateIncidenceAngle(tx_pos, intersec_point, target_pos)
%
% INPUT:
%   tx_pos         - Transmitter position [3 x N] or [3 x 1]
%   intersec_point - Intersection point on the surface [3 x N] or [3 x 1]
%   target_pos     - Target position [3 x 1]
%
% OUTPUT:
%   angulo_incidencia - Incidence angle in radians [1 x N]
%
% DESCRIPTION:
%   This function calculates the incidence angle of the ray from the
%   transmitter to the intersection point with respect to the normal of
%   the flat surface (which points upwards: [0, 0, 1])

% Surface normal vector (points upwards)
normal_superficie = [0; 0; 1];

% Determine the number of positions
if size(tx_pos, 2) > 1
    N = size(tx_pos, 2);
else
    N = 1;
end

% Inicializar salida
angulo_incidencia = zeros(1, N);

for i = 1:N
    % Obtener posiciones para esta iteración
    if size(tx_pos, 2) > 1
        tx_current = tx_pos(:, i);
        intersec_current = intersec_point(:, i);
    else
        tx_current = tx_pos;
        intersec_current = intersec_point;
    end
    
    % Vector del rayo incidente (desde transmisor hacia intersección)
    rayo_incidente = intersec_current - tx_current;
    
    % Normalizar el vector
    rayo_incidente_norm = rayo_incidente / norm(rayo_incidente);
    
    % Calcular el ángulo con la normal usando producto punto
    % cos(θ) = |a·b| / (|a||b|)
    % Como ambos vectores están normalizados: cos(θ) = |a·b|
    cos_angulo = abs(dot(rayo_incidente_norm, normal_superficie));
    
    % El ángulo de incidencia es el ángulo con respecto a la normal
    angulo_incidencia(i) = acos(cos_angulo);
end

end