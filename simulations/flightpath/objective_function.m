function neg_power = objective_function(xy, tg, Tx_pos, Rx_z, radPattern, az_Tx_ref, el_Tx_ref, ...
    Gt, Pt, sigma, lambda, n1, n2)
% OBJECTIVE_FUNCTION Calcula potencia bistática total usando transmitancias de potencia
%
% Esta función calcula la potencia total recibida usando las transmitancias
% de potencia correctas (no coeficientes de Fresnel directos) para conservar
% energía en la ecuación radar bistática.

% Posición candidata del receptor
Rx_pos = [xy(1); xy(2); Rx_z];

% Inicializar potencia total
total_power = 0;

% Iterar sobre TODOS los targets en el volumen de interés
for i = 1:size(tg, 2)
    target_pos = tg(:, i);
    
    % Tramo 1: Tx → Target
    intersec_Tx_Tg = calculateRefractionPointFermat(Tx_pos, target_pos, n1, n2);
    ang_inc_1 = calculateIncidenceAngle(Tx_pos, intersec_Tx_Tg, target_pos);
    [~, ~, ~, ~, transmitancia_Tx_Tg, ~] = calculateTMcoef(ang_inc_1, n1, n2);
    R_Tx_int1 = norm(Tx_pos - intersec_Tx_Tg, 2);
    R_int1_Tg = norm(intersec_Tx_Tg - target_pos, 2);
    
    % Tramo 2: Target → Rx
    intersec_Tg_Rx = calculateRefractionPointFermat(target_pos, Rx_pos, n2, n1);
    ang_inc_2 = calculateIncidenceAngle(target_pos, intersec_Tg_Rx, Rx_pos);
    [~, ~, ~, ~, transmitancia_Tg_Rx, ~] = calculateTMcoef(ang_inc_2, n2, n1);
    R_Tg_int2 = norm(target_pos - intersec_Tg_Rx, 2);
    R_int2_Rx = norm(intersec_Tg_Rx - Rx_pos, 2);
    
    % Ganancia del receptor para este target
    delta_Tg_Rx = target_pos - Rx_pos;
    [az_Rx, el_Rx, ~] = cart2sph(delta_Tg_Rx(1), delta_Tg_Rx(2), delta_Tg_Rx(3));
    az_Rx = mod(az_Rx + 2*pi, 2*pi);
    Gr = getAntennaGain(target_pos, Rx_pos, az_Rx, el_Rx, radPattern);
    
    % Potencia bistática para este target
    % Spreading loss calculado sobre la longitud total de cada tramo refractado:
    % R_T = R_Tx_int1 + R_int1_Tg  (camino completo Tx→target)
    % R_R = R_Tg_int2 + R_int2_Rx  (camino completo target→Rx)
    R_T = R_Tx_int1 + R_int1_Tg;
    R_R = R_Tg_int2 + R_int2_Rx;
    Pr_i = (Pt * Gt(i) * Gr * sigma * ...
        transmitancia_Tx_Tg * transmitancia_Tg_Rx * lambda^2) / ...
        ((4*pi)^3 * R_T^2 * R_R^2);
    
    % Acumular potencia total
    total_power = total_power + Pr_i;
end

% Retornar negativo para maximizar (fmincon minimiza)
neg_power = -total_power;
end
