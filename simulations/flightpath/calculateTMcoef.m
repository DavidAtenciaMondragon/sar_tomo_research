function [coef_Transm, coef_Reflex, angulo_brewster, distancia_brewster, transmitancia_potencia, reflectancia_potencia] = calculateTMcoef(angulo_incidencia, n1, n2)
% CALCULATE_TMCOEF Calculates transmission and reflection coefficients for TM waves (p-polarization)
%
% SYNOPSIS:
%   [coef_Transm, coef_Reflex, angulo_brewster, distancia_brewster, transmitancia_potencia, reflectancia_potencia] = calculateTMcoef(angulo_incidencia, n1, n2)
%
% INPUT:
%   angulo_incidencia - Incidence angle in radians [1 x N]
%   n1                - Refractive index of medium 1 (incident)
%   n2                - Refractive index of medium 2 (transmitted)
%
% OUTPUT:
%   coef_Transm        - TM transmission coefficients (electric field) [1 x N]
%   coef_Reflex        - TM reflection coefficients (electric field) [1 x N]
%   angulo_brewster    - Brewster angle in radians (scalar)
%   distancia_brewster - Angular distance to Brewster angle in degrees [1 x N]
%   transmitancia_potencia - Power transmittance [1 x N] (for radar equation)
%   reflectancia_potencia  - Power reflectance [1 x N] (for radar equation)
%
% DESCRIPTION:
%   This function calculates the transmission and reflection coefficients for 
%   TM (Transverse Magnetic, p-polarization) electromagnetic waves using 
%   the Fresnel equations at a planar interface between two lossless dielectric
%   media.
%
%   Ecuaciones de Fresnel para polarización TM (campos eléctricos):
%   r_TM = (n2*cos(θ1) - n1*cos(θ2)) / (n2*cos(θ1) + n1*cos(θ2))
%   t_TM = (2*n1*cos(θ1)) / (n2*cos(θ1) + n1*cos(θ2))
%
%   Transmitancia y Reflectancia de POTENCIA (para ecuaciones de radar):
%   T = (n2*cos(θ2))/(n1*cos(θ1)) * |t_TM|²
%   R = |r_TM|²
%   donde T + R = 1 (conservación de energía)
%
%   Ángulo de Brewster: θ_B = arctan(n2/n1)
%   En este ángulo, r_TM = 0 y t_TM es máximo.
%
%   donde θ1 es el ángulo de incidencia y θ2 es el ángulo de transmisión
%   calculado usando la ley de Snell: n1*sin(θ1) = n2*sin(θ2)

% Número de ángulos de entrada
N = length(angulo_incidencia);

% Calcular el ángulo de Brewster
% Para polarización TM: θ_B = arctan(n2/n1)
angulo_brewster = atan(n2/n1);

% Inicializar salidas
coef_Transm = zeros(1, N);
coef_Reflex = zeros(1, N);
distancia_brewster = zeros(1, N);
transmitancia_potencia = zeros(1, N);
reflectancia_potencia = zeros(1, N);

for i = 1:N
    theta1 = angulo_incidencia(i);
    
    % Calcular distancia al ángulo de Brewster
    distancia_brewster(i) = abs(theta1 - angulo_brewster) * 180/pi; % en grados
    
    % Verificar ángulo crítico para evitar reflexión total interna
    sin_theta_critico = n2 / n1;
    
    if sin(theta1) >= sin_theta_critico && n1 > n2
        % Reflexión total interna - no hay transmisión
        coef_Transm(i) = 0;
        coef_Reflex(i) = 1; % Reflexión total
        transmitancia_potencia(i) = 0;
        reflectancia_potencia(i) = 1;
        warning('Reflexión total interna detectada en el ángulo %.2f°', theta1*180/pi);
    else
        % Calcular ángulo de transmisión usando la ley de Snell
        sin_theta2 = (n1 / n2) * sin(theta1);
        
        % Verificar que sin_theta2 esté en el rango válido [-1, 1]
        if abs(sin_theta2) > 1
            coef_Transm(i) = 0;
            coef_Reflex(i) = 1;
            transmitancia_potencia(i) = 0;
            reflectancia_potencia(i) = 1;
            warning('Ángulo de transmisión inválido para θ1 = %.2f°', theta1*180/pi);
        else
            theta2 = asin(sin_theta2);
            cos_theta1 = cos(theta1);
            cos_theta2 = cos(theta2);
            
            % Denominador común en las ecuaciones de Fresnel (forma CORRECTA)
            denominador = n2 * cos_theta1 + n1 * cos_theta2;
            
            if abs(denominador) < eps
                coef_Transm(i) = 0;
                coef_Reflex(i) = 1;
                transmitancia_potencia(i) = 0;
                reflectancia_potencia(i) = 1;
                warning('Denominador muy pequeño en las ecuaciones de Fresnel');
            else
                % Coeficiente de reflexión TM usando ecuaciones de Fresnel (CORRECTA)
                % r_TM = (n2*cos(θ1) - n1*cos(θ2)) / (n2*cos(θ1) + n1*cos(θ2))
                numerador_r = n2 * cos_theta1 - n1 * cos_theta2;
                coef_Reflex(i) = numerador_r / denominador;
                
                % Coeficiente de transmisión TM usando ecuaciones de Fresnel (CORRECTA)
                % t_TM = (2*n1*cos(θ1)) / (n2*cos(θ1) + n1*cos(θ2))
                numerador_t = 2 * n1 * cos_theta1;
                coef_Transm(i) = numerador_t / denominador;
                
                % CÁLCULO DE TRANSMITANCIA Y REFLECTANCIA DE POTENCIA
                % Transmitancia de potencia: T = (n2*cos(θ2))/(n1*cos(θ1)) * |t_TM|²
                transmitancia_potencia(i) = (n2 * cos_theta2) / (n1 * cos_theta1) * abs(coef_Transm(i))^2;
                
                % Reflectancia de potencia: R = |r_TM|²
                reflectancia_potencia(i) = abs(coef_Reflex(i))^2;
                
                % Verificar conservación de energía: T + R ≈ 1
                conservacion_energia = transmitancia_potencia(i) + reflectancia_potencia(i);
                if abs(conservacion_energia - 1) > 0.01  % Tolerancia 1%
                    warning('Posible violación de conservación de energía: T + R = %.4f (ángulo %.1f°)', ...
                        conservacion_energia, theta1*180/pi);
                end
            end
        end
    end
    
end

% Verificar que los coeficientes estén en un rango físicamente realista
coef_Transm(coef_Transm < 0) = 0;
coef_Transm(coef_Transm > 2) = 2; % Límite teórico máximo para coef. de transmisión

% Los coeficientes de reflexión deben estar entre -1 y 1
coef_Reflex(coef_Reflex < -1) = -1;
coef_Reflex(coef_Reflex > 1) = 1;

% Verificar que las transmitancias/reflectancias de potencia estén en rango físico
transmitancia_potencia(transmitancia_potencia < 0) = 0;
transmitancia_potencia(transmitancia_potencia > 1) = 1; % Debe ser ≤ 1 por conservación de energía
reflectancia_potencia(reflectancia_potencia < 0) = 0;
reflectancia_potencia(reflectancia_potencia > 1) = 1;

% Resumen de resultados
% fprintf('\n=== ANÁLISIS DEL ÁNGULO DE BREWSTER ===\n');
% fprintf('Ángulo de Brewster: %.2f° (%.4f rad)\n', angulo_brewster*180/pi, angulo_brewster);
% fprintf('Rango de ángulos analizados: %.2f° a %.2f°\n', ...
%         min(angulo_incidencia)*180/pi, max(angulo_incidencia)*180/pi);

% Encontrar el ángulo más cercano al de Brewster
[min_dist, idx_min] = min(distancia_brewster);
% fprintf('Ángulo más cercano a Brewster: %.2f° (distancia: %.3f°)\n', ...
%         angulo_incidencia(idx_min)*180/pi, min_dist);
% fprintf('En este ángulo: r_TM = %.4f, t_TM = %.4f\n', ...
%         coef_Reflex(idx_min), coef_Transm(idx_min));

end