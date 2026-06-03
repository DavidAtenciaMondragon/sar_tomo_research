function [strReflexao, strRefraccoes] = calculaSlantRangeFermat_optimized(strDEM, Tx, Rx, P, n1, n2, bPlotVerbose, factor_interp)
% CALCULASLANTRANGEFERMAT_OPTIMIZED - Versión optimizada del cálculo de rango oblicuo usando el principio de Fermat
%
% Esta versión combina la precisión del método original con la eficiencia de la aproximación vectorizada
%
% Parámetros:
%   strDEM: Estructura con el modelo digital de elevación
%   Tx, Rx: Posiciones del transmisor y receptor
%   P: Punto de destino
%   n1, n2: Índices de refracción
%   bPlotVerbose: Flag para visualización
%   factor_interp: Factor de interpolación (opcional, default=2)

if nargin < 8
    factor_interp = 2; % Factor de interpolación por defecto
end

% Constantes físicas
c = physconst('lightspeed');
v1 = c/n1;
v2 = c/n2;

% Extraer datos del DEM
X_vec = strDEM.X_vec;
Y_vec = strDEM.Y_vec;
Z_DEM = strDEM.Z_DEM;
num_pares = size(Tx, 1);

% OPTIMIZACIÓN 1: Crear grilla interpolada más densa para mejor precisión
sizeGrid = length(X_vec);
x_interp = linspace(min(X_vec), max(X_vec), sizeGrid * factor_interp);
y_interp = linspace(min(Y_vec), max(Y_vec), sizeGrid * factor_interp);
[x_grid_interp, y_grid_interp] = meshgrid(x_interp, y_interp);

% Interpolación vectorizada del terreno
z_interp = interp2(X_vec, Y_vec, Z_DEM, x_grid_interp, y_grid_interp, 'linear');

% OPTIMIZACIÓN 2: Pre-calcular todas las posiciones de interfaz
pos_interface = [x_grid_interp(:), y_grid_interp(:), z_interp(:)];
valid_points = ~isnan(pos_interface(:,3)); % Filtrar puntos válidos
pos_interface = pos_interface(valid_points, :);

% OPTIMIZACIÓN 3: Pre-asignación de estructuras de salida
strReflexao = struct('P_reflex', zeros(3, num_pares), ...
                     'Ang_incidencia_deg', zeros(1, num_pares), ...
                     'Ang_reflexion_deg', zeros(1, num_pares));
strRefraccoes = struct('P_refrac_ida', zeros(3, num_pares), ...
                       'P_refrac_volta', zeros(3, num_pares), ...
                       'Ang_incidencia_ida_deg', zeros(1, num_pares), ...
                       'Ang_refrac_ida_deg', zeros(1, num_pares), ...
                       'Ang_incidencia_volta_deg', zeros(1, num_pares), ...
                       'Ang_refrac_volta_deg', zeros(1, num_pares));

% Preparar visualización si es necesario
if bPlotVerbose
    figure;
    hold on; grid on; axis equal; view(45, 30);
    surf(X_vec, Y_vec, Z_DEM, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    colormap('winter'); colorbar;
    title('Principio de Fermat Optimizado - Múltiples Trayectorias');
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Elevación Z (m)');
    plot3(P(1), P(2), P(3), 'p', 'MarkerSize', 15, 'MarkerFaceColor', 'k', 'DisplayName', 'Punto C (P)');
end

% OPTIMIZACIÓN 4: Procesamiento vectorizado por lotes
batch_size = min(num_pares, 100); % Procesar en lotes para control de memoria
for batch_start = 1:batch_size:num_pares
    batch_end = min(batch_start + batch_size - 1, num_pares);
    batch_indices = batch_start:batch_end;
    batch_size_actual = length(batch_indices);
    
    % Extraer lote actual
    Tx_batch = Tx(batch_indices, :);
    Rx_batch = Rx(batch_indices, :);
    
    % OPTIMIZACIÓN 5: Cálculo vectorizado de distancias
    % Para reflexión: Tx -> interfaz -> Rx
    dist_Tx_interface = pdist2(pos_interface, Tx_batch);  % [N_interface x batch_size]
    dist_interface_Rx = pdist2(pos_interface, Rx_batch);  % [N_interface x batch_size]
    tiempo_reflex = (dist_Tx_interface + dist_interface_Rx) / v1;
    
    % Para refracción Tx->P: Tx -> interfaz -> P  
    dist_interface_P = pdist2(pos_interface, P);  % [N_interface x 1] porque P es un solo punto
    % Expandir dist_interface_P para coincidir con batch_size
    dist_interface_P_expanded = repmat(dist_interface_P, 1, batch_size_actual);  % [N_interface x batch_size]
    tiempo_Tx_P = dist_Tx_interface / v1 + dist_interface_P_expanded / v2;
    
    % Para refracción P->Rx: P -> interfaz -> Rx
    % Corregir: P es una fila, pos_interface son filas también
    dist_P_interface = pdist2(P, pos_interface);  % [1 x N_interface] porque P es un solo punto
    % Transponer y expandir para coincidir con las otras matrices
    dist_P_interface_col = dist_P_interface';  % [N_interface x 1]
    dist_P_interface_expanded = repmat(dist_P_interface_col, 1, batch_size_actual);  % [N_interface x batch_size]
    tiempo_P_Rx = dist_P_interface_expanded / v2 + dist_interface_Rx / v1;
    
    % OPTIMIZACIÓN 6: Búsqueda vectorizada del mínimo
    [~, idx_reflex] = min(tiempo_reflex, [], 1);
    [~, idx_Tx_P] = min(tiempo_Tx_P, [], 1);
    [~, idx_P_Rx] = min(tiempo_P_Rx, [], 1);
    
    % Extraer puntos óptimos
    for i = 1:batch_size_actual
        batch_idx = batch_indices(i);
        
        % Puntos óptimos de intersección
        P_reflex = pos_interface(idx_reflex(i), :);
        Q_refract = pos_interface(idx_Tx_P(i), :);
        R_PRx = pos_interface(idx_P_Rx(i), :);
        
        % Almacenar resultados
        strReflexao.P_reflex(:, batch_idx) = P_reflex';
        strRefraccoes.P_refrac_ida(:, batch_idx) = Q_refract';
        strRefraccoes.P_refrac_volta(:, batch_idx) = R_PRx';
        
        % OPTIMIZACIÓN 7: Cálculo opcional de ángulos (solo si se necesitan)
        if nargout > 1 || bPlotVerbose
            % Calcular ángulos solo cuando sea necesario
            [ang_inc_refl, ang_refl] = calcular_angulos_reflexion(Tx_batch(i,:), Rx_batch(i,:), P_reflex, x_interp, y_interp, z_interp);
            [ang_inc_refr, ang_refr_QP, ang_inc_PR, ang_refr_RRx] = calcular_angulos_refraccion(Tx_batch(i,:), P, Rx_batch(i,:), Q_refract, R_PRx, x_interp, y_interp, z_interp);
            
            strReflexao.Ang_incidencia_deg(batch_idx) = ang_inc_refl;
            strReflexao.Ang_reflexion_deg(batch_idx) = ang_refl;
            strRefraccoes.Ang_incidencia_ida_deg(batch_idx) = ang_inc_refr;
            strRefraccoes.Ang_refrac_ida_deg(batch_idx) = ang_refr_QP;
            strRefraccoes.Ang_incidencia_volta_deg(batch_idx) = ang_inc_PR;
            strRefraccoes.Ang_refrac_volta_deg(batch_idx) = ang_refr_RRx;
        end
        
        % Visualización para el lote actual
        if bPlotVerbose
            plot3(Tx_batch(i,1), Tx_batch(i,2), Tx_batch(i,3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            plot3(Rx_batch(i,1), Rx_batch(i,2), Rx_batch(i,3), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
            
            % Rayos de reflexión
            plot3([Tx_batch(i,1), P_reflex(1)], [Tx_batch(i,2), P_reflex(2)], [Tx_batch(i,3), P_reflex(3)], 'r-', 'LineWidth', 1.5);
            plot3([P_reflex(1), Rx_batch(i,1)], [P_reflex(2), Rx_batch(i,2)], [P_reflex(3), Rx_batch(i,3)], 'g-', 'LineWidth', 1.5);
            
            % Rayos de refracción
            plot3([Tx_batch(i,1), Q_refract(1)], [Tx_batch(i,2), Q_refract(2)], [Tx_batch(i,3), Q_refract(3)], 'm--', 'LineWidth', 1.5);
            plot3([Q_refract(1), P(1)], [Q_refract(2), P(2)], [Q_refract(3), P(3)], 'c--', 'LineWidth', 1.5);
            plot3([P(1), R_PRx(1)], [P(2), R_PRx(2)], [P(3), R_PRx(3)], '-.', 'Color', [0 0.4 0.8], 'LineWidth', 1.5);
            plot3([R_PRx(1), Rx_batch(i,1)], [R_PRx(2), Rx_batch(i,2)], [R_PRx(3), Rx_batch(i,3)], '-.', 'Color', [0.1 0.8 0.1], 'LineWidth', 1.5);
        end
    end
end

end

% FUNCIÓN AUXILIAR: Cálculo eficiente de ángulos de reflexión
function [ang_inc_deg, ang_refl_deg] = calcular_angulos_reflexion(Tx, Rx, P_reflex, x_interp, y_interp, z_interp)
    % Calcular normal usando diferencias finitas optimizadas
    h = abs(x_interp(2) - x_interp(1)); % Usar resolución de la grilla
    
    % Interpolación local para gradiente
    F_interp = griddedInterpolant({y_interp, x_interp}, z_interp, 'linear', 'none');
    
    dzdx = (F_interp(P_reflex(2), P_reflex(1) + h) - F_interp(P_reflex(2), P_reflex(1) - h)) / (2*h);
    dzdy = (F_interp(P_reflex(2) + h, P_reflex(1)) - F_interp(P_reflex(2) - h, P_reflex(1))) / (2*h);
    
    N = [-dzdx, -dzdy, 1];
    N_unit = N / norm(N);
    
    % Vectores de rayo
    rayo_incidente = (Tx - P_reflex) / norm(Tx - P_reflex);
    rayo_reflejado = (Rx - P_reflex) / norm(Rx - P_reflex);
    
    % Ángulos
    ang_inc_deg = rad2deg(acos(abs(dot(rayo_incidente, N_unit))));
    ang_refl_deg = rad2deg(acos(abs(dot(rayo_reflejado, N_unit))));
end

% FUNCIÓN AUXILIAR: Cálculo eficiente de ángulos de refracción  
function [ang_inc_ida_deg, ang_refr_ida_deg, ang_inc_volta_deg, ang_refr_volta_deg] = calcular_angulos_refraccion(Tx, P, Rx, Q_refract, R_PRx, x_interp, y_interp, z_interp)
    h = abs(x_interp(2) - x_interp(1));
    F_interp = griddedInterpolant({y_interp, x_interp}, z_interp, 'linear', 'none');
    
    % Normal en Q (refracción ida)
    dzdx_Q = (F_interp(Q_refract(2), Q_refract(1) + h) - F_interp(Q_refract(2), Q_refract(1) - h)) / (2*h);
    dzdy_Q = (F_interp(Q_refract(2) + h, Q_refract(1)) - F_interp(Q_refract(2) - h, Q_refract(1))) / (2*h);
    N_Q = [-dzdx_Q, -dzdy_Q, 1] / norm([-dzdx_Q, -dzdy_Q, 1]);
    
    % Ángulos en Q
    rayo_inc_Q = (Tx - Q_refract) / norm(Tx - Q_refract);
    rayo_refr_Q = (P - Q_refract) / norm(P - Q_refract);
    ang_inc_ida_deg = rad2deg(acos(abs(dot(rayo_inc_Q, N_Q))));
    ang_refr_ida_deg = rad2deg(acos(abs(dot(rayo_refr_Q, -N_Q))));
    
    % Normal en R (refracción vuelta)
    dzdx_R = (F_interp(R_PRx(2), R_PRx(1) + h) - F_interp(R_PRx(2), R_PRx(1) - h)) / (2*h);
    dzdy_R = (F_interp(R_PRx(2) + h, R_PRx(1)) - F_interp(R_PRx(2) - h, R_PRx(1))) / (2*h);
    N_R = [-dzdx_R, -dzdy_R, 1] / norm([-dzdx_R, -dzdy_R, 1]);
    
    % Ángulos en R
    rayo_inc_R = (P - R_PRx) / norm(P - R_PRx);
    rayo_refr_R = (Rx - R_PRx) / norm(Rx - R_PRx);
    ang_inc_volta_deg = rad2deg(acos(abs(dot(rayo_inc_R, -N_R))));
    ang_refr_volta_deg = rad2deg(acos(abs(dot(rayo_refr_R, N_R))));
end