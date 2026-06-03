function [Gt_ajustado, Gr_ajustado] = calcular_ganancias_antena(Tx, Rx, P_targets, abertura_el, abertura_az)
    % CALCULAR_GANANCIAS_ANTENA - Calcula las ganancias de antena considerando patrón de radiación
    %
    % Entradas:
    %   Tx - Posición del transmisor [x, y, z]
    %   Rx - Posición del receptor [x, y, z] 
    %   P_targets - Matriz Nx3 con posiciones de objetivos [x, y, z]
    %   abertura_el - Abertura en elevación (grados) [default: 70°]
    %   abertura_az - Abertura en azimut (grados) [default: 20°]
    %
    % Salidas:
    %   Gt_ajustado - Vector con ganancias del transmisor para cada objetivo
    %   Gr_ajustado - Vector con ganancias del receptor para cada objetivo
    
    if nargin < 4
        abertura_el = 70;
    end
    if nargin < 5
        abertura_az = 20;
    end
    
    % CÁLCULO DE GANANCIAS DEL TRANSMISOR (Gt)
    
    % 1. Calcular ángulos azimut (radar siempre apunta hacia 0,0,0)
    ang_Tx_orig = atan2d(-Tx(2), -Tx(1)); % Ángulo original Tx->(0,0)   
    ang_Tx_all_targets = atan2d(P_targets(:,2)-Tx(2), P_targets(:,1)-Tx(1)); % Ángulo Tx->Pi
    
    % Calcular diferencia de ángulos (considerando ciclo completo 0-360)
    diff_ang_Tx = ang_Tx_all_targets - ang_Tx_orig;
    diff_ang_Tx(diff_ang_Tx > 180)  = diff_ang_Tx(diff_ang_Tx > 180) - 360;
    diff_ang_Tx(diff_ang_Tx < -180) = diff_ang_Tx(diff_ang_Tx < -180) + 360;

    % 2. Calcular ángulos elevación (radar siempre apunta hacia 0,0,0)
    elev_Tx_orig = asind(-Tx(3)/norm(Tx)); % Elevación original Tx->(0,0,0)
    elev_Tx_all_targets = asind((P_targets(:,3)-Tx(3))./vecnorm(P_targets - Tx, 2, 2)); % Elevación Tx->Pi
    diff_elev_Tx = elev_Tx_all_targets - elev_Tx_orig;
    
    % 3. Calcular patrón de radiación para Tx
    padrao = crearPadraoRadiacao(abertura_el, abertura_az);
    
    % Mapear ángulos al patrón de radiación
    [~,idx_az_Tx] = min(abs([-90:1:90].' - diff_ang_Tx.'));
    [~,idx_el_Tx] = min(abs([-90:1:90].' - diff_elev_Tx.'));
    
    idx_Gt = sub2ind(size(padrao), idx_el_Tx, idx_az_Tx);
    Gt_ajustado = padrao(idx_Gt); 
    
    % CÁLCULO DE GANANCIAS DEL RECEPTOR (Gr) - Similar pero para Rx
    
    % 1. Calcular ángulos azimut para receptor (apunta hacia 0,0,0)
    ang_Rx_orig = atan2d(-Rx(2), -Rx(1)); % Ángulo original Rx->(0,0)   
    ang_Rx_all_targets = atan2d(P_targets(:,2)-Rx(2), P_targets(:,1)-Rx(1)); % Ángulo Rx->Pi
    
    % Calcular diferencia de ángulos
    diff_ang_Rx = ang_Rx_all_targets - ang_Rx_orig;
    diff_ang_Rx(diff_ang_Rx > 180)  = diff_ang_Rx(diff_ang_Rx > 180) - 360;
    diff_ang_Rx(diff_ang_Rx < -180) = diff_ang_Rx(diff_ang_Rx < -180) + 360;

    % 2. Calcular ángulos elevación para receptor
    elev_Rx_orig = asind(-Rx(3)/norm(Rx)); % Elevación original Rx->(0,0,0)
    elev_Rx_all_targets = asind((P_targets(:,3)-Rx(3))./vecnorm(P_targets - Rx, 2, 2)); % Elevación Rx->Pi
    diff_elev_Rx = elev_Rx_all_targets - elev_Rx_orig;
    
    % 3. Mapear ángulos al patrón de radiación para Rx
    [~,idx_az_Rx] = min(abs([-90:1:90].' - diff_ang_Rx.'));
    [~,idx_el_Rx] = min(abs([-90:1:90].' - diff_elev_Rx.'));
    
    idx_Gr = sub2ind(size(padrao), idx_el_Rx, idx_az_Rx);
    Gr_ajustado = padrao(idx_Gr); 

end