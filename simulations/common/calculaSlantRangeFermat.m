function [strReflexao, strRefraccoes] = calculaSlantRangeFermat(strDEM, Tx ,Rx, P, n1, n2, bPlotVerbose)

c = physconst('lightspeed');

v1 = c/n1;
v2 = c/n2;

X_vec = strDEM.X_vec;
Y_vec = strDEM.Y_vec;
Z_DEM = strDEM.Z_DEM;

num_pares = size(Tx, 1);

P_REF = [];
Q_REF = [];
R_PRX = [];

% --- Pre-asignación de los arrays de estructuras de salida ---
strReflexao   = struct('P_reflex', [], 'Ang_incidencia_deg', [], 'Ang_reflexion_deg', []);
strRefraccoes = struct('P_refrac_ida', [], 'P_refrac_volta', [], ...
    'Ang_incidencia_ida_deg', [], 'Ang_refrac_ida_deg', [], ...
    'Ang_incidencia_volta_deg', [], 'Ang_refrac_volta_deg', []);

% --- CREACIÓN DEL OBJETO DE INTERPOLACIÓN ---
interfaz_interpolada = griddedInterpolant({Y_vec, X_vec}, double(Z_DEM), 'linear', 'none');


if bPlotVerbose
    figure;
    hold on; grid on; axis equal; view(45, 30);
    
    % Dibujar la superficie una sola vez
    surf(X_vec, Y_vec, Z_DEM, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    colormap('winter'); colorbar;
    title('Principio de Fermat en un DEM (Múltiples Trayectorias)');
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Elevación Z (m)');
    
    % Dibujar el punto P (es común para todos)
    plot3(P(1), P(2), P(3), 'p', 'MarkerSize', 15, 'MarkerFaceColor', 'k', 'DisplayName', 'Punto C (P)');
end

% BUCLE PRINCIPAL: Iterar sobre cada par Tx/Rx
for i = 1:num_pares

%     % Porcentaje de progreso cada 10%
%     if mod(i, round(num_pares/10)) == 0
%         fprintf('--- PROCESANDO PAR %d de %d ---\n', i, num_pares);
%     end

    Tx_actual = Tx(i, :);
    Rx_actual = Rx(i, :);
        
    % Punto inicial para la búsqueda (centro del DEM)
    punto_inicial = [mean(X_vec), mean(Y_vec)];
    
    % --- REFLEXIÓN (A -> S -> B) ---
    fun_reflex = @(Pxy) tiempo_total_reflex(Tx_actual, Rx_actual, Pxy, interfaz_interpolada, v1);
    Pxy_reflex = fminsearch(fun_reflex, punto_inicial);
    z_reflex = interfaz_interpolada(Pxy_reflex(2), Pxy_reflex(1));
    Preflex = [Pxy_reflex, z_reflex];
    
    % --- REFRACCIÓN (A -> Q -> P) ---
    fun_refract = @(Qxy) tiempo_total_refract(Tx_actual, P, Qxy, interfaz_interpolada, v1, v2);
    Qxy_refract = fminsearch(fun_refract, punto_inicial);
    z_refract = interfaz_interpolada(Qxy_refract(2), Qxy_refract(1));
    Qrefract = [Qxy_refract, z_refract];
    
    % --- REFRACCIÓN INVERSA (P -> R -> B) ---
    fun_refract_PRx = @(Rxy) tiempo_total_refract_inversa(P, Rx_actual, Rxy, interfaz_interpolada, v2, v1);
    Rxy_PRx = fminsearch(fun_refract_PRx, punto_inicial);
    z_PRx = interfaz_interpolada(Rxy_PRx(2), Rxy_PRx(1));
    R_PRx = [Rxy_PRx, z_PRx];
    
    %% Cálculo de Ángulos (en grados)
    % Reflexión
    N_reflex = calcular_normal(Pxy_reflex, interfaz_interpolada);
    rayo_incidente_refl = (Tx_actual - Preflex) / norm(Tx_actual - Preflex);
    rayo_reflejado = (Rx_actual - Preflex) / norm(Rx_actual - Preflex);
    ang_inc_refl_deg = rad2deg(acos(abs(dot(rayo_incidente_refl, N_reflex))));
    ang_refl_deg = rad2deg(acos(abs(dot(rayo_reflejado, N_reflex))));
    
    % Refracción Ida
    N_refract = calcular_normal(Qxy_refract, interfaz_interpolada);
    rayo_incidente_refr = (Tx_actual - Qrefract) / norm(Tx_actual - Qrefract);
    rayo_refractado_QP = (P - Qrefract) / norm(P - Qrefract);
    ang_inc_refr_deg = rad2deg(acos(abs(dot(rayo_incidente_refr, N_refract))));
    ang_refr_QP_deg = rad2deg(acos(abs(dot(rayo_refractado_QP, -N_refract))));
    
    % Refracción Vuelta
    N_refract_inv = calcular_normal(Rxy_PRx, interfaz_interpolada);
    rayo_incidente_PR = (P - R_PRx) / norm(P - R_PRx);
    rayo_refractado_RRx = (Rx_actual - R_PRx) / norm(Rx_actual - R_PRx);
    ang_inc_PR_deg = rad2deg(acos(abs(dot(rayo_incidente_PR, -N_refract_inv))));
    ang_refr_RRx_deg = rad2deg(acos(abs(dot(rayo_refractado_RRx, N_refract_inv))));
    
    % Append posicciones 
    P_REF = [P_REF, Preflex.'];
    Q_REF = [Q_REF, Qrefract.'];
    R_PRX = [R_PRX, R_PRx.'];
    
    if bPlotVerbose
        % Puntos Tx y Rx del par actual
        plot3(Tx_actual(1), Tx_actual(2), Tx_actual(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        plot3(Rx_actual(1), Rx_actual(2), Rx_actual(3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
        
        % Rayos
        plot3([Tx_actual(1), Preflex(1)], [Tx_actual(2), Preflex(2)], [Tx_actual(3), Preflex(3)], 'r-', 'LineWidth', 2);
        plot3([Preflex(1), Rx_actual(1)], [Preflex(2), Rx_actual(2)], [Preflex(3), Rx_actual(3)], 'g-', 'LineWidth', 2);
        plot3([Tx_actual(1), Qrefract(1)], [Tx_actual(2), Qrefract(2)], [Tx_actual(3), Qrefract(3)], 'm--', 'LineWidth', 2);
        plot3([Qrefract(1), P(1)], [Qrefract(2), P(2)], [Qrefract(3), P(3)], 'c--', 'LineWidth', 2);
        plot3([P(1), R_PRx(1)], [P(2), R_PRx(2)], [P(3), R_PRx(3)], '-.', 'Color', [0 0.4 0.8], 'LineWidth', 2);
        plot3([R_PRx(1), Rx_actual(1)], [R_PRx(2), Rx_actual(2)], [R_PRx(3), Rx_actual(3)], '-.', 'Color', [0.1 0.8 0.1], 'LineWidth', 2);
        
        % Opcional: Vectores Normales (puede saturar el gráfico)
        escala_normal = 40;
        quiver3(Preflex(1), Preflex(2), Preflex(3), N_reflex(1), N_reflex(2), N_reflex(3), escala_normal, 'k', 'LineWidth', 1);
        quiver3(Qrefract(1), Qrefract(2), Qrefract(3), N_refract(1), N_refract(2), N_refract(3), escala_normal, 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
        quiver3(R_PRx(1), R_PRx(2), R_PRx(3), -N_refract_inv(1), -N_refract_inv(2), -N_refract_inv(3), escala_normal, 'k', 'LineWidth', 2, 'HandleVisibility', 'off');

        
        % Imprimir resultados en consola
        fprintf('Punto óptimo (S): (%.2f, %.2f, %.2f)\n', Preflex);
        fprintf('Ángulo incidencia: %.2f° | Ángulo reflexión: %.2f°\n', ang_inc_refl_deg, ang_refl_deg);
        fprintf('Verificación Snell (Ida): n1*sin(θ1)=%.3f | n2*sin(θ2)=%.3f\n', n1*sind(ang_inc_refr_deg), n2*sind(ang_refr_QP_deg));
        fprintf('Verificación Snell (Vuelta): n2*sin(θ2)=%.3f | n1*sin(θ1)=%.3f\n\n', n2*sind(ang_inc_PR_deg), n1*sind(ang_refr_RRx_deg));
    end
    
%     disp(strcat("[",num2str(Tx_actual),"] -> Refract: ","[",num2str(Qrefract),"] Target: ","[",num2str(P),"]"))
%     disp(strcat("[",num2str(Rx_actual),"] -> Refract: ","[",num2str(R_PRx),"] Target: ","[",num2str(P),"]"))
%     
%     r1Tx = sqrt(sum((Tx_actual - Qrefract).^2));
%     r2Tx = sqrt(sum((P - Qrefract).^2));
%     r1Rx = sqrt(sum((Rx_actual - Qrefract).^2));
%     r2Rx = sqrt(sum((P - R_PRx).^2));
%     
%     disp(strcat("Ranges: r1Tx=",num2str(r1Tx),", r2Tx=",num2str(r2Tx),", r1Rx=",num2str(r1Rx),", r2Rx=",num2str(r2Rx)))
%     disp("----------------------------------------------------------------------------------")


end

disp("----------------------------------------------------------------------------------")

%% Almacenamiento de Resultados en las Estructuras
strReflexao.P_reflex = P_REF;
strReflexao.Ang_incidencia_deg = 0; % ang_inc_refl_deg;
strReflexao.Ang_reflexion_deg = 0;  % ang_refl_deg;

strRefraccoes.P_refrac_ida = Q_REF;
strRefraccoes.P_refrac_volta = R_PRX;
strRefraccoes.Ang_incidencia_ida_deg = 0;   % ang_inc_refr_deg;
strRefraccoes.Ang_refrac_ida_deg = 0;       % ang_refr_QP_deg;
strRefraccoes.Ang_incidencia_volta_deg = 0; % ang_inc_PR_deg;
strRefraccoes.Ang_refrac_volta_deg = 0;     % ang_refr_RRx_deg;

end