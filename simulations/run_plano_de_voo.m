clc;
clear;
close all;
%
addpath(genpath('gs'))
addpath(genpath('proc'))
addpath(genpath('tools'))
addpath(genpath('flightpath'))

% Parametros ambiente
systemJSON = json2struct(strcat('parametros',filesep,'system_espiral_plano_voo.json'));
strSystem  = systemJSON.system; clear systemJSON;

% Parametros del transmisor
radarJSON  = json2struct(strcat('parametros',filesep,'radarTx_espiral_plano_voo.json'));
strRadarTx = radarJSON.radar; clear radarJSON;

targetJSON = json2struct(strcat('parametros',filesep,'target_espiral_plano_voo.json'));
strTarget  = targetJSON.target; clear targetJSON;

% Transmissor
[PxT_, PyT_, PzT_, t] = funcao_espiral(strRadarTx.NumVoltasEsp, strRadarTx.RaioMenorEsp, strRadarTx.RaioMaiorEsp,...
                                    strRadarTx.AltMaiorEsp,strRadarTx.AltMenorEsp, strRadarTx.Vt, strRadarTx.PRF, strRadarTx.NorthOffset);

                                
% Decimar posiciones del transmisor para acelerar optimización
decimationFactor = 200;
PxT = PxT_(1:decimationFactor:end);
PyT = PyT_(1:decimationFactor:end);
PzT = PzT_(1:decimationFactor:end);
t_  = t(1:decimationFactor:end);
                                    
Tx_pos = [PxT; PyT; PzT];

gridTarget = strTarget.grid;
tg         = createGridTarget(gridTarget.xSize, gridTarget.ySize, gridTarget.zMin, gridTarget.zMax, gridTarget.nx, gridTarget.ny, gridTarget.nz).';

% Parámetros de radar
Pt     = strRadarTx.PotenciaTx;                            % Potencia transmitida [W]
lambda = physconst('lightspeed')/strRadarTx.FreqPortadora; % Longitud de onda [m] (400 MHz)
sigma  = strTarget.rcs;                                    % RCS del target [m²]

% Antena
radPattern = createRadiationPattern(strRadarTx.AperturaElev, strRadarTx.AperturaAzimut);

%% Inicializar resultados
Potencia_resultados = struct();
Potencia_resultados.Tx_positions = Tx_pos;
Potencia_resultados.targets = tg;
Potencia_resultados.parametros.n1 = 1; % Siempre aire
Potencia_resultados.parametros.n2 = strSystem.IndiceRefracaoSolo;
Potencia_resultados.parametros.Pt = Pt;
Potencia_resultados.parametros.lambda = lambda;
Potencia_resultados.parametros.sigma = sigma;
%% Optimización del receptor bistático para CADA posición del transmisor

% Ángulo de Brewster para orientar la búsqueda
angulo_brewster = atan(strSystem.IndiceRefracaoSolo/1);
fprintf('\nÁngulo de Brewster: %.1f°\n', angulo_brewster*180/pi);

fprintf('Optimizando posición del receptor para cada posición. Número de posiciones Tx: %d\n', length(PxT));
fprintf('Número de targets en el volumen: %d\n', size(tg,2));

% Configuración del optimizador
lb = [-strSystem.searchRadius, -strSystem.searchRadius]; % Límites inferiores
ub = [strSystem.searchRadius, strSystem.searchRadius];   % Límites superiores

% Estrategia de puntos iniciales inteligente
% 1. Punto cerca del transmisor
% 2. Punto en dirección opuesta a los targets (principio de bistático)
% 3. Punto cerca del ángulo de Brewster teórico
% 4. Puntos aleatorios para explorar

n_starts = 5;  % Total de puntos iniciales por posición del transmisor

% Opciones de optimización
options = optimoptions('fmincon', 'Display', 'none', ...
    'MaxIterations', 200, 'FunctionTolerance', 1e-10, ...
    'Algorithm', 'interior-point');

% Arrays para almacenar resultados
Rx_opt_all    = zeros(3, length(PxT)); % Posiciones óptimas del receptor
max_power_all = zeros(1, length(PxT)); % Potencias máximas

% Optimización para cada posición del transmisor
for tx_idx = 1:length(PxT)
    
    fprintf('\nOptimizando para posición Tx %d/%d: [%.1f, %.1f, %.1f]...\n', tx_idx, length(PxT), PxT(tx_idx), PyT(tx_idx), PzT(tx_idx));

    % Posición actual del transmisor
    Tx_pos_current = [PxT(tx_idx); PyT(tx_idx); PzT(tx_idx)];
    Rx_z           = PzT(tx_idx); % Misma altitud que esta posición del transmisor

    % 1) Calcular ganancias del transmisor para todos los targets
    [az_Tx_ref_current, el_Tx_ref_current, Gt_current] = calculateTxGainsForTargets(Tx_pos_current, tg, radPattern);
    
    % Estrategia de puntos iniciales inteligente
    results  = cell(n_starts, 1);
    fvals    = zeros(n_starts, 1);
    
    % 2) Preparar geometría para elegir posiciones candidatas de Rx
    search_geo = buildRxSearchGeometry(Tx_pos_current, Rx_z, tg, angulo_brewster);
   
    for i = 1:n_starts
        % 3) Elegir posición candidata de inicio
        x0 = chooseCandidateRxPos(i, search_geo.target_center, search_geo.targets_distance, ...
            search_geo.opposite_dir, search_geo.perp_dir, search_geo.distancia_brewster, search_geo.azimut_opuesto);
        
        % Asegurar que esté dentro de los límites
        x0 = max(lb, min(ub, x0));
          
        % 4) Optimizar desde la posición candidata
        [x_opt, fval] = fmincon(@(xy) objective_function(xy, tg, Tx_pos_current, Rx_z, radPattern, ...
            az_Tx_ref_current, el_Tx_ref_current, Gt_current, Pt, sigma, lambda, 1, strSystem.IndiceRefracaoSolo), ...
            x0, [], [], [], [], lb, ub, [], options);
        
        results{i} = x_opt;
        fvals(i) = fval;
           
    end
    
    % Seleccionar mejor resultado para esta posición Tx
    [~, best_idx] = min(fvals);
    Rx_opt_all(:, tx_idx) = [results{best_idx}(1); results{best_idx}(2); Rx_z];
    max_power_all(tx_idx) = -fvals(best_idx);

end

fprintf('\n=== RESULTADOS COMPLETOS DE OPTIMIZACIÓN ===\n');
for tx_idx = 1:length(PxT)
    fprintf('Tx %d: [%.1f,%.1f,%.1f] → Rx óptimo: [%.2f,%.2f,%.1f] - Potencia: %.1f dBm\n', ...
        tx_idx, PxT(tx_idx), PyT(tx_idx), PzT(tx_idx), ...
        Rx_opt_all(:, tx_idx), 10*log10(max_power_all(tx_idx)*1000));
end

% Agregar resultados de optimización a la estructura de potencia
Potencia_resultados.optimizacion.Rx_optimo = Rx_opt_all;
Potencia_resultados.optimizacion.potencia_maxima = max_power_all;
Potencia_resultados.optimizacion.potencia_maxima_dBm = 10*log10(max_power_all*1000);
Potencia_resultados.fecha_calculo = datestr(now);

%% Guardar resultados en archivos
fprintf('\n=== GUARDANDO RESULTADOS ===\n');

outputDir = fullfile('io', 'plan_vuelo');

% Crear carpeta de resultados si no existe
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Nombre del archivo con timestamp
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename_base = fullfile(outputDir, sprintf('potencia_calculo_%s', timestamp));

% Guardar estructura completa en .mat
save([filename_base '.mat'], 'Potencia_resultados');
fprintf('Datos completos guardados en: %s.mat\n', filename_base);

% Guardar resultados de optimización en CSV separado
fid_opt = fopen([filename_base '_optimizacion.csv'], 'w');
fprintf(fid_opt, 'Tx_ID,Tx_X,Tx_Y,Tx_Z,Rx_opt_X,Rx_opt_Y,Rx_opt_Z,Potencia_maxima_W,Potencia_maxima_dBm\n');
for tx_idx = 1:length(PxT)
    fprintf(fid_opt, '%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.6e,%.3f\n', ...
        tx_idx, ...
        Tx_pos(1,tx_idx), Tx_pos(2,tx_idx), Tx_pos(3,tx_idx), ...
        Rx_opt_all(1,tx_idx), Rx_opt_all(2,tx_idx), Rx_opt_all(3,tx_idx), ...
        max_power_all(tx_idx), ...
        10*log10(max_power_all(tx_idx)*1000));
end
fclose(fid_opt);
fprintf('Optimización CSV guardado en: %s_optimizacion.csv\n', filename_base);

fprintf('\n=== ARCHIVOS GENERADOS ===\n');
fprintf('1. %s.mat - Estructura completa con todos los datos\n', filename_base);
fprintf('2. %s_optimizacion.csv - Posiciones óptimas de receptor y potencias máximas\n', filename_base);


%% Visualización final de la configuración bistática optimizada
plot_bistatic_configuration(PxT, PyT, PzT, tg, Rx_opt_all, 1, strSystem.IndiceRefracaoSolo, angulo_brewster);

%% Grabar trayectoria completa (extrapolado), extrapolar R_opt_all, tiempos iniciales t_, tiempos a los cuales interpolar: t
R_opt_x = interp1(t_, Rx_opt_all(1,:), t, 'pchip', 'extrap');
R_opt_y = interp1(t_, Rx_opt_all(2,:), t, 'pchip', 'extrap');
R_opt_z = interp1(t_, Rx_opt_all(3,:), t, 'pchip', 'extrap');

R_opt_all_extrap = [R_opt_x; R_opt_y; R_opt_z];
Tx_all           = [PxT_; PyT_; PzT_];

%% Salvar en -mat (/io/plan_vuelo)
save(fullfile(outputDir, 'R_opt_all_extrap.mat'),'R_opt_all_extrap');
save(fullfile(outputDir, 'Tx_all.mat'),'Tx_all');
