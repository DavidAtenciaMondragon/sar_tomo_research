clc;
clear;
close all;
%
addpath(genpath('gs'))
addpath(genpath('proc'))
addpath(genpath('tools'))
addpath(genpath('flightpath'))

%% Cargar parámetros
systemJSON = json2struct(strcat('parametros',filesep,'system_espiral_plano_voo.json'));
strSystem  = systemJSON.system; clear systemJSON;

radarJSON  = json2struct(strcat('parametros',filesep,'radarTx_espiral_plano_voo.json'));
strRadarTx = radarJSON.radar; clear radarJSON;

targetJSON = json2struct(strcat('parametros',filesep,'target_espiral_plano_voo.json'));
strTarget  = targetJSON.target; clear targetJSON;

%% Trayectoria del transmisor
[PxT_, PyT_, PzT_, t] = funcao_espiral(strRadarTx.NumVoltasEsp, strRadarTx.RaioMenorEsp, strRadarTx.RaioMaiorEsp,...
                                    strRadarTx.AltMaiorEsp, strRadarTx.AltMenorEsp, strRadarTx.Vt, strRadarTx.PRF, strRadarTx.NorthOffset);

% Decimar posiciones del transmisor para acelerar optimización
decimationFactor = 200;
PxT = PxT_(1:decimationFactor:end);
PyT = PyT_(1:decimationFactor:end);
PzT = PzT_(1:decimationFactor:end);
t_  = t(1:decimationFactor:end);

Tx_pos = [PxT; PyT; PzT];

%% Grid de targets subsuperficiales
gridTarget = strTarget.grid;
tg = createGridTarget(gridTarget.xSize, gridTarget.ySize, gridTarget.zMin, gridTarget.zMax, ...
    gridTarget.nx, gridTarget.ny, gridTarget.nz).';

target_center_2d = mean(tg(1:2,:), 2);   % centroide horizontal del volumen [2x1]

%% Parámetros del sistema
Pt     = strRadarTx.PotenciaTx;
lambda = physconst('lightspeed') / strRadarTx.FreqPortadora;   % 0.75 m @ 400 MHz
sigma  = strTarget.rcs;
n1     = 1;
n2     = strSystem.IndiceRefracaoSolo;

% Ancho de banda del chirp [Hz]
B = strRadarTx.FreqMayor - strRadarTx.FreqMenor;

% Geometría de la hélice cónica (para el cálculo de resolución vertical)
Delta_z    = strRadarTx.AltMaiorEsp  - strRadarTx.AltMenorEsp;    % variación de altura [m]
Delta_rho  = strRadarTx.RaioMaiorEsp - strRadarTx.RaioMenorEsp;   % variación de radio [m]
B_helix    = sqrt(Delta_z^2 + Delta_rho^2);                        % longitud arco espiral [m]
beta_helix = atan2(Delta_z, Delta_rho);                            % ángulo inclinación [rad]

fprintf('Longitud espiral B_helix = %.2f m,  beta = %.1f°\n', B_helix, beta_helix*180/pi);
fprintf('Ancho de banda B = %.0f MHz,  lambda = %.3f m\n', B/1e6, lambda);

%% Patrón de antena
radPattern = createRadiationPattern(strRadarTx.AperturaElev, strRadarTx.AperturaAzimut);

%% Parámetro de tradeoff potencia–resolución
% alpha_res = 0   → solo potencia (comportamiento original)
% alpha_res = 1   → solo resolución
% alpha_res = 0.3 → 70% peso en potencia, 30% en resolución (valor por defecto)
alpha_res = 0.3;

%% Ángulo de Brewster
angulo_brewster = atan(n2 / n1);
fprintf('\nÁngulo de Brewster: %.1f°\n', angulo_brewster*180/pi);
fprintf('Optimizando %d posiciones Tx  |  %d targets  |  alpha_res = %.2f\n', ...
    length(PxT), size(tg,2), alpha_res);

%% Configuración del optimizador
lb = [-strSystem.searchRadius, -strSystem.searchRadius];
ub = [ strSystem.searchRadius,  strSystem.searchRadius];
n_starts = 5;

options = optimoptions('fmincon', 'Display', 'none', ...
    'MaxIterations', 200, 'FunctionTolerance', 1e-10, ...
    'Algorithm', 'interior-point');

%% Arrays de resultados
Rx_opt_all     = zeros(3, length(PxT));
max_power_all  = zeros(1, length(PxT));
delta_xy_all   = zeros(1, length(PxT));
delta_z_all    = zeros(1, length(PxT));
cost_opt_all   = zeros(1, length(PxT));

%% Bucle de optimización por posición del transmisor
for tx_idx = 1:length(PxT)

    fprintf('\nTx %d/%d: [%.1f, %.1f, %.1f]...\n', tx_idx, length(PxT), PxT(tx_idx), PyT(tx_idx), PzT(tx_idx));

    Tx_pos_current = [PxT(tx_idx); PyT(tx_idx); PzT(tx_idx)];
    Rx_z           = PzT(tx_idx);

    % Ganancias del transmisor hacia todos los targets (precalculadas, fijas para este Tx)
    [~, ~, Gt_current] = calculateTxGainsForTargets(Tx_pos_current, tg, radPattern);

    % Geometría de búsqueda guiada por el ángulo de Brewster
    search_geo = buildRxSearchGeometry(Tx_pos_current, Rx_z, tg, angulo_brewster);

    results = cell(n_starts, 1);
    fvals   = zeros(n_starts, 1);

    for i = 1:n_starts
        x0 = chooseCandidateRxPos(i, search_geo.target_center, search_geo.targets_distance, ...
            search_geo.opposite_dir, search_geo.perp_dir, search_geo.distancia_brewster, search_geo.azimut_opuesto);
        x0 = max(lb, min(ub, x0));

        [x_opt, fval] = fmincon(...
            @(xy) objective_function(xy, tg, Tx_pos_current, Rx_z, radPattern, ...
                Gt_current, Pt, sigma, lambda, n1, n2, ...
                alpha_res, B, B_helix, beta_helix, target_center_2d), ...
            x0, [], [], [], [], lb, ub, [], options);

        results{i} = x_opt;
        fvals(i)   = fval;
    end

    % Mejor resultado (menor costo combinado)
    [~, best_idx]       = min(fvals);
    best_xy             = results{best_idx};
    Rx_opt_all(:,tx_idx) = [best_xy(1); best_xy(2); Rx_z];
    cost_opt_all(tx_idx) = fvals(best_idx);

    % Evaluar potencia y resolución en el punto óptimo (para reporting)
    c_power = objective_function(best_xy, tg, Tx_pos_current, Rx_z, radPattern, ...
        Gt_current, Pt, sigma, lambda, n1, n2, ...
        0, B, B_helix, beta_helix, target_center_2d);          % alpha=0 → cost = -log10(Pr)
    max_power_all(tx_idx) = 10^(-c_power);

    [delta_xy_all(tx_idx), delta_z_all(tx_idx)] = calculateBistaticResolution(...
        Tx_pos_current, Rx_opt_all(:,tx_idx), target_center_2d, n2, lambda, B, B_helix, beta_helix);

end

%% Resumen de resultados
fprintf('\n=== RESULTADOS DE OPTIMIZACIÓN (alpha_res = %.2f) ===\n', alpha_res);
fprintf('%-5s  %-20s  %-20s  %-12s  %-12s  %-12s\n', ...
    'Tx', 'Rx_opt [m]', 'Pot [dBm]', 'δxy [cm]', 'δz [cm]', 'Costo J');
for tx_idx = 1:length(PxT)
    fprintf('Tx%2d  [%6.1f,%6.1f,%5.1f]  %8.1f dBm  %8.2f cm  %8.2f cm  %8.4f\n', ...
        tx_idx, ...
        Rx_opt_all(1,tx_idx), Rx_opt_all(2,tx_idx), Rx_opt_all(3,tx_idx), ...
        10*log10(max_power_all(tx_idx)*1000), ...
        delta_xy_all(tx_idx)*100, ...
        delta_z_all(tx_idx)*100, ...
        cost_opt_all(tx_idx));
end

%% Guardar resultados
fprintf('\n=== GUARDANDO RESULTADOS ===\n');

outputDir = fullfile('io', 'plan_vuelo');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

timestamp     = datestr(now, 'yyyymmdd_HHMMSS');
filename_base = fullfile(outputDir, sprintf('optimizacion_%s', timestamp));

% Estructura completa
Resultados.Tx_positions   = Tx_pos;
Resultados.targets        = tg;
Resultados.target_center  = target_center_2d;
Resultados.parametros.n1        = n1;
Resultados.parametros.n2        = n2;
Resultados.parametros.Pt        = Pt;
Resultados.parametros.lambda    = lambda;
Resultados.parametros.B         = B;
Resultados.parametros.B_helix   = B_helix;
Resultados.parametros.beta_helix = beta_helix;
Resultados.parametros.alpha_res = alpha_res;
Resultados.parametros.sigma     = sigma;
Resultados.optimizacion.Rx_optimo          = Rx_opt_all;
Resultados.optimizacion.potencia_maxima_W  = max_power_all;
Resultados.optimizacion.potencia_maxima_dBm = 10*log10(max_power_all*1000);
Resultados.optimizacion.delta_xy_m         = delta_xy_all;
Resultados.optimizacion.delta_z_m          = delta_z_all;
Resultados.optimizacion.costo_combinado    = cost_opt_all;
Resultados.fecha_calculo = datestr(now);

save([filename_base '.mat'], 'Resultados');
fprintf('Datos .mat guardados en: %s.mat\n', filename_base);

% CSV con potencia y resolución
fid = fopen([filename_base '.csv'], 'w');
fprintf(fid, 'Tx_ID,Tx_X,Tx_Y,Tx_Z,Rx_X,Rx_Y,Rx_Z,Pot_W,Pot_dBm,delta_xy_m,delta_z_m,costo_J\n');
for tx_idx = 1:length(PxT)
    fprintf(fid, '%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.6e,%.3f,%.6f,%.6f,%.6f\n', ...
        tx_idx, ...
        Tx_pos(1,tx_idx), Tx_pos(2,tx_idx), Tx_pos(3,tx_idx), ...
        Rx_opt_all(1,tx_idx), Rx_opt_all(2,tx_idx), Rx_opt_all(3,tx_idx), ...
        max_power_all(tx_idx), 10*log10(max_power_all(tx_idx)*1000), ...
        delta_xy_all(tx_idx), delta_z_all(tx_idx), ...
        cost_opt_all(tx_idx));
end
fclose(fid);
fprintf('CSV guardado en: %s.csv\n', filename_base);

%% Visualización
plot_bistatic_configuration(PxT, PyT, PzT, tg, Rx_opt_all, n1, n2, angulo_brewster);

%% Interpolar trayectoria completa del receptor (PCHIP sobre puntos decimados)
R_opt_x = interp1(t_, Rx_opt_all(1,:), t, 'pchip', 'extrap');
R_opt_y = interp1(t_, Rx_opt_all(2,:), t, 'pchip', 'extrap');
R_opt_z = interp1(t_, Rx_opt_all(3,:), t, 'pchip', 'extrap');

R_opt_all_extrap = [R_opt_x; R_opt_y; R_opt_z];
Tx_all           = [PxT_; PyT_; PzT_];

save(fullfile(outputDir, 'R_opt_all_extrap.mat'), 'R_opt_all_extrap');
save(fullfile(outputDir, 'Tx_all.mat'),           'Tx_all');
fprintf('Trayectorias guardadas en: %s\n', outputDir);
