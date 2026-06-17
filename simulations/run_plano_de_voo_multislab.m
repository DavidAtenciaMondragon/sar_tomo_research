clc;
clear;
close all;

addpath(genpath('gs'))
addpath(genpath('proc'))
addpath(genpath('tools'))
addpath(genpath('flightpath'))
addpath(genpath('common'))

%% Cargar parámetros
systemJSON = json2struct(strcat('parametros',filesep,'system_multislab_plano_voo.json'));
strSystem  = systemJSON.system; clear systemJSON;

radarJSON  = json2struct(strcat('parametros',filesep,'radarTx_multislab_plano_voo.json'));
strRadarTx = radarJSON.radar; clear radarJSON;

targetJSON = json2struct(strcat('parametros',filesep,'target_multislab_plano_voo.json'));
strTarget  = targetJSON.target; clear targetJSON;

%% Estructura de capas de suelo
% n_layers(i), d_layers(i) : indice y espesor de la capa i (superficie → fondo)
n_layers = [strSystem.CapasSuelo.n];
d_layers = [strSystem.CapasSuelo.d];
M        = numel(n_layers);

fprintf('Medio multicapa (%d capas):\n', M);
cumD = 0;
for k = 1:M
    fprintf('  Capa %d: n = %.1f,  d = %.2f m  (z = %.2f .. %.2f m)\n', ...
        k, n_layers(k), d_layers(k), -cumD, -(cumD+d_layers(k)));
    cumD = cumD + d_layers(k);
end

%% Trayectoria del transmisor
[PxT_, PyT_, PzT_, t] = funcao_espiral(strRadarTx.NumVoltasEsp, strRadarTx.RaioMenorEsp, ...
    strRadarTx.RaioMaiorEsp, strRadarTx.AltMaiorEsp, strRadarTx.AltMenorEsp, ...
    strRadarTx.Vt, strRadarTx.PRF, strRadarTx.NorthOffset);

% Decimar posiciones del transmisor para acelerar la optimización
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

target_center_2d = mean(tg(1:2,:), 2);   % centroide horizontal [2x1]
target_center_z  = mean(tg(3,:));         % profundidad media

% Capa que contiene el centroide del volumen de targets
depthCentroid = -target_center_z;
cumDepths     = cumsum(d_layers);
L_centroid    = find(cumDepths >= depthCentroid, 1, 'first');
n_eff         = n_layers(L_centroid);     % indice efectivo del target (para resolución)

fprintf('\nCentroide de targets: z = %.2f m → capa %d (n = %.1f)\n', ...
    target_center_z, L_centroid, n_eff);

%% Parámetros del sistema
f      = strRadarTx.FreqPortadora;
lambda = physconst('lightspeed') / f;
Pt     = strRadarTx.PotenciaTx;
sigma  = strTarget.rcs;
B      = strRadarTx.FreqMayor - strRadarTx.FreqMenor;

Delta_z   = strRadarTx.AltMaiorEsp  - strRadarTx.AltMenorEsp;
Delta_rho = strRadarTx.RaioMaiorEsp - strRadarTx.RaioMenorEsp;
B_helix   = sqrt(Delta_z^2 + Delta_rho^2);
beta_helix = atan2(Delta_z, Delta_rho);

fprintf('\nf = %.0f MHz,  lambda = %.3f m\n', f/1e6, lambda);
fprintf('B = %.0f MHz,  B_helix = %.2f m,  beta_helix = %.1f°\n', ...
    B/1e6, B_helix, beta_helix*180/pi);

%% Ángulo de incidencia óptimo (generalizado para multicapa)
% Para una interfaz única, el máximo de T_TM ocurre en el ángulo de Brewster
% theta_B = atan(n2/n1). Para una pila multicapa no existe solución analítica
% cerrada: se barre theta en [0, 90°) y se identifica el máximo numérico.
% Se usan los slabs que la onda atraviesa hasta la capa del centroide.

n_slabs_cent = n_layers(1:L_centroid-1);
d_slabs_cent = d_layers(1:L_centroid-1);

theta_sweep = linspace(0, pi/2 - 1e-4, 2000);
[~, ~, ~, T_sweep] = multislabTMcoef(f, theta_sweep, 1.0, n_slabs_cent, d_slabs_cent, n_layers(L_centroid));
[T_max, idx_opt]   = max(T_sweep);
angulo_opt         = theta_sweep(idx_opt);

% Ángulo de Brewster de la primera interfaz (referencia)
angulo_brewster_1 = atan(n_layers(1) / 1.0);

fprintf('\nÁngulo de Brewster (interfaz 1, aire→n=%.1f): %.1f°\n', ...
    n_layers(1), angulo_brewster_1*180/pi);
fprintf('Ángulo óptimo TM multicapa (máximo T_TM=%.4f):  %.1f°\n', ...
    T_max, angulo_opt*180/pi);

%% Patrón de antena
radPattern = createRadiationPattern(strRadarTx.AperturaElev, strRadarTx.AperturaAzimut);

%% Parámetro de tradeoff potencia–resolución
% alpha_res = 0  → solo potencia
% alpha_res = 1  → solo resolución
% alpha_res ∈ (0,1) → combinado; default: 0.3
alpha_res = 0.3;

%% Tolerancia del bisector (más laxa en optimización para reducir tiempo de cómputo)
threshold_opt = 1e-4;   % m — para el optimizador (fmincon)

%% Configuración del optimizador
lb = [-strSystem.searchRadius, -strSystem.searchRadius];
ub = [ strSystem.searchRadius,  strSystem.searchRadius];
n_starts = 5;

options = optimoptions('fmincon', 'Display', 'none', ...
    'MaxIterations', 200, 'FunctionTolerance', 1e-10, ...
    'Algorithm', 'interior-point');

%% Arrays de resultados
Rx_opt_all    = zeros(3, length(PxT));
max_power_all = zeros(1, length(PxT));
delta_xy_all  = zeros(1, length(PxT));
delta_z_all   = zeros(1, length(PxT));
cost_opt_all  = zeros(1, length(PxT));

fprintf('\nOptimizando %d posiciones Tx  |  %d targets  |  alpha_res = %.2f\n', ...
    length(PxT), size(tg,2), alpha_res);

%% Bucle de optimización por posición del transmisor
for tx_idx = 1:length(PxT)

    fprintf('\nTx %d/%d: [%.1f, %.1f, %.1f]...\n', ...
        tx_idx, length(PxT), PxT(tx_idx), PyT(tx_idx), PzT(tx_idx));

    Tx_pos_current = [PxT(tx_idx); PyT(tx_idx); PzT(tx_idx)];
    Rx_z           = PzT(tx_idx);

    % Ganancias del transmisor hacia todos los targets
    [~, ~, Gt_current] = calculateTxGainsForTargets(Tx_pos_current, tg, radPattern);

    % Datos Tx→Target (precalculados, constantes para este Tx)
    [R_T_pre, T1_pre] = precomputeTxDataMultislab( ...
        Tx_pos_current, tg, n_layers, d_layers, f, threshold_opt);

    % Geometría de búsqueda guiada por el ángulo óptimo multicapa
    search_geo = buildRxSearchGeometry(Tx_pos_current, Rx_z, tg, angulo_opt);

    results = cell(n_starts, 1);
    fvals   = zeros(n_starts, 1);

    for i = 1:n_starts
        x0 = chooseCandidateRxPos(i, search_geo.target_center, search_geo.targets_distance, ...
            search_geo.opposite_dir, search_geo.perp_dir, ...
            search_geo.distancia_brewster, search_geo.azimut_opuesto);
        x0 = max(lb, min(ub, x0));

        [x_opt, fval] = fmincon(...
            @(xy) objective_function_multislab(xy, tg, Tx_pos_current, Rx_z, radPattern, ...
                Gt_current, Pt, sigma, lambda, ...
                n_layers, d_layers, f, threshold_opt, ...
                alpha_res, B, B_helix, beta_helix, target_center_2d, ...
                n_eff, R_T_pre, T1_pre), ...
            x0, [], [], [], [], lb, ub, [], options);

        results{i} = x_opt;
        fvals(i)   = fval;
    end

    % Mejor resultado entre los n_starts arranques
    [~, best_idx]         = min(fvals);
    best_xy               = results{best_idx};
    Rx_opt_all(:, tx_idx) = [best_xy(1); best_xy(2); Rx_z];
    cost_opt_all(tx_idx)  = fvals(best_idx);

    % Evaluar potencia en el punto óptimo (alpha=0 → -log10(Pr))
    c_power = objective_function_multislab(best_xy, tg, Tx_pos_current, Rx_z, radPattern, ...
        Gt_current, Pt, sigma, lambda, ...
        n_layers, d_layers, f, threshold_opt, ...
        0, B, B_helix, beta_helix, target_center_2d, ...
        n_eff, R_T_pre, T1_pre);
    max_power_all(tx_idx) = 10^(-c_power);

    [delta_xy_all(tx_idx), delta_z_all(tx_idx)] = calculateBistaticResolution( ...
        Tx_pos_current, Rx_opt_all(:, tx_idx), target_center_2d, ...
        n_eff, lambda, B, B_helix, beta_helix);
end

%% Resumen de resultados
fprintf('\n=== RESULTADOS DE OPTIMIZACIÓN (alpha_res = %.2f) ===\n', alpha_res);
fprintf('Capas: ');
for k = 1:M
    fprintf('n%d=%.1f(d=%.1fm) ', k, n_layers(k), d_layers(k));
end
fprintf('\n');
fprintf('%-5s  %-22s  %-12s  %-10s  %-10s  %-10s\n', ...
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

outputDir = fullfile('io', 'plan_vuelo_multislab');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

timestamp     = datestr(now, 'yyyymmdd_HHMMSS');
filename_base = fullfile(outputDir, sprintf('optimizacion_multislab_%s', timestamp));

Resultados.Tx_positions   = Tx_pos;
Resultados.targets        = tg;
Resultados.target_center  = target_center_2d;
Resultados.parametros.n_layers    = n_layers;
Resultados.parametros.d_layers    = d_layers;
Resultados.parametros.n_eff       = n_eff;
Resultados.parametros.L_centroid  = L_centroid;
Resultados.parametros.f           = f;
Resultados.parametros.Pt          = Pt;
Resultados.parametros.lambda      = lambda;
Resultados.parametros.B           = B;
Resultados.parametros.B_helix     = B_helix;
Resultados.parametros.beta_helix  = beta_helix;
Resultados.parametros.alpha_res   = alpha_res;
Resultados.parametros.sigma       = sigma;
Resultados.parametros.angulo_opt  = angulo_opt;
Resultados.optimizacion.Rx_optimo           = Rx_opt_all;
Resultados.optimizacion.potencia_maxima_W   = max_power_all;
Resultados.optimizacion.potencia_maxima_dBm = 10*log10(max_power_all*1000);
Resultados.optimizacion.delta_xy_m          = delta_xy_all;
Resultados.optimizacion.delta_z_m           = delta_z_all;
Resultados.optimizacion.costo_combinado     = cost_opt_all;
Resultados.fecha_calculo = datestr(now);

save([filename_base '.mat'], 'Resultados');
fprintf('Datos .mat guardados en: %s.mat\n', filename_base);

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
plot_bistatic_configuration(PxT, PyT, PzT, tg, Rx_opt_all, 1, n_eff, angulo_opt);

%% Interpolar trayectoria completa del receptor (PCHIP)
R_opt_x = interp1(t_, Rx_opt_all(1,:), t, 'pchip', 'extrap');
R_opt_y = interp1(t_, Rx_opt_all(2,:), t, 'pchip', 'extrap');
R_opt_z = interp1(t_, Rx_opt_all(3,:), t, 'pchip', 'extrap');

R_opt_all_extrap = [R_opt_x; R_opt_y; R_opt_z];
Tx_all           = [PxT_; PyT_; PzT_];

save(fullfile(outputDir, 'R_opt_all_extrap.mat'), 'R_opt_all_extrap');
save(fullfile(outputDir, 'Tx_all.mat'),           'Tx_all');
fprintf('Trayectorias guardadas en: %s\n', outputDir);
