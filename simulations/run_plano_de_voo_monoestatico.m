clc;
clear;
close all;
%
% RUN_PLANO_DE_VOO_MONOESTATICO  Línea base monoestática (Rx = Tx) para comparar
% contra la optimización bistática de run_plano_de_voo.m.
%
% Usa la misma trayectoria helicoidal, los mismos parámetros JSON y las mismas
% funciones de potencia/resolución (objective_function.m, calculateBistaticResolution.m)
% que run_plano_de_voo.m. No hay optimización: para cada posición del Tx a lo
% largo de la hélice, la única posición de Rx es la del propio Tx, así que se
% evalúa la función de costo una sola vez por posición (sin fmincon).
%
% El struct Resultados y el CSV de salida tienen exactamente los mismos campos
% y columnas que run_plano_de_voo.m, para poder compararlos directamente
% (p.ej. ganancia de potencia y resolución del biestático optimizado vs. este
% monoestático).

addpath(genpath('gs'))
addpath(genpath('proc'))
addpath(genpath('tools'))
addpath(genpath('flightpath'))

%% Cargar parámetros (mismos JSON que run_plano_de_voo.m, para comparabilidad directa)
systemJSON = json2struct(strcat('parametros',filesep,'system_espiral_plano_voo.json'));
strSystem  = systemJSON.system; clear systemJSON;

radarJSON  = json2struct(strcat('parametros',filesep,'radarTx_espiral_plano_voo.json'));
strRadarTx = radarJSON.radar; clear radarJSON;

targetJSON = json2struct(strcat('parametros',filesep,'target_espiral_plano_voo.json'));
strTarget  = targetJSON.target; clear targetJSON;

%% Trayectoria del transmisor (idéntica a run_plano_de_voo.m)
[PxT_, PyT_, PzT_, t] = funcao_espiral(strRadarTx.NumVoltasEsp, strRadarTx.RaioMenorEsp, strRadarTx.RaioMaiorEsp,...
                                    strRadarTx.AltMaiorEsp, strRadarTx.AltMenorEsp, strRadarTx.Vt, strRadarTx.PRF, strRadarTx.NorthOffset);

% Decimar posiciones del transmisor (mismo factor que run_plano_de_voo.m, para
% que la posición i de este script corresponda a la posición i del optimizado)
decimationFactor = 200;
PxT = PxT_(1:decimationFactor:end);
PyT = PyT_(1:decimationFactor:end);
PzT = PzT_(1:decimationFactor:end);
t_  = t(1:decimationFactor:end);

Tx_pos = [PxT; PyT; PzT];

%% Configuración monoestática: el receptor coincide con el transmisor en todo instante
Rx_opt_all = Tx_pos;

%% Grid de targets subsuperficiales
gridTarget = strTarget.grid;
tg = createGridTarget(gridTarget.xSize, gridTarget.ySize, gridTarget.zMin, gridTarget.zMax, ...
    gridTarget.nx, gridTarget.ny, gridTarget.nz).';

target_center_2d = mean(tg(1:2,:), 2);   % centroide horizontal del volumen [2x1]

%% Parámetros del sistema (idénticos a run_plano_de_voo.m)
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

%% Parámetro de tradeoff potencia–resolución (mismo valor que run_plano_de_voo.m,
%  para que el costo J sea comparable entre ambos escenarios)
alpha_res = 0.3;

%% Ángulo de Brewster (solo como referencia/anotación de la figura; en el caso
%  monoestático no guía ninguna búsqueda porque no hay optimización)
angulo_brewster = atan(n2 / n1);
fprintf('\nÁngulo de Brewster: %.1f°\n', angulo_brewster*180/pi);
fprintf('Evaluando %d posiciones Tx=Rx (monoestático)  |  %d targets  |  alpha_res = %.2f\n', ...
    length(PxT), size(tg,2), alpha_res);

%% Arrays de resultados (mismos nombres que run_plano_de_voo.m)
max_power_all  = zeros(1, length(PxT));
delta_xy_all   = zeros(1, length(PxT));
delta_z_all    = zeros(1, length(PxT));
cost_opt_all   = zeros(1, length(PxT));

%% Evaluación directa por posición del transmisor
% Sin fmincon: Rx = Tx es la única posición posible, así que la función de
% costo se evalúa una sola vez por posición (en vez de 5 arranques x fmincon).
for tx_idx = 1:length(PxT)

    fprintf('\nTx %d/%d: [%.1f, %.1f, %.1f]...\n', tx_idx, length(PxT), PxT(tx_idx), PyT(tx_idx), PzT(tx_idx));

    Tx_pos_current = [PxT(tx_idx); PyT(tx_idx); PzT(tx_idx)];
    Rx_z           = PzT(tx_idx);
    xy_mono        = [PxT(tx_idx); PyT(tx_idx)];   % posición horizontal del Rx = la del Tx

    % Ganancias del transmisor hacia todos los targets
    [~, ~, Gt_current] = calculateTxGainsForTargets(Tx_pos_current, tg, radPattern);

    % Datos de refracción Tx→Target
    [R_T_pre, T1_pre] = precomputeTxData(Tx_pos_current, tg, n1, n2);

    % Costo combinado J en el punto monoestático (misma función y mismo alpha_res
    % que usa el optimizador bistático, para que J sea directamente comparable)
    cost_opt_all(tx_idx) = objective_function(xy_mono, tg, Tx_pos_current, Rx_z, radPattern, ...
        Gt_current, Pt, sigma, lambda, n1, n2, ...
        alpha_res, B, B_helix, beta_helix, target_center_2d, R_T_pre, T1_pre);

    % Potencia recibida pura (alpha=0 → cost = -log10(Pr)), igual que en run_plano_de_voo.m
    c_power = objective_function(xy_mono, tg, Tx_pos_current, Rx_z, radPattern, ...
        Gt_current, Pt, sigma, lambda, n1, n2, ...
        0, B, B_helix, beta_helix, target_center_2d, R_T_pre, T1_pre);
    max_power_all(tx_idx) = 10^(-c_power);

    [delta_xy_all(tx_idx), delta_z_all(tx_idx)] = calculateBistaticResolution(...
        Tx_pos_current, Rx_opt_all(:,tx_idx), target_center_2d, n2, lambda, B, B_helix, beta_helix);

end

%% Resumen de resultados
fprintf('\n=== RESULTADOS MONOESTÁTICOS (alpha_res = %.2f) ===\n', alpha_res);
fprintf('%-5s  %-20s  %-20s  %-12s  %-12s  %-12s\n', ...
    'Tx', 'Rx=Tx [m]', 'Pot [dBm]', 'δxy [cm]', 'δz [cm]', 'Costo J');
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

outputDir = fullfile('io', 'plan_vuelo_monoestatico');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

timestamp     = datestr(now, 'yyyymmdd_HHMMSS');
filename_base = fullfile(outputDir, sprintf('monoestatico_%s', timestamp));

% Estructura completa (mismos campos que run_plano_de_voo.m)
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

% CSV con potencia y resolución (mismas columnas que run_plano_de_voo.m)
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

%% Visualización (Rx se traza superpuesto al Tx: es el caso Rx=Tx)
plot_bistatic_configuration(PxT, PyT, PzT, tg, Rx_opt_all, n1, n2, angulo_brewster, outputDir);

%% Interpolar trayectoria completa del receptor (PCHIP sobre puntos decimados)
R_opt_x = interp1(t_, Rx_opt_all(1,:), t, 'pchip', 'extrap');
R_opt_y = interp1(t_, Rx_opt_all(2,:), t, 'pchip', 'extrap');
R_opt_z = interp1(t_, Rx_opt_all(3,:), t, 'pchip', 'extrap');

R_opt_all_extrap = [R_opt_x; R_opt_y; R_opt_z];
Tx_all           = [PxT_; PyT_; PzT_];

save(fullfile(outputDir, 'R_opt_all_extrap.mat'), 'R_opt_all_extrap');
save(fullfile(outputDir, 'Tx_all.mat'),           'Tx_all');
fprintf('Trayectorias guardadas en: %s\n', outputDir);
