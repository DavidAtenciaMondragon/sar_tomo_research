clc;
clear;
close all;
%
% OTIMIZA_TX_BREWSTER_PURO  Pipeline de dos etapas que optimiza *solo* la
% transmitancia TM de cada tramo (T1 en la Etapa 1, T2 en la Etapa 2) --
% NO la ecuación de radar completa (sin ganancia de antena, sin 1/R^2, sin
% peso de resolución). Es la contraparte "aislada" de run_otimiza_tx.m: sirve
% para ver dónde pone Tx y Rx un criterio de Brewster puro, sin que la
% pérdida de propagación domine la decisión (como se vio en la sesión
% anterior: con potencia completa, el óptimo real gana al candidato Brewster
% por +4.7 dB, casi todo por 1/R^2, no por T2).
%
% ETAPA 1 (Tx): igual que run_otimiza_tx.m con gamma=0 fijo -- minimiza
% -log10(T1_sim), T1 simulado con geometría real (Fermat + grid completo +
% varias posiciones de la hélice candidata), búsqueda de sección áurea sobre
% psi0.
%
% ETAPA 2 (Rx): para cada posición decimada del Tx, maximiza T2 promedio
% (sobre todo el grid de alvos) variando SOLO la distancia radial de Rx al
% centroide del grid, a lo largo de un azimut FIJO (opuesto al Tx, misma
% convención que buildRxSearchGeometry.m). Esto es deliberado, no una
% simplificación arbitraria: T2 depende únicamente del ángulo de incidencia
% en el suelo, que a su vez depende solo de la distancia radial Rx-alvo, NO
% del azimut -- así que el conjunto de posiciones que maximizan T2 es un
% CÍRCULO completo alrededor del alvo, no un punto único. Fijar el azimut
% (en vez de dejar que un fmincon 2D converja a un punto arbitrario de ese
% círculo, dependiente del arranque) es la única forma honesta de resolver
% esa degeneración.
%
% Potencia, delta_xy y delta_z SÍ se calculan y guardan, pero solo como
% diagnóstico posterior -- nunca son parte del costo que se optimiza aquí.
%
% Salida: Resultados.mat + CSVs en io/otimiza_tx_brewster_puro/, con la misma
% estructura de campos que run_otimiza_tx.m (Tx_positions, targets,
% target_center, optimizacion.*) para poder reusar
% verifica_otimiza_tx_vs_planvoo.m manualmente si se quiere comparar.

addpath(genpath('gs'))
addpath(genpath('tools'))
addpath(genpath('flightpath'))

%% Cargar parámetros
systemJSON = json2struct(strcat('parametros',filesep,'system_espiral_plano_voo.json'));
strSystem  = systemJSON.system; clear systemJSON;

radarJSON  = json2struct(strcat('parametros',filesep,'radarTx_espiral_plano_voo.json'));
strRadarTx = radarJSON.radar; clear radarJSON;

targetJSON = json2struct(strcat('parametros',filesep,'target_espiral_plano_voo.json'));
strTarget  = targetJSON.target; clear targetJSON;

optJSON    = json2struct(strcat('parametros',filesep,'otimiza_tx_brewster_puro.json'));
strOpt     = optJSON.otimiza_tx_brewster_puro; clear optJSON;

%% Grid de alvos subsuperficiales
gridTarget = strTarget.grid;
tg = createGridTarget(gridTarget.xSize, gridTarget.ySize, gridTarget.zMin, gridTarget.zMax, ...
    gridTarget.nx, gridTarget.ny, gridTarget.nz).';
target_center_2d = mean(tg(1:2,:), 2);

%% Parámetros físicos
n1     = 1;
n2     = strSystem.IndiceRefracaoSolo;
c      = physconst('lightspeed');
lambda = c / strRadarTx.FreqPortadora;
B      = strRadarTx.FreqMayor - strRadarTx.FreqMenor;

theta_B_air    = atan(n2 / n1);   % Brewster lado aire  (Tx->Tg, y Tg->Rx visto desde el aire)
theta_B_ground = atan(n1 / n2);   % Brewster lado suelo (Tg->Rx visto desde el suelo)

fprintf('theta_B,aire = %.2f°   theta_B,suelo = %.2f°\n', theta_B_air*180/pi, theta_B_ground*180/pi);
fprintf('n1=%.1f, n2=%.1f, f0=%.0f MHz, B=%.0f MHz, lambda=%.3f m\n', ...
    n1, n2, strRadarTx.FreqPortadora/1e6, B/1e6, lambda);

%% Geometría de referencia (hélice actualmente en uso)
rho0_nom      = (strRadarTx.RaioMenorEsp + strRadarTx.RaioMaiorEsp) / 2;
z0_nom        = (strRadarTx.AltMenorEsp  + strRadarTx.AltMaiorEsp)  / 2;
psi0_nom      = atan(rho0_nom / z0_nom);
R0_nom        = hypot(rho0_nom, z0_nom);
Delta_rho_nom = strRadarTx.RaioMaiorEsp - strRadarTx.RaioMenorEsp;
Delta_z_nom   = strRadarTx.AltMaiorEsp  - strRadarTx.AltMenorEsp;
B_helix       = hypot(Delta_rho_nom, Delta_z_nom);

fprintf('\nSistema actual: psi0=%.2f°, R0=%.2f m, z0=%.2f m, B_helix=%.2f m\n', ...
    psi0_nom*180/pi, R0_nom, z0_nom, B_helix);

%% ============================================================
%% ETAPA 1 (Tx): maximizar SOLO T1 -- varredura + óptimo por sección áurea
%% ============================================================
psi0_sweep_deg = strOpt.psiMinDeg:strOpt.psiStepDeg:strOpt.psiMaxDeg;
Npsi      = numel(psi0_sweep_deg);
scenarios = strOpt.scenarios;
Nscen     = numel(scenarios);
Nsamples  = strOpt.azimuthSamples;

T1_sim = nan(Nscen, Npsi);

fprintf('\n=== ETAPA 1: varredura de psi0 (solo T1, geometria real) ===\n');
for s = 1:Nscen
    scenario = scenarios{s};
    fprintf('-- Escenario: %s --\n', scenario);
    for k = 1:Npsi
        T1_sim(s,k) = simulateT1Design(psi0_sweep_deg(k), scenario, R0_nom, z0_nom, B_helix, ...
            strRadarTx, tg, n1, n2, Nsamples);
    end
end

tolDeg = strOpt.goldenSectionTolDeg;
psi0_star_sim = nan(1, Nscen);
T1_star_sim   = nan(1, Nscen);

fprintf('\n=== ETAPA 1: óptimo psi0* (maximiza T1, sección áurea) ===\n');
for s = 1:Nscen
    scenario = scenarios{s};
    costFun = @(psi0_deg) costT1Sim(psi0_deg, scenario, R0_nom, z0_nom, B_helix, ...
        strRadarTx, tg, n1, n2, Nsamples);
    psi0_star_sim(s) = goldenSectionMin(costFun, strOpt.psiMinDeg, strOpt.psiMaxDeg, tolDeg);
    T1_star_sim(s)   = simulateT1Design(psi0_star_sim(s), scenario, R0_nom, z0_nom, B_helix, ...
        strRadarTx, tg, n1, n2, Nsamples);
    fprintf('  %-10s  psi0*_sim = %6.2f°   T1=%.4f (%.3f dB)   |diff a theta_B,aire| = %.2f°\n', ...
        scenario, psi0_star_sim(s), T1_star_sim(s), 10*log10(T1_star_sim(s)), ...
        abs(psi0_star_sim(s) - theta_B_air*180/pi));
end

%% ============================================================
%% ETAPA 1 -> ETAPA 2: hélice final del Tx
%% ============================================================
scenarioFinal = strOpt.scenarioFinal;
idxFinal      = find(strcmp(scenarios, scenarioFinal), 1);
if isempty(idxFinal)
    error('otimiza_tx_brewster_puro:scenarioFinal', ...
        'strOpt.scenarioFinal="%s" no está en strOpt.scenarios.', scenarioFinal);
end

psi0_star_final = psi0_star_sim(idxFinal);
psi0_star_rad   = deg2rad(psi0_star_final);

switch scenarioFinal
    case 'R0fixed'
        rho0_star = R0_nom * sin(psi0_star_rad);
        z0_star   = R0_nom * cos(psi0_star_rad);
    case 'z0fixed'
        z0_star   = z0_nom;
        rho0_star = z0_nom * tan(psi0_star_rad);
end

Delta_rho_star = B_helix * cos(psi0_star_rad);
Delta_z_star   = B_helix * sin(psi0_star_rad);

RaioMenor_star = rho0_star - Delta_rho_star/2;
RaioMaior_star = rho0_star + Delta_rho_star/2;
AltMenor_star  = z0_star   - Delta_z_star/2;
AltMaior_star  = z0_star   + Delta_z_star/2;

fprintf('\n=== ETAPA 1 -> ETAPA 2: hélice final del Tx (escenario %s) ===\n', scenarioFinal);
fprintf('psi0* = %.2f°   rho0=%.2f m   z0=%.2f m\n', psi0_star_final, rho0_star, z0_star);
fprintf('RaioMenor=%.2f m  RaioMaior=%.2f m  AltMenor=%.2f m  AltMaior=%.2f m\n', ...
    RaioMenor_star, RaioMaior_star, AltMenor_star, AltMaior_star);

if RaioMenor_star <= 1.0 || AltMenor_star <= 1.0
    error('otimiza_tx_brewster_puro:geometriaInvalida', ...
        'La hélice final es geométricamente inválida (radio o altura mínima <= 0).');
end

[PxT_, PyT_, PzT_, t] = funcao_espiral(strRadarTx.NumVoltasEsp, RaioMaior_star, RaioMenor_star, ...
    AltMaior_star, AltMenor_star, strRadarTx.Vt, strRadarTx.PRF, strRadarTx.NorthOffset);

decimationFactor = strOpt.decimationFactor;
PxT = PxT_(1:decimationFactor:end);
PyT = PyT_(1:decimationFactor:end);
PzT = PzT_(1:decimationFactor:end);
t_  = t(1:decimationFactor:end);
Tx_pos = [PxT; PyT; PzT];
N = length(PxT);

%% ============================================================
%% ETAPA 2 (Rx): maximizar SOLO T2 promedio, distancia radial a lo largo de
%% un azimut fijo (opuesto al Tx) -- el resto (potencia, delta_xy, delta_z)
%% se calcula solo como diagnóstico, no entra en el costo.
%% ============================================================
Pt    = strRadarTx.PotenciaTx;
sigma = strTarget.rcs;
radPattern = createRadiationPattern(strRadarTx.AperturaElev, strRadarTx.AperturaAzimut);
beta_helix = psi0_star_rad;

fprintf('\n=== ETAPA 2: %d posiciones Tx (decimationFactor=%d), maximizando T2 puro ===\n', ...
    N, decimationFactor);

target_center_3d = mean(tg, 2);   % centroide 3D, para el diagnóstico de ángulo en el suelo

Rx_opt_all     = zeros(3, N);
T1_all         = zeros(1, N);
T2_all         = zeros(1, N);
psi_Rx_air_all = zeros(1, N);   % parámetro de búsqueda: ángulo AIRE de colocación de Rx (recta Rx->centroide horizontal)
ang2_ground_all = zeros(1, N);  % diagnóstico: ángulo de incidencia real en el SUELO (vía Fermat), el que compara con theta_B_ground
max_power_all  = zeros(1, N);   % diagnóstico
delta_xy_all   = zeros(1, N);   % diagnóstico
delta_z_all    = zeros(1, N);   % diagnóstico
cost_opt_all   = zeros(1, N);   % -log10(T2_avg) -- el costo que sí se optimizó

for tx_idx = 1:N
    Tx_pos_current = [PxT(tx_idx); PyT(tx_idx); PzT(tx_idx)];
    Rx_z           = PzT(tx_idx);

    search_geo = buildRxSearchGeometry(Tx_pos_current, Rx_z, tg, theta_B_air);
    azimuth    = search_geo.azimut_opuesto;   % dirección fija, opuesta al Tx

    % psiRx_deg parametriza la COLOCACIÓN de Rx (ángulo, medido en el aire,
    % de la recta Rx->centroide horizontal) -- NO es el ángulo de incidencia
    % en el suelo. El óptimo de este parámetro cae cerca de theta_B,aire
    % (63.4°) porque esa es, por Snell, la colocación en el aire que produce
    % el ángulo de incidencia real en el suelo (vía refracción de Fermat)
    % más cercano a theta_B,suelo (26.6°) -- son la MISMA geometría de rayo
    % de Brewster, vista desde cada lado.
    costFun = @(psiRx_deg) costT2Sim(psiRx_deg, azimuth, target_center_2d, Rx_z, tg, n1, n2);
    psiRx_star_deg = goldenSectionMin(costFun, strOpt.psiMinDeg, strOpt.psiMaxDeg, tolDeg);
    psi_Rx_air_all(tx_idx) = psiRx_star_deg;

    rho_Rx = Rx_z * tan(deg2rad(psiRx_star_deg));
    Rx_xy  = target_center_2d(:)' + rho_Rx * [cos(azimuth), sin(azimuth)];
    Rx_opt_all(:, tx_idx) = [Rx_xy(1); Rx_xy(2); Rx_z];
    Rx_pos_current = Rx_opt_all(:, tx_idx);

    % -- Diagnóstico (NO optimizado): T1, T2, potencia completa, resolución,
    %    y el ángulo real de incidencia en el suelo (vía Fermat, en el
    %    centroide) para comparar directamente contra theta_B,suelo --
    [~, ~, Gt_current]  = calculateTxGainsForTargets(Tx_pos_current, tg, radPattern);
    [R_T_pre, T1_pre]   = precomputeTxData(Tx_pos_current, tg, n1, n2);
    T1_all(tx_idx)      = mean(T1_pre);
    T2_all(tx_idx)       = 10^(-costFun(psiRx_star_deg));   % T2 promedio en el óptimo
    cost_opt_all(tx_idx) = costFun(psiRx_star_deg);

    intersec2_c = calculateRefractionPointFermat(target_center_3d, Rx_pos_current, n2, n1);
    ang2_ground_all(tx_idx) = calculateIncidenceAngle(target_center_3d, intersec2_c, Rx_pos_current);

    c_power = objective_function([Rx_xy(1); Rx_xy(2)], tg, Tx_pos_current, Rx_z, radPattern, ...
        Gt_current, Pt, sigma, lambda, n1, n2, 0, B, B_helix, beta_helix, target_center_2d, R_T_pre, T1_pre);
    max_power_all(tx_idx) = 10^(-c_power);

    [delta_xy_all(tx_idx), delta_z_all(tx_idx)] = calculateBistaticResolution(...
        Tx_pos_current, Rx_opt_all(:,tx_idx), target_center_2d, n2, lambda, B, B_helix, beta_helix);

    if mod(tx_idx, 20) == 0 || tx_idx == N
        fprintf('  Tx %3d/%3d: psi_Rx,aire*=%5.1f°  ang2,suelo=%5.1f°  T1=%.4f  T2=%.4f  Pot=%6.1f dBm  dz=%6.1f cm\n', ...
            tx_idx, N, psiRx_star_deg, ang2_ground_all(tx_idx)*180/pi, T1_all(tx_idx), T2_all(tx_idx), ...
            10*log10(max_power_all(tx_idx)*1000), delta_z_all(tx_idx)*100);
    end
end

fprintf('\n=== RESUMEN ETAPA 2 (Brewster puro) ===\n');
fprintf('psi_Rx,aire* medio  = %6.2f°  (theta_B,aire=%.2f°,   |diff| medio=%.2f°)\n', ...
    mean(psi_Rx_air_all), theta_B_air*180/pi, mean(abs(psi_Rx_air_all - theta_B_air*180/pi)));
fprintf('ang2,suelo medio    = %6.2f°  (theta_B,suelo=%.2f°,  |diff| medio=%.2f°)\n', ...
    mean(ang2_ground_all)*180/pi, theta_B_ground*180/pi, mean(abs(ang2_ground_all*180/pi - theta_B_ground*180/pi)));
fprintf('T1 medio = %.4f   T2 medio = %.4f\n', mean(T1_all), mean(T2_all));
fprintf('Potencia media (diagnóstico) = %.2f dBm\n', 10*log10(mean(max_power_all)*1000));
fprintf('delta_z medio (diagnóstico)  = %.2f cm\n', mean(delta_z_all)*100);

%% Guardar resultados
outputDir = fullfile('io', 'otimiza_tx_brewster_puro');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
timestamp     = datestr(now, 'yyyymmdd_HHMMSS');
filename_base = fullfile(outputDir, sprintf('otimiza_tx_brewster_puro_%s', timestamp));

Resultados.parametros.n1        = n1;
Resultados.parametros.n2        = n2;
Resultados.parametros.lambda    = lambda;
Resultados.parametros.B         = B;
Resultados.parametros.B_helix   = B_helix;
Resultados.parametros.theta_B_air    = theta_B_air;
Resultados.parametros.theta_B_ground = theta_B_ground;
Resultados.parametros.psi0_nom  = psi0_nom;
Resultados.parametros.R0_nom    = R0_nom;
Resultados.parametros.z0_nom    = z0_nom;
Resultados.parametros.Pt        = Pt;
Resultados.parametros.sigma     = sigma;
Resultados.parametros.alpha_res = 0;   % no se pesa resolución -- solo diagnóstico

Resultados.etapa1.psi0_deg_sweep = psi0_sweep_deg;
Resultados.etapa1.T1_sim         = T1_sim;
Resultados.etapa1.scenarios      = {scenarios{:}};
Resultados.etapa1.psi0_star_sim  = psi0_star_sim;
Resultados.etapa1.T1_star_sim    = T1_star_sim;

Resultados.etapa1_final.scenario   = scenarioFinal;
Resultados.etapa1_final.psi0_star  = psi0_star_final;
Resultados.etapa1_final.rho0_star  = rho0_star;
Resultados.etapa1_final.z0_star    = z0_star;
Resultados.etapa1_final.beta_helix = beta_helix;
Resultados.etapa1_final.RaioMenor  = RaioMenor_star;
Resultados.etapa1_final.RaioMaior  = RaioMaior_star;
Resultados.etapa1_final.AltMenor   = AltMenor_star;
Resultados.etapa1_final.AltMaior   = AltMaior_star;

Resultados.Tx_positions  = Tx_pos;
Resultados.targets       = tg;
Resultados.target_center = target_center_2d;
Resultados.optimizacion.Rx_optimo           = Rx_opt_all;
Resultados.optimizacion.psi_Rx_aire_deg     = psi_Rx_air_all;    % parámetro de colocación (aire)
Resultados.optimizacion.ang2_suelo_deg      = ang2_ground_all*180/pi;  % diagnóstico: ángulo real (suelo, Fermat)
Resultados.optimizacion.T1                  = T1_all;
Resultados.optimizacion.T2                  = T2_all;
Resultados.optimizacion.potencia_maxima_W   = max_power_all;      % diagnóstico
Resultados.optimizacion.potencia_maxima_dBm = 10*log10(max_power_all*1000);  % diagnóstico
Resultados.optimizacion.delta_xy_m          = delta_xy_all;       % diagnóstico
Resultados.optimizacion.delta_z_m           = delta_z_all;        % diagnóstico
Resultados.optimizacion.costo_combinado     = cost_opt_all;       % -log10(T2), lo que sí se optimizó

Resultados.fecha_calculo = datestr(now);

save([filename_base '.mat'], 'Resultados');
fprintf('\nDatos .mat guardados en: %s.mat\n', filename_base);

fid = fopen([filename_base '.csv'], 'w');
fprintf(fid, ['Tiempo_s,Tx_X,Tx_Y,Tx_Z,Rx_X,Rx_Y,Rx_Z,psi_Rx_aire_deg,ang2_suelo_deg,' ...
    'T1,T2,Pot_dBm_diag,delta_xy_m_diag,delta_z_m_diag\n']);
for tx_idx = 1:N
    fprintf(fid, '%.6f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.4f,%.4f,%.6f,%.6f,%.3f,%.6f,%.6f\n', ...
        t_(tx_idx), Tx_pos(1,tx_idx), Tx_pos(2,tx_idx), Tx_pos(3,tx_idx), ...
        Rx_opt_all(1,tx_idx), Rx_opt_all(2,tx_idx), Rx_opt_all(3,tx_idx), ...
        psi_Rx_air_all(tx_idx), ang2_ground_all(tx_idx)*180/pi, T1_all(tx_idx), T2_all(tx_idx), ...
        10*log10(max_power_all(tx_idx)*1000), delta_xy_all(tx_idx), delta_z_all(tx_idx));
end
fclose(fid);
fprintf('CSV guardado en: %s.csv\n', filename_base);

%% Visualización
fig = figure('Position', [100, 100, 900, 500]);
plot(psi0_sweep_deg, T1_sim(idxFinal,:), '-', 'LineWidth', 2, 'Color', [0.17 0.48 0.71]);
hold on;
xline(theta_B_air*180/pi, ':', 'Color', [0.88 0.51 0.08], 'LineWidth', 1.5, 'DisplayName', '\theta_{B,aire}');
xline(psi0_star_final, '-.', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.5, 'DisplayName', '\psi_0^*');
xlabel('\psi_0 [deg]'); ylabel('T_1');
title(sprintf('Etapa 1 (Tx, Brewster puro) -- escenario %s', scenarioFinal), 'Interpreter', 'none');
legend('T_1(\psi_0)', '\theta_{B,aire}', '\psi_0^*', 'Location', 'southwest');
grid on; hold off;
fig_path = fullfile(outputDir, sprintf('etapa1_T1_%s.jpg', timestamp));
saveas(fig, fig_path);
fprintf('Figura (Etapa 1) guardada en: %s\n', fig_path);

plot_bistatic_configuration(PxT, PyT, PzT, tg, Rx_opt_all, n1, n2, theta_B_air, outputDir);

fprintf('\nListo. psi0*=%.2f° (Tx), ang2_suelo medio=%.2f° (Rx, vs theta_B,suelo=%.2f°) -- ambos\n', ...
    psi0_star_final, mean(ang2_ground_all)*180/pi, theta_B_ground*180/pi);
fprintf('transmitancia pura, sin considerar ganancia de antena ni 1/R^2.\n');

%% ==================== Funciones locales ====================

function T1_avg = simulateT1Design(psi0_deg, scenario, R0_nom, z0_nom, B_helix, strRadarTx, ...
    tg, n1, n2, Nsamples)
% Igual que simulateTxDesign de run_otimiza_tx.m, pero solo T1 (sin delta_z).

psi0_rad = deg2rad(psi0_deg);

switch scenario
    case 'R0fixed'
        rho0 = R0_nom * sin(psi0_rad);
        z0   = R0_nom * cos(psi0_rad);
    case 'z0fixed'
        z0   = z0_nom;
        rho0 = z0_nom * tan(psi0_rad);
    otherwise
        error('simulateT1Design:escenario', 'Escenario desconocido: %s', scenario);
end

Delta_rho = B_helix * cos(psi0_rad);
Delta_z   = B_helix * sin(psi0_rad);

RaioMenor = rho0 - Delta_rho/2;
RaioMaior = rho0 + Delta_rho/2;
AltMenor  = z0   - Delta_z/2;
AltMaior  = z0   + Delta_z/2;

if RaioMenor <= 1.0 || AltMenor <= 1.0
    T1_avg = NaN;
    return;
end

[Px_, Py_, Pz_, t] = funcao_espiral(strRadarTx.NumVoltasEsp, RaioMaior, RaioMenor, ...
    AltMaior, AltMenor, strRadarTx.Vt, strRadarTx.PRF, strRadarTx.NorthOffset);

idx = unique(round(linspace(1, numel(t), Nsamples)));
T1_pts = zeros(1, numel(idx));
for j = 1:numel(idx)
    Tx_pos = [Px_(idx(j)); Py_(idx(j)); Pz_(idx(j))];
    [~, T1_pre] = precomputeTxData(Tx_pos, tg, n1, n2);
    T1_pts(j) = mean(T1_pre);
end
T1_avg = mean(T1_pts);

end

function J = costT1Sim(psi0_deg, scenario, R0_nom, z0_nom, B_helix, strRadarTx, tg, n1, n2, Nsamples)
T1_avg = simulateT1Design(psi0_deg, scenario, R0_nom, z0_nom, B_helix, strRadarTx, tg, n1, n2, Nsamples);
if ~isfinite(T1_avg) || T1_avg <= 0
    J = 1e6;
else
    J = -log10(T1_avg);   % pura transmitancia -- sin peso de resolución
end
end

function J = costT2Sim(psiRx_deg, azimuth, target_center_2d, Rx_z, tg, n1, n2)
% Costo de la Etapa 2: -log10(T2 promedio), variando solo la distancia
% radial de Rx (a lo largo de un azimut fijo) via psi_Rx = atan(rho/Rx_z).

psiRx_rad = deg2rad(psiRx_deg);
rho = Rx_z * tan(psiRx_rad);
Rx_xy  = target_center_2d(:)' + rho * [cos(azimuth), sin(azimuth)];
Rx_pos = [Rx_xy(1); Rx_xy(2); Rx_z];

Ntg = size(tg, 2);
T2vals = zeros(1, Ntg);
for i = 1:Ntg
    target_pos = tg(:, i);
    intersec2  = calculateRefractionPointFermat(target_pos, Rx_pos, n2, n1);
    ang2       = calculateIncidenceAngle(target_pos, intersec2, Rx_pos);
    [~,~,~,~, T2vals(i), ~] = calculateTMcoef(ang2, n2, n1);
end
T2_avg = mean(T2vals);

if ~isfinite(T2_avg) || T2_avg <= 0
    J = 1e6;
else
    J = -log10(T2_avg);
end

end

function xStar = goldenSectionMin(costFun, xMin, xMax, tol)
gr = (sqrt(5) - 1) / 2;
a = xMin; b = xMax;
cc = b - gr*(b - a);
dd = a + gr*(b - a);
fc = costFun(cc);
fd = costFun(dd);
while abs(b - a) > tol
    if fc < fd
        b  = dd; dd = cc; fd = fc;
        cc = b - gr*(b - a);
        fc = costFun(cc);
    else
        a  = cc; cc = dd; fc = fd;
        dd = a + gr*(b - a);
        fd = costFun(dd);
    end
end
xStar = (a + b) / 2;
end
