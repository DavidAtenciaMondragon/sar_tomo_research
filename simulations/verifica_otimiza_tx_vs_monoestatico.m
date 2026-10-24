clc;
clear;
close all;
%
% VERIFICA_OTIMIZA_TX_VS_MONOESTATICO  Compara el optimizador biestático de
% dos etapas (run_otimiza_tx.m: Etapa 1 -- Tx, Etapa 2 -- Rx) contra una
% línea base MONOESTÁTICA (Rx = Tx) construida sobre EXACTAMENTE LA MISMA
% geometría, para que la única diferencia entre ambos escenarios sea la
% estrategia del receptor.
%
% MOTIVACIÓN (por qué no simplemente usar io/plan_vuelo_monoestatico/*.mat):
% run_plano_de_voo_monoestatico.m genera su propia hélice de Tx a partir del
% psi0 NOMINAL del JSON (system_espiral_plano_voo.json), que es distinto del
% psi0* que la Etapa 1 de run_otimiza_tx.m encuentra óptimo. Comparar esos
% dos .mat directamente mezclaría dos efectos a la vez -- "mejor Tx" (Etapa 1)
% y "mejor Rx" (Etapa 2) -- y no aislaría el aporte de cada uno. Aquí, en
% cambio, el monoestático se recalcula tomando literalmente Resultados.Tx_positions
% del propio run_otimiza_tx.m (la hélice psi0* ya optimizada) y forzando
% Rx=Tx en cada una de esas posiciones, con los mismos n1/n2/Pt/sigma/lambda/
% B/B_helix/beta_helix/alpha_res que usó la Etapa 2. Así:
%
%   - Trayectoria del Tx:      IDÉNTICA en ambos escenarios (la de la Etapa 1)
%   - Grid de targets:         IDÉNTICO (Resultados.targets)
%   - Parámetros físicos:      IDÉNTICOS (Resultados.parametros / etapa1_final)
%   - Única diferencia:        posición del Rx -- optimizado (Etapa 2) vs Rx=Tx
%
% lo que hace la comparación punto a punto (misma N, misma trayectoria Tx)
% tan rigurosa como la de verifica_plano_de_voo.m, pero aplicada al método de
% dos etapas en vez del optimizador de solo-Rx.
%
% Requiere: haber corrido antes run_otimiza_tx.m (io/otimiza_tx/otimiza_tx_*.mat).
% Por defecto toma el .mat más reciente de esa carpeta.

addpath(genpath('gs'))
addpath(genpath('proc'))
addpath(genpath('tools'))
addpath(genpath('flightpath'))

%% ==================== Cargar la corrida del optimizador de dos etapas ====================
biDir  = fullfile('io', 'otimiza_tx');
biFile = '';   % dejar vacío para usar el .mat más reciente de biDir

if isempty(biFile)
    biFile = latestMatFile(biDir, 'otimiza_tx_*.mat');
end

fprintf('=== ARCHIVO CARGADO ===\n');
fprintf('Optimizador de dos etapas (Etapa 1 + Etapa 2):  %s\n', biFile);

Sbi    = load(biFile, 'Resultados');
Res_bi = Sbi.Resultados;

%% Geometría y parámetros físicos (reusados TAL CUAL para el monoestático)
Tx_pos_all       = Res_bi.Tx_positions;          % [3 x N] -- hélice psi0* de la Etapa 1
tg               = Res_bi.targets;               % [3 x Ntg] -- mismo grid
target_center_2d = Res_bi.target_center;         % [2x1]
N                = size(Tx_pos_all, 2);

n1        = Res_bi.parametros.n1;
n2        = Res_bi.parametros.n2;
lambda    = Res_bi.parametros.lambda;
B         = Res_bi.parametros.B;
B_helix   = Res_bi.parametros.B_helix;
beta_helix = Res_bi.etapa1_final.beta_helix;     % hélice REAL de la Etapa 1 (beta=psi0*)
Pt        = Res_bi.parametros.Pt;
sigma     = Res_bi.parametros.sigma;
alpha_res = Res_bi.parametros.alpha_res;         % mismo peso potencia/resolución -> J comparable
theta_B   = Res_bi.parametros.theta_B;
psi0_star = Res_bi.etapa1_final.psi0_star;

% El patrón de radiación no se guarda en el .mat (es una matriz derivada de
% la apertura de antena) -- se reconstruye desde el mismo JSON que usó
% run_otimiza_tx.m. Asume que radarTx_espiral_plano_voo.json no cambió desde
% esa corrida (igual que hace run_plano_de_voo_monoestatico.m).
radarJSON  = json2struct(strcat('parametros', filesep, 'radarTx_espiral_plano_voo.json'));
strRadarTx = radarJSON.radar; clear radarJSON;
radPattern = createRadiationPattern(strRadarTx.AperturaElev, strRadarTx.AperturaAzimut);

fprintf('\n=== ESCENARIO COMPARTIDO (Etapa 1 -> hélice psi0*) ===\n');
fprintf('psi0* = %.2f°   theta_B (Brewster) = %.2f°   escenario Etapa 1 = %s\n', ...
    psi0_star, theta_B*180/pi, Res_bi.etapa1_final.scenario);
fprintf('n1=%.1f, n2=%.1f, lambda=%.3f m, B=%.0f MHz, B_helix=%.2f m, alpha_res=%.2f\n', ...
    n1, n2, lambda, B/1e6, B_helix, alpha_res);
fprintf('N posiciones Tx/Rx comparadas: %d  |  Ntargets = %d\n', N, size(tg,2));

%% ==================== Recalcular la línea base monoestática (Rx = Tx) ====================
% Misma metodología que run_plano_de_voo_monoestatico.m (objective_function +
% calculateBistaticResolution, sin fmincon porque Rx=Tx es la única posición
% posible), pero evaluada sobre la hélice psi0* de la Etapa 1 en vez de la
% hélice nominal del JSON.
max_power_mono = zeros(1, N);
delta_xy_mono  = zeros(1, N);
delta_z_mono   = zeros(1, N);
costo_mono      = zeros(1, N);
psi0_mono      = zeros(1, N);   % ángulo de mira Tx=Rx hacia el centroide (diagnóstico)

% Diagnóstico del lado biestático: reconstruir psi0/DeltaPhi del Rx óptimo de
% la Etapa 2 (no se guardaron en el .mat, solo delta_xy/delta_z resultantes)
Rx_opt_bi   = Res_bi.optimizacion.Rx_optimo;
psi0_bi     = zeros(1, N);
deltaPhi_bi = zeros(1, N);

fprintf('\nRecalculando monoestático (Rx=Tx) sobre la MISMA hélice Tx de la Etapa 1...\n');
for k = 1:N
    Tx_pos_current = Tx_pos_all(:, k);
    Rx_z           = Tx_pos_current(3);
    xy_mono        = Tx_pos_current(1:2);

    [~, ~, Gt_current]  = calculateTxGainsForTargets(Tx_pos_current, tg, radPattern);
    [R_T_pre, T1_pre]   = precomputeTxData(Tx_pos_current, tg, n1, n2);

    costo_mono(k) = objective_function(xy_mono, tg, Tx_pos_current, Rx_z, radPattern, ...
        Gt_current, Pt, sigma, lambda, n1, n2, ...
        alpha_res, B, B_helix, beta_helix, target_center_2d, R_T_pre, T1_pre);

    c_power = objective_function(xy_mono, tg, Tx_pos_current, Rx_z, radPattern, ...
        Gt_current, Pt, sigma, lambda, n1, n2, ...
        0, B, B_helix, beta_helix, target_center_2d, R_T_pre, T1_pre);
    max_power_mono(k) = 10^(-c_power);

    [delta_xy_mono(k), delta_z_mono(k), psi0_mono(k)] = calculateBistaticResolution(...
        Tx_pos_current, Tx_pos_current, target_center_2d, n2, lambda, B, B_helix, beta_helix);

    [~, ~, psi0_bi(k), deltaPhi_bi(k)] = calculateBistaticResolution(...
        Tx_pos_current, Rx_opt_bi(:,k), target_center_2d, n2, lambda, B, B_helix, beta_helix);
end

%% Métricas del optimizador de dos etapas (ya calculadas por run_otimiza_tx.m)
Pot_dBm_bi   = Res_bi.optimizacion.potencia_maxima_dBm(:)';
delta_xy_bi  = Res_bi.optimizacion.delta_xy_m(:)' * 100;   % cm
delta_z_bi   = Res_bi.optimizacion.delta_z_m(:)'  * 100;   % cm
costo_bi     = Res_bi.optimizacion.costo_combinado(:)';

Pot_dBm_mono  = 10*log10(max_power_mono * 1000);
delta_xy_mono = delta_xy_mono * 100;   % cm
delta_z_mono  = delta_z_mono  * 100;   % cm

%% Ganancias del optimizador de dos etapas sobre el monoestático (>0 = mejor)
ganancia_potencia_dB = Pot_dBm_bi - Pot_dBm_mono;
mejora_delta_xy_cm   = delta_xy_mono - delta_xy_bi;
mejora_delta_z_cm    = delta_z_mono  - delta_z_bi;
mejora_costo_J       = costo_mono     - costo_bi;

%% ==================== Resumen de diagnóstico ====================
fprintf('\n=== RESUMEN: DOS ETAPAS (Tx+Rx) vs. MONOESTÁTICO (Rx=Tx), MISMA HÉLICE Tx ===\n');
fprintf('Posiciones comparadas: %d\n\n', N);

print_metric_summary('Potencia [dB]',   ganancia_potencia_dB, N);
print_metric_summary('Área δxy·δz [%]', 100 * (1 - (delta_xy_bi.*delta_z_bi) ./ (delta_xy_mono.*delta_z_mono)), N);
print_metric_summary('Costo J',         mejora_costo_J, N);

fprintf('\n%-16s  %10s  %10s\n', 'Métrica', 'Dos etapas', 'Monoest.');
fprintf('%-16s  %10.2f  %10.2f\n', 'Pot media [dBm]', mean(Pot_dBm_bi),  mean(Pot_dBm_mono));
fprintf('%-16s  %10.2f  %10.2f\n', 'δxy media [cm]',  mean(delta_xy_bi), mean(delta_xy_mono));
fprintf('%-16s  %10.2f  %10.2f\n', 'δz media [cm]',   mean(delta_z_bi),  mean(delta_z_mono));
fprintf('%-16s  %10.4f  %10.4f\n', 'J media',         mean(costo_bi),    mean(costo_mono));

fprintf('\n=== ÁNGULO DE MIRA DEL Rx vs. BREWSTER (%.2f°) ===\n', theta_B*180/pi);
fprintf('%-24s  %10s  %10s\n', '', 'Dos etapas', 'Monoest.');
fprintf('%-24s  %10.2f  %10.2f  [°, medio]\n', 'psi0 del Rx', mean(psi0_bi)*180/pi, mean(psi0_mono)*180/pi);
dist_bi_brewster   = mean(abs(psi0_bi - theta_B))*180/pi;
dist_mono_brewster = mean(abs(psi0_mono - theta_B))*180/pi;
fprintf('%-24s  %10.2f  %10.2f  [°, dist. media a Brewster]\n', '|psi0 - theta_B|', ...
    dist_bi_brewster, dist_mono_brewster);

if alpha_res > 0 && dist_bi_brewster < dist_mono_brewster
    fprintf(['\nEn el monoestático, el mismo ángulo psi0 debe servir tanto para T1 (Tx->target) como\n' ...
             'para T2 (target->Rx): no puede acercarse a Brewster por ambos lados a la vez. La Etapa 2\n' ...
             'mueve el Rx de forma independiente del Tx y consigue un psi0 más cercano a Brewster para T2,\n' ...
             'lo que en parte explica la ganancia de potencia medida arriba.\n']);
elseif alpha_res == 0
    fprintf(['\nalpha_res=0 en esta corrida: la Etapa 2 optimizó potencia PURA, sin penalizar la\n' ...
             'resolución. El Rx puede terminar lejos del ángulo de Brewster (como aquí) porque le\n' ...
             'conviene más reducir la distancia geométrica R_R que acercar T2 a su máximo -- la ganancia\n' ...
             'de potencia frente al monoestático viene entonces del 1/R^2 y la ganancia de antena, no\n' ...
             'solo de Brewster. Para una comparación centrada en el efecto Brewster, re-correr\n' ...
             'run_otimiza_tx.m con alpha_res>0 (p.ej. 0.3, el valor por defecto de otimiza_tx.json).\n']);
else
    fprintf(['\nAquí el Rx del monoestático (forzado a psi0 del Tx) queda casualmente más cerca de\n' ...
             'Brewster que el Rx optimizado de la Etapa 2; la ganancia de potencia medida arriba viene\n' ...
             'entonces de otros términos de la ecuación de radar (1/R^2, ganancia de antena), no de T2.\n']);
end

%% ==================== Guardar resultados ====================
outputDir = fullfile('io', 'otimiza_tx_vs_monoestatico');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% .mat del monoestático recalculado (mismos nombres de campo que los demás
% scripts run_*.m, para poder reusar herramientas existentes si se desea)
Resultados.Tx_positions   = Tx_pos_all;
Resultados.targets        = tg;
Resultados.target_center  = target_center_2d;
Resultados.parametros     = Res_bi.parametros;
Resultados.parametros.beta_helix = beta_helix;
Resultados.optimizacion.Rx_optimo           = Tx_pos_all;   % Rx = Tx
Resultados.optimizacion.potencia_maxima_W   = max_power_mono;
Resultados.optimizacion.potencia_maxima_dBm = Pot_dBm_mono;
Resultados.optimizacion.delta_xy_m          = delta_xy_mono / 100;
Resultados.optimizacion.delta_z_m           = delta_z_mono  / 100;
Resultados.optimizacion.costo_combinado     = costo_mono;
Resultados.fecha_calculo   = datestr(now);
Resultados.fuente_geometria = biFile;   % de dónde vino la hélice Tx compartida

mono_mat_path = fullfile(outputDir, sprintf('monoestatico_matched_%s.mat', timestamp));
save(mono_mat_path, 'Resultados');
fprintf('\n.mat del monoestático (hélice compartida) guardado en: %s\n', mono_mat_path);

% CSV comparativo punto a punto
csv_path = fullfile(outputDir, sprintf('comparacion_%s.csv', timestamp));
fid = fopen(csv_path, 'w');
fprintf(fid, ['Tx_ID,Pot_dBm_dosEtapas,Pot_dBm_mono,ganancia_potencia_dB,' ...
              'delta_xy_cm_dosEtapas,delta_xy_cm_mono,mejora_delta_xy_cm,' ...
              'delta_z_cm_dosEtapas,delta_z_cm_mono,mejora_delta_z_cm,' ...
              'costo_J_dosEtapas,costo_J_mono,mejora_costo_J,' ...
              'psi0_rx_dosEtapas_deg,psi0_rx_mono_deg\n']);
for i = 1:N
    fprintf(fid, '%d,%.3f,%.3f,%.3f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.3f,%.3f\n', ...
        i, Pot_dBm_bi(i), Pot_dBm_mono(i), ganancia_potencia_dB(i), ...
        delta_xy_bi(i), delta_xy_mono(i), mejora_delta_xy_cm(i), ...
        delta_z_bi(i), delta_z_mono(i), mejora_delta_z_cm(i), ...
        costo_bi(i), costo_mono(i), mejora_costo_J(i), ...
        psi0_bi(i)*180/pi, psi0_mono(i)*180/pi);
end
fclose(fid);
fprintf('CSV comparativo guardado en: %s\n', csv_path);

%% ==================== Figura 1: comparación por posición Tx ====================
tx_idx = 1:N;

fig = figure('Position', [100, 100, 1400, 900]);
set(fig, 'WindowState', 'maximized');

subplot(3,2,1); hold on
plot(tx_idx, Pot_dBm_bi,   '-o', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', 'Dos etapas (Tx+Rx)');
plot(tx_idx, Pot_dBm_mono, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'Monoestático (Rx=Tx)');
xlabel('Posición Tx/Rx'); ylabel('Potencia [dBm]'); title('Potencia recibida');
legend('Location','best'); grid on; hold off;

subplot(3,2,2); hold on
bar(tx_idx, ganancia_potencia_dB, 'FaceColor', [0.47 0.67 0.19]);
plot(xlim, [0 0], 'k--', 'LineWidth', 1);
xlabel('Posición Tx/Rx'); ylabel('Ganancia [dB]'); title('Ganancia de potencia: dos etapas − monoestático');
grid on; hold off;

subplot(3,2,3); hold on
plot(tx_idx, delta_xy_bi,   '-o', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', 'Dos etapas');
plot(tx_idx, delta_xy_mono, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'Monoestático');
xlabel('Posición Tx/Rx'); ylabel('\delta_{xy} [cm]'); title('Resolución horizontal');
legend('Location','best'); grid on; hold off;

subplot(3,2,4); hold on
plot(tx_idx, delta_z_bi,   '-o', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', 'Dos etapas');
plot(tx_idx, delta_z_mono, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'Monoestático');
xlabel('Posición Tx/Rx'); ylabel('\delta_z [cm]'); title('Resolución vertical');
legend('Location','best'); grid on; hold off;

subplot(3,2,5); hold on
plot(tx_idx, costo_bi,   '-o', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', 'Dos etapas');
plot(tx_idx, costo_mono,  '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'Monoestático');
xlabel('Posición Tx/Rx'); ylabel('Costo J'); title('Costo combinado (menor = mejor)');
legend('Location','best'); grid on; hold off;

subplot(3,2,6);
histogram(ganancia_potencia_dB, min(15, max(5, round(N/3))), ...
    'FaceColor', [0.70 0.70 1.00], 'EdgeColor', 'black');
xlabel('Ganancia de potencia [dB]'); ylabel('Frecuencia');
title('Distribución de la ganancia (dos etapas vs. monoestático)'); grid on;

sgtitle('Optimizador de dos etapas (Tx+Rx) vs. Monoestático (Rx=Tx) -- misma hélice Tx', 'FontWeight', 'bold');
fig_path = fullfile(outputDir, sprintf('comparacion_%s.jpg', timestamp));
saveas(fig, fig_path);
fprintf('Figura comparativa guardada en: %s\n', fig_path);

%% ==================== Figura 2: trayectorias 3D + planta ====================
fig2 = figure('Position', [100, 100, 1500, 700]);
set(fig2, 'WindowState', 'maximized');

col_tx   = [0.00 0.45 0.74];
col_bi   = [0.85 0.33 0.10];

ax1 = subplot(1,2,1); hold(ax1, 'on');
plot3(ax1, Tx_pos_all(1,:), Tx_pos_all(2,:), Tx_pos_all(3,:), '-o', 'Color', col_tx, ...
    'LineWidth', 1.3, 'MarkerSize', 3, 'DisplayName', 'Tx (= Rx monoestático)');
plot3(ax1, Rx_opt_bi(1,:), Rx_opt_bi(2,:), Rx_opt_bi(3,:), '--s', 'Color', col_bi, ...
    'LineWidth', 1.1, 'MarkerSize', 3, 'MarkerFaceColor', col_bi, 'DisplayName', 'Rx óptimo (dos etapas)');
for i = 1:N
    plot3(ax1, [Tx_pos_all(1,i), Rx_opt_bi(1,i)], [Tx_pos_all(2,i), Rx_opt_bi(2,i)], ...
        [Tx_pos_all(3,i), Rx_opt_bi(3,i)], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.4, ...
        'HandleVisibility', 'off');
end
plot3(ax1, tg(1,:), tg(2,:), tg(3,:), 'd', 'Color', 'black', 'MarkerSize', 4, ...
    'LineWidth', 1, 'MarkerFaceColor', 'red', 'DisplayName', 'Blancos subsuperficiales');

ext = max([abs(Tx_pos_all(1:2,:)), abs(Rx_opt_bi(1:2,:))], [], 'all');
[Xs, Ys] = meshgrid(-ext:10:ext, -ext:10:ext);
surf(ax1, Xs, Ys, zeros(size(Xs)), 'FaceAlpha', 0.20, 'FaceColor', [0.8 0.8 0.6], ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');

xlabel(ax1, 'X [m]'); ylabel(ax1, 'Y [m]'); zlabel(ax1, 'Z [m]');
title(ax1, 'Vista 3D'); legend(ax1, 'Location', 'northeast', 'FontSize', 8);
grid(ax1, 'on'); view(ax1, 45, 20); axis(ax1, 'tight'); box(ax1, 'on'); hold(ax1, 'off');

ax2 = subplot(1,2,2); hold(ax2, 'on');
plot(ax2, Tx_pos_all(1,:), Tx_pos_all(2,:), '-o', 'Color', col_tx, 'LineWidth', 1.3, ...
    'MarkerSize', 3, 'DisplayName', 'Tx (= Rx monoestático)');
plot(ax2, Rx_opt_bi(1,:), Rx_opt_bi(2,:), '--s', 'Color', col_bi, 'LineWidth', 1.1, ...
    'MarkerSize', 3, 'MarkerFaceColor', col_bi, 'DisplayName', 'Rx óptimo (dos etapas)');
plot(ax2, tg(1,:), tg(2,:), 'd', 'Color', 'black', 'MarkerSize', 4, ...
    'LineWidth', 1, 'MarkerFaceColor', 'red', 'DisplayName', 'Blancos');

xlabel(ax2, 'X [m]'); ylabel(ax2, 'Y [m]');
title(ax2, 'Planta (vista superior)'); legend(ax2, 'Location', 'best', 'FontSize', 8);
grid(ax2, 'on'); axis(ax2, 'equal'); box(ax2, 'on'); hold(ax2, 'off');

sgtitle('Trayectorias: Tx compartido, Rx dos etapas vs. Rx=Tx (monoestático)', 'FontWeight', 'bold');
fig2_path = fullfile(outputDir, sprintf('comparacion_3d_%s.jpg', timestamp));
saveas(fig2, fig2_path);
fprintf('Figura (3D + planta) guardada en: %s\n', fig2_path);

%% ==================== Figura 3: números para congreso ====================
col_mono = [0.00 0.45 0.74];
col_bi2  = [0.85 0.33 0.10];
col_stem = [0.65 0.65 0.65];

fig3 = figure('Position', [100, 100, 1500, 650]);
set(fig3, 'WindowState', 'maximized');

ax1b = subplot('Position', [0.06 0.15 0.58 0.75]);
hold on;
for i = 1:N
    plot([tx_idx(i) tx_idx(i)], [Pot_dBm_mono(i) Pot_dBm_bi(i)], '-', ...
        'Color', col_stem, 'LineWidth', 1.2, 'HandleVisibility', 'off');
end
plot(tx_idx, Pot_dBm_mono, 'o', 'MarkerSize', 8, 'MarkerFaceColor', col_mono, ...
    'MarkerEdgeColor', 'white', 'LineWidth', 1, 'DisplayName', 'Monoestático (Rx=Tx)');
plot(tx_idx, Pot_dBm_bi, 'o', 'MarkerSize', 8, 'MarkerFaceColor', col_bi2, ...
    'MarkerEdgeColor', 'white', 'LineWidth', 1, 'DisplayName', 'Dos etapas (Tx+Rx)');

[max_gain, idx_max_gain] = max(ganancia_potencia_dB);
text(tx_idx(idx_max_gain), Pot_dBm_bi(idx_max_gain) + 1.5, ...
    sprintf('+%.1f dB', max_gain), 'FontSize', 10, 'FontWeight', 'bold', ...
    'Color', col_bi2, 'HorizontalAlignment', 'center');

xlabel('Posición Tx/Rx');
ylabel('Potencia recibida [dBm]');
title('Potencia por posición: monoestático \rightarrow dos etapas');
legend('Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
grid on; box on; hold off;

ax2b = subplot('Position', [0.70 0.06 0.27 0.88]);
axis(ax2b, 'off'); xlim(ax2b, [0 1]); ylim(ax2b, [0 1]); hold(ax2b, 'on');

mean_gain   = mean(ganancia_potencia_dB);
min_gain    = min(ganancia_potencia_dB);
ratio_lin   = 10^(mean_gain/10);
casos_mejor = sum(ganancia_potencia_dB > 0);

text(0.5, 0.86, sprintf('+%.1f dB', mean_gain), 'FontSize', 44, 'FontWeight', 'bold', ...
    'Color', col_bi2, 'HorizontalAlignment', 'center', 'Parent', ax2b);
text(0.5, 0.72, 'ganancia media de potencia', 'FontSize', 12, ...
    'Color', [0.30 0.30 0.30], 'HorizontalAlignment', 'center', 'Parent', ax2b);
text(0.5, 0.67, sprintf('dos etapas vs. monoestático (N=%d posiciones, misma hélice Tx)', N), ...
    'FontSize', 9, 'Color', [0.50 0.50 0.50], 'HorizontalAlignment', 'center', 'Parent', ax2b);

text(0.5, 0.50, sprintf('\\times%.1f', ratio_lin), 'FontSize', 26, 'FontWeight', 'bold', ...
    'Color', [0.20 0.20 0.20], 'HorizontalAlignment', 'center', 'Parent', ax2b);
text(0.5, 0.42, 'veces más potencia (escala lineal)', 'FontSize', 10, ...
    'Color', [0.40 0.40 0.40], 'HorizontalAlignment', 'center', 'Parent', ax2b);

text(0.5, 0.27, sprintf('%.1f dB', min_gain), 'FontSize', 18, 'FontWeight', 'bold', ...
    'Color', [0.40 0.40 0.40], 'HorizontalAlignment', 'center', 'Parent', ax2b);
text(0.5, 0.20, 'ganancia mínima (peor caso)', 'FontSize', 9, ...
    'Color', [0.50 0.50 0.50], 'HorizontalAlignment', 'center', 'Parent', ax2b);

text(0.5, 0.06, sprintf('%d de %d posiciones (%.0f%%) con ganancia positiva', ...
    casos_mejor, N, 100*casos_mejor/N), 'FontSize', 10, 'Color', [0.40 0.40 0.40], ...
    'HorizontalAlignment', 'center', 'Parent', ax2b);
hold(ax2b, 'off');

sgtitle('Ganancia de potencia: Dos etapas (Tx+Rx) vs. Monoestático', ...
    'FontWeight', 'bold', 'FontSize', 14);

fig3_path = fullfile(outputDir, sprintf('comparacion_numerica_%s.jpg', timestamp));
saveas(fig3, fig3_path);
fprintf('Figura numérica (para congreso) guardada en: %s\n', fig3_path);

%% ==================== Funciones locales ====================

function filepath = latestMatFile(directory, pattern)
    files = dir(fullfile(directory, pattern));
    if isempty(files)
        error(['No se encontró ningún archivo "%s" en %s.\n' ...
               'Correr primero run_otimiza_tx.m.'], pattern, directory);
    end
    [~, idx] = max([files.datenum]);
    filepath = fullfile(directory, files(idx).name);
end

function print_metric_summary(label, delta, N)
    casos_mejor = sum(delta > 0);
    fprintf('%-16s  media=%8.3f  min=%8.3f  max=%8.3f  | dos etapas mejor en %d/%d (%.1f%%)\n', ...
        label, mean(delta), min(delta), max(delta), casos_mejor, N, 100*casos_mejor/N);
end
