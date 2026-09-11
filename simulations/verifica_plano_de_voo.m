clc;
clear;
close all;
%
% VERIFICA_PLANO_DE_VOO  Compara el desempeño del receptor bistático optimizado
% (run_plano_de_voo.m) contra la línea base monoestática (Rx = Tx, generada por
% run_plano_de_voo_monoestatico.m).
%
% Ambos scripts guardan un struct `Resultados` con los mismos campos
% (Resultados.optimizacion.potencia_maxima_dBm / delta_xy_m / delta_z_m /
% costo_combinado), usando la misma trayectoria helicoidal, el mismo grid de
% targets y el mismo alpha_res — así que se pueden comparar posición a posición
% sin reinterpolar ni reescalar nada.
%
% Por defecto toma el .mat MÁS RECIENTE de cada carpeta de salida. Para
% comparar corridas específicas, asignar explícitamente bistaticFile /
% monostaticFile más abajo.

%% Selección de archivos de entrada
bistaticDir  = fullfile('io', 'plan_vuelo');
monostaticDir = fullfile('io', 'plan_vuelo_monoestatico');

bistaticFile   = '';   % dejar vacío para usar el .mat más reciente de bistaticDir
monostaticFile = '';   % dejar vacío para usar el .mat más reciente de monostaticDir

if isempty(bistaticFile)
    bistaticFile = latestMatFile(bistaticDir, 'optimizacion_*.mat');
end
if isempty(monostaticFile)
    monostaticFile = latestMatFile(monostaticDir, 'monoestatico_*.mat');
end

fprintf('=== ARCHIVOS CARGADOS ===\n');
fprintf('Bistático (optimizado):  %s\n', bistaticFile);
fprintf('Monoestático (Rx = Tx):  %s\n', monostaticFile);

Sbi   = load(bistaticFile,   'Resultados'); Res_bi   = Sbi.Resultados;
Smono = load(monostaticFile, 'Resultados'); Res_mono = Smono.Resultados;

%% Verificar que ambas corridas son comparables (misma trayectoria Tx)
N_bi   = size(Res_bi.Tx_positions, 2);
N_mono = size(Res_mono.Tx_positions, 2);
if N_bi ~= N_mono
    error(['El número de posiciones Tx difiere entre corridas (bistático=%d, monoestático=%d).\n' ...
           'Volver a correr ambos scripts con los mismos JSON de parámetros antes de comparar.'], N_bi, N_mono);
end

tx_mismatch = max(vecnorm(Res_bi.Tx_positions - Res_mono.Tx_positions, 2, 1));
if tx_mismatch > 1e-3
    warning(['Las trayectorias Tx difieren hasta %.3f m entre las dos corridas. ' ...
             'Verificar que ambos scripts usaron los mismos parámetros JSON.'], tx_mismatch);
end

N = N_bi;
tx_idx = (1:N)';

%% Extraer métricas (mismos campos en ambos structs)
Pot_dBm_bi   = Res_bi.optimizacion.potencia_maxima_dBm(:);
Pot_dBm_mono = Res_mono.optimizacion.potencia_maxima_dBm(:);

delta_xy_bi   = Res_bi.optimizacion.delta_xy_m(:)   * 100;   % cm
delta_xy_mono = Res_mono.optimizacion.delta_xy_m(:) * 100;   % cm

delta_z_bi   = Res_bi.optimizacion.delta_z_m(:)   * 100;     % cm
delta_z_mono = Res_mono.optimizacion.delta_z_m(:) * 100;     % cm

costo_bi   = Res_bi.optimizacion.costo_combinado(:);
costo_mono = Res_mono.optimizacion.costo_combinado(:);

%% Ganancias del bistático sobre el monoestático (positivo = bistático mejor)
ganancia_potencia_dB = Pot_dBm_bi - Pot_dBm_mono;         % dBm - dBm = dB (ganancia directa)
mejora_delta_xy_cm   = delta_xy_mono - delta_xy_bi;         % >0: bistático tiene menor (mejor) δxy
mejora_delta_z_cm    = delta_z_mono  - delta_z_bi;          % >0: bistático tiene menor (mejor) δz
mejora_costo_J       = costo_mono - costo_bi;               % >0: bistático tiene menor (mejor) J

%% Resumen de diagnóstico
fprintf('\n=== RESUMEN: BISTÁTICO OPTIMIZADO vs. MONOESTÁTICO (Rx=Tx) ===\n');
fprintf('Posiciones Tx comparadas: %d\n\n', N);

print_metric_summary('Potencia [dB]',      ganancia_potencia_dB, N);
print_metric_summary('Área δxy·δz [%]',    100 * (1 - (delta_xy_bi.*delta_z_bi) ./ (delta_xy_mono.*delta_z_mono)), N);
print_metric_summary('Costo J',            mejora_costo_J, N);

fprintf('\n%-16s  %10s  %10s\n', 'Métrica', 'Bistático', 'Monoest.');
fprintf('%-16s  %10.2f  %10.2f\n', 'Pot media [dBm]', mean(Pot_dBm_bi),   mean(Pot_dBm_mono));
fprintf('%-16s  %10.2f  %10.2f\n', 'δxy media [cm]',  mean(delta_xy_bi),  mean(delta_xy_mono));
fprintf('%-16s  %10.2f  %10.2f\n', 'δz media [cm]',   mean(delta_z_bi),   mean(delta_z_mono));
fprintf('%-16s  %10.4f  %10.4f\n', 'J media',         mean(costo_bi),     mean(costo_mono));

%% Guardar CSV comparativo
outputDir = fullfile('io', 'plan_vuelo_comparacion');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
csv_path  = fullfile(outputDir, sprintf('comparacion_%s.csv', timestamp));

fid = fopen(csv_path, 'w');
fprintf(fid, ['Tx_ID,Pot_dBm_bi,Pot_dBm_mono,ganancia_potencia_dB,' ...
              'delta_xy_cm_bi,delta_xy_cm_mono,mejora_delta_xy_cm,' ...
              'delta_z_cm_bi,delta_z_cm_mono,mejora_delta_z_cm,' ...
              'costo_J_bi,costo_J_mono,mejora_costo_J\n']);
for i = 1:N
    fprintf(fid, '%d,%.3f,%.3f,%.3f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f\n', ...
        i, Pot_dBm_bi(i), Pot_dBm_mono(i), ganancia_potencia_dB(i), ...
        delta_xy_bi(i), delta_xy_mono(i), mejora_delta_xy_cm(i), ...
        delta_z_bi(i), delta_z_mono(i), mejora_delta_z_cm(i), ...
        costo_bi(i), costo_mono(i), mejora_costo_J(i));
end
fclose(fid);
fprintf('\nCSV comparativo guardado en: %s\n', csv_path);

%% Visualización
fig = figure('Position', [100, 100, 1400, 900]);
set(fig, 'WindowState', 'maximized');

subplot(3,2,1); hold on
plot(tx_idx, Pot_dBm_bi,   '-o', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', 'Bistatic optimized');
plot(tx_idx, Pot_dBm_mono, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'Monostatic (Rx=Tx)');
xlabel('Tx/Rx Position Number'); ylabel('Power [dBm]'); title('Received power');
legend('Location','best'); grid on; hold off;

subplot(3,2,2); hold on
bar(tx_idx, ganancia_potencia_dB, 'FaceColor', [0.47 0.67 0.19]);
plot(xlim, [0 0], 'k--', 'LineWidth', 1);
xlabel('Tx/Rx Position Number'); ylabel('Gain [dB]'); title('Power gain: bistatic − monostatic');
grid on; hold off;

subplot(3,2,3); hold on
plot(tx_idx, delta_xy_bi,   '-o', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', 'Bistatic');
plot(tx_idx, delta_xy_mono, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'Monostatic');
xlabel('Tx/Rx Position Number'); ylabel('\delta_{xy} [cm]'); title('Horizontal resolution');
legend('Location','best'); grid on; hold off;

subplot(3,2,4); hold on
plot(tx_idx, delta_z_bi,   '-o', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', 'Bistatic');
plot(tx_idx, delta_z_mono, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'Monostatic');
xlabel('Tx/Rx Position Number'); ylabel('\delta_z [cm]'); title('Vertical resolution');
legend('Location','best'); grid on; hold off;

subplot(3,2,5); hold on
plot(tx_idx, costo_bi,   '-o', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', 'Bistatic');
plot(tx_idx, costo_mono, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'Monostatic');
xlabel('Tx/Rx Position Number'); ylabel('Cost J'); title('Combined cost (lower = better)');
legend('Location','best'); grid on; hold off;

subplot(3,2,6);
histogram(ganancia_potencia_dB, min(15, max(5, round(N/3))), ...
    'FaceColor', [0.70 0.70 1.00], 'EdgeColor', 'black');
xlabel('Power gain [dB]'); ylabel('Frequency');
title('Distribution of the bistatic vs. monostatic gain'); grid on;

sgtitle('Bistatic optimized vs. Monostatic (Rx = Tx) over the same Tx helix', 'FontWeight','bold');

fig_path = fullfile(outputDir, sprintf('comparacion_%s.jpg', timestamp));
saveas(fig, fig_path);
fprintf('Figura comparativa guardada en: %s\n', fig_path);

%% Figura independiente: potencia recibida (mismo contenido que el subplot 3,2,1)
fig1b = figure('Position', [100, 100, 900, 600]);
hold on;
plot(tx_idx, Pot_dBm_bi,   '-o', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', 'Bistatic optimized');
plot(tx_idx, Pot_dBm_mono, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'Monostatic (Rx=Tx)');
xlabel('Tx/Rx Position Number'); ylabel('Power [dBm]'); title('Received power');
legend('Location','best'); grid on; hold off;

fig1b_path = fullfile(outputDir, sprintf('comparacion_potencia_%s.jpg', timestamp));
saveas(fig1b, fig1b_path);
fprintf('Figura de potencia (independiente) guardada en: %s\n', fig1b_path);

%% Figura independiente: ganancia de potencia (mismo contenido que el subplot 3,2,2)
fig2b = figure('Position', [100, 100, 900, 600]);
hold on;
bar(tx_idx, ganancia_potencia_dB, 'FaceColor', [0.47 0.67 0.19]);
plot(xlim, [0 0], 'k--', 'LineWidth', 1);
xlabel('Tx/Rx Position Number'); ylabel('Gain [dB]'); title('Power gain: bistatic − monostatic');
grid on; hold off;

fig2b_path = fullfile(outputDir, sprintf('comparacion_ganancia_%s.jpg', timestamp));
saveas(fig2b, fig2b_path);
fprintf('Figura de ganancia (independiente) guardada en: %s\n', fig2b_path);

%% Figura independiente: resolución horizontal (mismo contenido que el subplot 3,2,3)
fig3b = figure('Position', [100, 100, 900, 600]);
hold on;
plot(tx_idx, delta_xy_bi,   '-o', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', 'Bistatic');
plot(tx_idx, delta_xy_mono, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'Monostatic');
xlabel('Tx/Rx Position Number'); ylabel('\delta_{xy} [cm]'); title('Horizontal resolution');
legend('Location','best'); grid on; hold off;

fig3b_path = fullfile(outputDir, sprintf('comparacion_deltaxy_%s.jpg', timestamp));
saveas(fig3b, fig3b_path);
fprintf('Figura de resolución horizontal (independiente) guardada en: %s\n', fig3b_path);

%% Figura independiente: resolución vertical (mismo contenido que el subplot 3,2,4)
fig4b = figure('Position', [100, 100, 900, 600]);
hold on;
plot(tx_idx, delta_z_bi,   '-o', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', 'Bistatic');
plot(tx_idx, delta_z_mono, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'Monostatic');
xlabel('Tx/Rx Position Number'); ylabel('\delta_z [cm]'); title('Vertical resolution');
legend('Location','best'); grid on; hold off;

fig4b_path = fullfile(outputDir, sprintf('comparacion_deltaz_%s.jpg', timestamp));
saveas(fig4b, fig4b_path);
fprintf('Figura de resolución vertical (independiente) guardada en: %s\n', fig4b_path);

%% Figura independiente: costo combinado (mismo contenido que el subplot 3,2,5)
fig5b = figure('Position', [100, 100, 900, 600]);
hold on;
plot(tx_idx, costo_bi,   '-o', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', 'Bistatic');
plot(tx_idx, costo_mono, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'Monostatic');
xlabel('Tx/Rx Position Number'); ylabel('Cost J'); title('Combined cost (lower = better)');
legend('Location','best'); grid on; hold off;

fig5b_path = fullfile(outputDir, sprintf('comparacion_costo_%s.jpg', timestamp));
saveas(fig5b, fig5b_path);
fprintf('Figura de costo (independiente) guardada en: %s\n', fig5b_path);

%% Figura independiente: histograma de ganancia (mismo contenido que el subplot 3,2,6)
fig6b = figure('Position', [100, 100, 900, 600]);
histogram(ganancia_potencia_dB, min(15, max(5, round(N/3))), ...
    'FaceColor', [0.70 0.70 1.00], 'EdgeColor', 'black');
xlabel('Power gain [dB]'); ylabel('Frequency');
title('Distribution of the bistatic vs. monostatic gain'); grid on;

fig6b_path = fullfile(outputDir, sprintf('comparacion_histograma_%s.jpg', timestamp));
saveas(fig6b, fig6b_path);
fprintf('Figura de histograma (independiente) guardada en: %s\n', fig6b_path);

%% Figura "vendible": ganancia de potencia en número, para congreso
% Panel izquierdo: dumbbell (antes→después) por posición Tx — monoestático vs
% bistático, sin saturar de etiquetas. Panel derecho: los números que cierran
% el argumento (ganancia media en dB, ratio lineal, min/max), grandes y
% directos, para una diapositiva/póster.

col_mono = [0.00 0.45 0.74];   % mismo azul que el resto del script
col_bi   = [0.85 0.33 0.10];   % mismo naranja que el resto del script
col_stem = [0.65 0.65 0.65];   % gris recesivo, solo conecta — no es dato

fig2 = figure('Position', [100, 100, 1500, 650]);
set(fig2, 'WindowState', 'maximized');

% --- Panel izquierdo: dumbbell chart de potencia por posición Tx ---
ax1 = subplot('Position', [0.06 0.15 0.58 0.75]);
hold on;
for i = 1:N
    plot([tx_idx(i) tx_idx(i)], [Pot_dBm_mono(i) Pot_dBm_bi(i)], '-', ...
        'Color', col_stem, 'LineWidth', 1.2, 'HandleVisibility', 'off');
end
plot(tx_idx, Pot_dBm_mono, 'o', 'MarkerSize', 8, 'MarkerFaceColor', col_mono, ...
    'MarkerEdgeColor', 'white', 'LineWidth', 1, 'DisplayName', 'Monostatic (Rx=Tx)');
plot(tx_idx, Pot_dBm_bi, 'o', 'MarkerSize', 8, 'MarkerFaceColor', col_bi, ...
    'MarkerEdgeColor', 'white', 'LineWidth', 1, 'DisplayName', 'Bistatic optimized');

% Label only the point with the largest gain (avoid cluttering the chart with
% one number per point — the rest of the values are in the CSV/table)
[max_gain, idx_max_gain] = max(ganancia_potencia_dB);
text(tx_idx(idx_max_gain), Pot_dBm_bi(idx_max_gain) + 1.5, ...
    sprintf('+%.1f dB', max_gain), 'FontSize', 10, 'FontWeight', 'bold', ...
    'Color', col_bi, 'HorizontalAlignment', 'center');

xlabel('Tx/Rx Position Number');
ylabel('Received power [dBm]');
title('Received power per position: monostatic \rightarrow bistatic optimized');
legend('Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
grid on; box on; hold off;

% --- Right panel: headline numbers (hero figure + stat tiles) ---
ax2 = subplot('Position', [0.70 0.06 0.27 0.88]);
axis(ax2, 'off'); xlim(ax2, [0 1]); ylim(ax2, [0 1]); hold(ax2, 'on');

mean_gain = mean(ganancia_potencia_dB);
min_gain  = min(ganancia_potencia_dB);
ratio_lin = 10^(mean_gain/10);
casos_mejor = sum(ganancia_potencia_dB > 0);

% Hero figure: mean gain in dB (the number that carries the argument)
text(0.5, 0.86, sprintf('+%.1f dB', mean_gain), 'FontSize', 44, 'FontWeight', 'bold', ...
    'Color', col_bi, 'HorizontalAlignment', 'center', 'Parent', ax2);
text(0.5, 0.72, 'mean power gain', 'FontSize', 12, ...
    'Color', [0.30 0.30 0.30], 'HorizontalAlignment', 'center', 'Parent', ax2);
text(0.5, 0.67, sprintf('bistatic optimized vs. monostatic (N=%d positions)', N), ...
    'FontSize', 9, 'Color', [0.50 0.50 0.50], 'HorizontalAlignment', 'center', 'Parent', ax2);

% Secondary stat tile: same number in linear scale (intuitive reading)
text(0.5, 0.50, sprintf('\\times%.1f', ratio_lin), 'FontSize', 26, 'FontWeight', 'bold', ...
    'Color', [0.20 0.20 0.20], 'HorizontalAlignment', 'center', 'Parent', ax2);
text(0.5, 0.42, 'times more power (linear scale)', 'FontSize', 10, ...
    'Color', [0.40 0.40 0.40], 'HorizontalAlignment', 'center', 'Parent', ax2);

% Tertiary stat tiles: worst case and coverage
text(0.5, 0.27, sprintf('%.1f dB', min_gain), 'FontSize', 18, 'FontWeight', 'bold', ...
    'Color', [0.40 0.40 0.40], 'HorizontalAlignment', 'center', 'Parent', ax2);
text(0.5, 0.20, 'minimum gain (worst case)', 'FontSize', 9, ...
    'Color', [0.50 0.50 0.50], 'HorizontalAlignment', 'center', 'Parent', ax2);

text(0.5, 0.06, sprintf('%d of %d positions (%.0f%%) with positive gain', ...
    casos_mejor, N, 100*casos_mejor/N), 'FontSize', 10, 'Color', [0.40 0.40 0.40], ...
    'HorizontalAlignment', 'center', 'Parent', ax2);
hold(ax2, 'off');

sgtitle('Power gain: Bistatic optimized vs. Monostatic', ...
    'FontWeight', 'bold', 'FontSize', 14);

fig2_path = fullfile(outputDir, sprintf('comparacion_numerica_%s.jpg', timestamp));
saveas(fig2, fig2_path);
fprintf('Figura numérica (para congreso) guardada en: %s\n', fig2_path);

%% ==================== Funciones locales ====================

function filepath = latestMatFile(directory, pattern)
    files = dir(fullfile(directory, pattern));
    if isempty(files)
        error(['No se encontró ningún archivo "%s" en %s.\n' ...
               'Correr primero el script correspondiente (run_plano_de_voo.m o ' ...
               'run_plano_de_voo_monoestatico.m).'], pattern, directory);
    end
    [~, idx] = max([files.datenum]);
    filepath = fullfile(directory, files(idx).name);
end

function print_metric_summary(label, delta, N)
    casos_mejor = sum(delta > 0);
    fprintf('%-16s  media=%8.3f  min=%8.3f  max=%8.3f  | bistático mejor en %d/%d (%.1f%%)\n', ...
        label, mean(delta), min(delta), max(delta), casos_mejor, N, 100*casos_mejor/N);
end
