clc;
clear;
close all;
%
% VERIFICA_OTIMIZA_TX_VS_PLANVOO  Compara N enfoques de planificación de
% vuelo biestático (N>=2), cada uno el resultado de uno de los scripts
% run_*.m / otimiza_tx*.m del proyecto. Por defecto compara los tres
% enfoques desarrollados hasta ahora:
%
%   1) Solo Rx (línea base)          -- run_plano_de_voo.m: Tx fijo en el
%      psi0 del JSON, solo Rx optimizado por la ecuación de radar completa.
%   2) Tx+Rx, ecuación de radar      -- run_otimiza_tx.m: Etapa 1 (Tx) +
%      Etapa 2 (Rx), ambas optimizan la ecuación de radar completa
%      (potencia y/o resolución, según alpha_res/gamma).
%   3) Tx+Rx, transmitancia pura     -- otimiza_tx_brewster_puro.m: Etapa 1
%      + Etapa 2, ambas optimizan SOLO T1/T2 (ángulo de Brewster puro), sin
%      ganancia de antena ni 1/R^2.
%
% Las trayectorias de Tx de cada enfoque tienen geometría distinta (psi0
% distinto), así que el número de posiciones decimadas puede diferir entre
% enfoques -- la comparación se hace sobre estadísticas agregadas (media,
% desviación estándar) y sobre la fracción normalizada de la trayectoria
% (0->1), no posición a posición (a diferencia de verifica_plano_de_voo.m,
% que sí compara punto a punto porque ahí la trayectoria de Tx es idéntica).
%
% Por defecto toma el .mat MÁS RECIENTE de cada carpeta. Para comparar
% corridas específicas, fijar el campo 'file' de la entrada correspondiente
% en RUNS más abajo. Para comparar solo 2 enfoques (o agregar un cuarto),
% simplemente editar la lista RUNS -- el resto del script es genérico.

addpath(genpath('tools'))
addpath(genpath('flightpath'))

%% ==================== CONFIGURACIÓN: enfoques a comparar ====================
% Cada entrada: label (para leyendas/tablas), dir+pattern (para tomar el
% .mat más reciente) o file (ruta explícita, tiene prioridad sobre dir+pattern).
RUNS = {
    struct('label', 'Solo Rx (línea base)',      'dir', fullfile('io','plan_vuelo'),               'pattern', 'optimizacion_*.mat',              'file', '')
    struct('label', 'Tx+Rx (ecuación radar)',    'dir', fullfile('io','otimiza_tx'),                'pattern', 'otimiza_tx_*.mat',                'file', '')
    struct('label', 'Tx+Rx (transmitancia pura)','dir', fullfile('io','otimiza_tx_brewster_puro'),  'pattern', 'otimiza_tx_brewster_puro_*.mat',  'file', '')
};
%% ==============================================================================

M = numel(RUNS);
if M < 2
    error('verifica_otimiza_tx_vs_planvoo:pocosRuns', 'RUNS debe tener al menos 2 enfoques.');
end

% Paleta de colores (se recicla si hay más de 6 enfoques)
palette = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; ...
           0.49 0.18 0.56; 0.93 0.69 0.13; 0.30 0.75 0.93];

%% Cargar cada enfoque
fprintf('=== ARCHIVOS CARGADOS ===\n');
data = struct([]);
for k = 1:M
    r = RUNS{k};
    if isempty(r.file)
        matFile = latestMatFile(r.dir, r.pattern);
    else
        matFile = r.file;
    end
    fprintf('%-28s  %s\n', r.label, matFile);
    S = load(matFile, 'Resultados');
    data = [data, extractRunData(S.Resultados, r.label, matFile)]; %#ok<AGROW>
end

n1 = data(1).n1; n2 = data(1).n2;
theta_B_air = atan(n2/n1);

%% ==================== Resumen numérico ====================
fprintf('\n=== GEOMETRÍA DEL Tx (psi0) Y MÉTRICAS AGREGADAS ===\n');
fprintf('%-28s  %8s  %5s  %10s  %10s  %10s  %14s\n', ...
    'Enfoque', 'psi0 [°]', 'N', 'Pot [dBm]', 'δxy [cm]', 'δz [cm]', 'Costo (ver nota)');
for k = 1:M
    fprintf('%-28s  %8.2f  %5d  %10.2f  %10.2f  %10.2f  %14.4f\n', ...
        data(k).label, data(k).psi0_deg, data(k).N, data(k).mean_pot, ...
        data(k).mean_dxy, data(k).mean_dz, data(k).mean_cost);
end
fprintf(['\nNota: el "costo" no está en la misma base entre enfoques -- para "ecuación radar"\n' ...
         'es J=-(1-a)log10(Pr)+a*log10(dxy*dz); para "transmitancia pura" es -log10(T2). No comparar\n' ...
         'el valor numérico del costo entre esos dos, solo dentro del mismo enfoque a lo largo del vuelo.\n']);

fprintf('\n=== GANANCIA DE CADA ENFOQUE FRENTE A LA LÍNEA BASE (enfoque #1: %s) ===\n', data(1).label);
fprintf('%-28s  %10s  %12s  %12s\n', 'Enfoque', 'ΔPot [dB]', 'Δδxy [%]', 'Δδz [%]');
for k = 2:M
    dPot = data(k).mean_pot - data(1).mean_pot;
    dDxy = 100*(data(k).mean_dxy - data(1).mean_dxy)/data(1).mean_dxy;
    dDz  = 100*(data(k).mean_dz  - data(1).mean_dz)/data(1).mean_dz;
    fprintf('%-28s  %+10.2f  %+11.1f%%  %+11.1f%%\n', data(k).label, dPot, dDxy, dDz);

    % Ganancia teórica predicha solo por T1(psi0) [Proposición 5.1, aditividad]
    T1_base_theory = fresnelT1(deg2rad(data(1).psi0_deg), n1, n2);
    T1_k_theory    = fresnelT1(deg2rad(data(k).psi0_deg), n1, n2);
    ganancia_T1_dB = 10*log10(T1_k_theory / T1_base_theory);
    fprintf('  -> ganancia predicha solo por T1(psi0) [Prop. 5.1]: %+.3f dB  (medida: %+.3f dB)\n', ...
        ganancia_T1_dB, dPot);
end

%% Guardar CSV resumen
outputDir = fullfile('io', 'otimiza_tx_vs_planvoo');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

csv_path = fullfile(outputDir, sprintf('comparacion_%s.csv', timestamp));
fid = fopen(csv_path, 'w');
fprintf(fid, 'enfoque,archivo,psi0_deg,N,potencia_media_dBm,delta_xy_media_cm,delta_z_media_cm,costo_medio\n');
for k = 1:M
    fprintf(fid, '%s,%s,%.4f,%d,%.4f,%.4f,%.4f,%.6f\n', ...
        data(k).label, data(k).file, data(k).psi0_deg, data(k).N, ...
        data(k).mean_pot, data(k).mean_dxy, data(k).mean_dz, data(k).mean_cost);
end
fclose(fid);
fprintf('\nCSV resumen guardado en: %s\n', csv_path);

%% ==================== Figura 1: series vs. fracción de trayectoria ====================
fig = figure('Position', [100, 100, 1400, 850]);
set(fig, 'WindowState', 'maximized');

metricNames  = {'pot', 'dz', 'dxy', 'cost'};
metricLabels = {'Potencia [dBm]', '\delta_z [cm]', '\delta_{xy} [cm]', 'Costo (no comparable entre enfoques)'};
metricTitles = {'Potencia recibida', 'Resolución vertical', 'Resolución horizontal', 'Costo optimizado (propio de cada enfoque)'};

for m = 1:4
    subplot(2,2,m); hold on;
    for k = 1:M
        frac = (0:data(k).N-1)' / max(data(k).N-1, 1);
        col  = palette(mod(k-1, size(palette,1))+1, :);
        plot(frac, data(k).(metricNames{m}), '-', 'Color', col, 'LineWidth', 1.1, ...
            'Marker', 'o', 'MarkerSize', 2.5, 'DisplayName', data(k).label);
    end
    xlabel('Fracción de la trayectoria (0\rightarrow1)'); ylabel(metricLabels{m});
    title(metricTitles{m}); legend('Location', 'best', 'FontSize', 8); grid on; hold off;
end

sgtitle('Comparación de enfoques a lo largo de la trayectoria', 'FontWeight', 'bold');
fig_path = fullfile(outputDir, sprintf('comparacion_%s.jpg', timestamp));
saveas(fig, fig_path);
fprintf('Figura guardada en: %s\n', fig_path);

%% ==================== Figura 2: barras (media +/- desv. estándar) ====================
fig2 = figure('Position', [100, 100, 1300, 500]);

subplot(1,3,1);
bar(1:M, [data.mean_pot]); hold on;
errorbar(1:M, [data.mean_pot], [data.std_pot], 'k.', 'LineWidth', 1.2);
set(gca, 'XTickLabel', {data.label}, 'XTickLabelRotation', 20, 'FontSize', 7.5);
ylabel('Potencia [dBm]'); title('Potencia media'); grid on; hold off;

subplot(1,3,2);
bar(1:M, [data.mean_dz]); hold on;
errorbar(1:M, [data.mean_dz], [data.std_dz], 'k.', 'LineWidth', 1.2);
set(gca, 'XTickLabel', {data.label}, 'XTickLabelRotation', 20, 'FontSize', 7.5);
ylabel('\delta_z [cm]'); title('Resolución vertical media'); grid on; hold off;

subplot(1,3,3);
bar(1:M, [data.psi0_deg]); hold on;
yline(theta_B_air*180/pi, ':', 'Color', [0.88 0.51 0.08], 'LineWidth', 1.5, 'DisplayName', '\theta_{B,aire}');
set(gca, 'XTickLabel', {data.label}, 'XTickLabelRotation', 20, 'FontSize', 7.5);
ylabel('\psi_0 [deg]'); title('Ángulo de look del Tx'); legend('Location','best'); grid on; hold off;

sgtitle('Medias +/- desviación estándar', 'FontWeight', 'bold');
fig2_path = fullfile(outputDir, sprintf('comparacion_barras_%s.jpg', timestamp));
saveas(fig2, fig2_path);
fprintf('Figura (barras) guardada en: %s\n', fig2_path);

%% ==================== Figura 3: trayectorias 3D superpuestas ====================
tg_plot = data(1).targets;

fig3 = figure('Position', [100, 100, 1500, 700]);
set(fig3, 'WindowState', 'maximized');

ax1 = subplot(1,2,1); hold(ax1, 'on');
maxExtent = 0;
for k = 1:M
    colTx = palette(mod(k-1, size(palette,1))+1, :);
    colRx = colTx;   % mismo color que Tx -- se distinguen por el estilo de línea (- vs --)
    plot3(ax1, data(k).Tx(1,:), data(k).Tx(2,:), data(k).Tx(3,:), '-o', 'Color', colTx, ...
        'LineWidth', 1.1, 'MarkerSize', 2.5, 'DisplayName', sprintf('Tx: %s (\\psi_0=%.1f°)', data(k).label, data(k).psi0_deg));
    plot3(ax1, data(k).Rx(1,:), data(k).Rx(2,:), data(k).Rx(3,:), '--s', 'Color', colRx, ...
        'LineWidth', 0.9, 'MarkerSize', 2.5, 'DisplayName', sprintf('Rx: %s', data(k).label));
    localExtent = max([abs(data(k).Tx(1:2,:)), abs(data(k).Rx(1:2,:))], [], 'all');
    maxExtent = max(maxExtent, localExtent);
end
plot3(ax1, tg_plot(1,:), tg_plot(2,:), tg_plot(3,:), 'd', 'Color', 'black', 'MarkerSize', 4, ...
    'LineWidth', 1, 'MarkerFaceColor', 'red', 'DisplayName', 'Blancos subsuperficiales');

[Xs, Ys] = meshgrid(-maxExtent:10:maxExtent, -maxExtent:10:maxExtent);
surf(ax1, Xs, Ys, zeros(size(Xs)), 'FaceAlpha', 0.20, 'FaceColor', [0.8 0.8 0.6], ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');

xlabel(ax1, 'X [m]'); ylabel(ax1, 'Y [m]'); zlabel(ax1, 'Z [m]');
title(ax1, 'Vista 3D'); legend(ax1, 'Location', 'northeast', 'FontSize', 7);
grid(ax1, 'on'); view(ax1, 45, 20); axis(ax1, 'tight'); box(ax1, 'on'); hold(ax1, 'off');

ax2 = subplot(1,2,2); hold(ax2, 'on');
for k = 1:M
    colTx = palette(mod(k-1, size(palette,1))+1, :);
    colRx = colTx;   % mismo color que Tx -- se distinguen por el estilo de línea (- vs --)
    plot(ax2, data(k).Tx(1,:), data(k).Tx(2,:), '-o', 'Color', colTx, 'LineWidth', 1.1, ...
        'MarkerSize', 2.5, 'DisplayName', sprintf('Tx: %s', data(k).label));
    plot(ax2, data(k).Rx(1,:), data(k).Rx(2,:), '--s', 'Color', colRx, 'LineWidth', 0.9, ...
        'MarkerSize', 2.5, 'DisplayName', sprintf('Rx: %s', data(k).label));
end
plot(ax2, tg_plot(1,:), tg_plot(2,:), 'd', 'Color', 'black', 'MarkerSize', 4, ...
    'LineWidth', 1, 'MarkerFaceColor', 'red', 'DisplayName', 'Blancos');

xlabel(ax2, 'X [m]'); ylabel(ax2, 'Y [m]');
title(ax2, 'Planta (vista superior)'); legend(ax2, 'Location', 'best', 'FontSize', 7);
grid(ax2, 'on'); axis(ax2, 'equal'); box(ax2, 'on'); hold(ax2, 'off');

sgtitle('Trayectorias Tx/Rx -- comparación de enfoques', 'FontWeight', 'bold');
fig3_path = fullfile(outputDir, sprintf('comparacion_3d_%s.jpg', timestamp));
saveas(fig3, fig3_path);
fprintf('Figura (3D + planta) guardada en: %s\n', fig3_path);

%% ==================== Funciones locales ====================

function filepath = latestMatFile(directory, pattern)
    files = dir(fullfile(directory, pattern));
    if isempty(files)
        error(['No se encontró ningún archivo "%s" en %s.\n' ...
               'Correr primero el script correspondiente (run_plano_de_voo.m, ' ...
               'run_otimiza_tx.m u otimiza_tx_brewster_puro.m).'], pattern, directory);
    end
    [~, idx] = max([files.datenum]);
    filepath = fullfile(directory, files(idx).name);
end

function d = extractRunData(R, label, matFile)
% EXTRACTRUNDATA  Normaliza el struct Resultados de cualquiera de los tres
% scripts (run_plano_de_voo.m, run_otimiza_tx.m, otimiza_tx_brewster_puro.m)
% a un formato común para comparar.

d.label = label;
d.file  = matFile;
d.n1 = R.parametros.n1;
d.n2 = R.parametros.n2;

d.Tx = R.Tx_positions;
d.Rx = R.optimizacion.Rx_optimo;
d.targets = R.targets;
d.N = size(d.Tx, 2);

if isfield(R, 'etapa1_final') && isfield(R.etapa1_final, 'psi0_star')
    d.psi0_deg = R.etapa1_final.psi0_star;   % Tx optimizado por la Etapa 1
else
    % Tx fijo (línea base): estimar psi0 a partir de la geometría real de la hélice
    d.psi0_deg = atan(mean(vecnorm(d.Tx(1:2,:), 2, 1)) / mean(d.Tx(3,:))) * 180/pi;
end

d.pot  = R.optimizacion.potencia_maxima_dBm(:);
d.dxy  = R.optimizacion.delta_xy_m(:) * 100;
d.dz   = R.optimizacion.delta_z_m(:) * 100;
d.cost = R.optimizacion.costo_combinado(:);

d.mean_pot = mean(d.pot); d.std_pot = std(d.pot);
d.mean_dxy = mean(d.dxy); d.std_dxy = std(d.dxy);
d.mean_dz  = mean(d.dz);  d.std_dz  = std(d.dz);
d.mean_cost = mean(d.cost);

end

function T1 = fresnelT1(theta_i_rad, n1, n2)
    sin_t = min(n1/n2 * sin(theta_i_rad), 1.0);
    cos_i = cos(theta_i_rad);
    cos_t = sqrt(max(1 - sin_t^2, 0));
    r = (n2*cos_i - n1*cos_t) / (n2*cos_i + n1*cos_t);
    T1 = 1 - r^2;
end
