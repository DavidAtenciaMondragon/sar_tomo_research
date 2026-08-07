clc;
clear;
close all;

addpath(genpath('gs'))
addpath(genpath('proc'))
addpath(genpath('tools'))
addpath(genpath('flightpath'))
addpath(genpath('common'))

% ANALISIS_PLAN_VUELO_MULTISLAB  Compara los resultados guardados por
% run_plano_de_voo_multislab.m para distintos tamaños de grid de targets
% y distintos valores de alpha_res, con el fin de visualizar el tradeoff
% potencia-resolución y sacar conclusiones sobre el efecto de ambos
% parámetros.
%
% IMPORTANTE: costo_J NO es comparable directamente entre corridas con
% alpha_res distinto, porque J = -(1-alpha)*log10(Pr) + alpha*log10(delta_xy*delta_z)
% cambia de definición con alpha (ver objective_function_multislab.m). Las
% conclusiones de tradeoff se basan en potencia_dBm y delta_xy*delta_z,
% que sí son cantidades físicas comparables entre corridas.

%% Preparar carpeta de salida y log de esta corrida
baseDir = fullfile('io', 'plan_vuelo_multislab');
outDir  = fullfile(baseDir, 'analisis_comparativo');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

logFile = fullfile(outDir, sprintf('analisis_log_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
diary(logFile);
cleanupDiary = onCleanup(@() diary('off')); %#ok<NASGU>

%% Escanear todas las corridas disponibles
runDirs = dir(baseDir);
runDirs = runDirs([runDirs.isdir] & ~ismember({runDirs.name}, {'.', '..'}));

escenarios = struct('carpeta', {}, 'alpha_res', {}, 'grid_size_m', {}, ...
    'n_targets', {}, 'potencia_dBm', {}, 'delta_xy_m', {}, 'delta_z_m', {}, 'costo_J', {}, ...
    'fecha_num', {});

for k = 1:numel(runDirs)
    matFile = fullfile(baseDir, runDirs(k).name, 'optimizacion_multislab.mat');
    if ~isfile(matFile)
        continue
    end

    data = load(matFile, 'Resultados');
    R = data.Resultados;
    tg = R.config.targets;

    escenarios(end+1).carpeta    = runDirs(k).name;                 %#ok<SAGROW>
    escenarios(end).alpha_res    = R.config.alpha_res;
    escenarios(end).grid_size_m  = max(abs(tg(1, :)));               % semi-ancho en X (= xSize)
    escenarios(end).n_targets    = size(tg, 2);
    escenarios(end).potencia_dBm = R.log.potencia_dBm;
    escenarios(end).delta_xy_m   = R.log.delta_xy_m;
    escenarios(end).delta_z_m    = R.log.delta_z_m;
    escenarios(end).costo_J      = R.log.costo_J;
    escenarios(end).fecha_num    = datenum(R.fecha_calculo);
end

if isempty(escenarios)
    error('No se encontraron resultados (optimizacion_multislab.mat) en %s', baseDir);
end

fprintf('Se cargaron %d escenarios desde %s\n\n', numel(escenarios), baseDir);

%% Deduplicar: si hay varias corridas para el mismo (alpha_res, grid) —
% por ejemplo, un barrido anterior generado con una versión desactualizada
% del optimizador, junto con un barrido nuevo — quedarse solo con la más
% reciente por combinación, para no mezclar resultados obsoletos con
% resultados vigentes en la tabla/figuras.
claves = arrayfun(@(e) sprintf('%.4f_%.1f', e.alpha_res, round(e.grid_size_m)), ...
    escenarios, 'UniformOutput', false);
[clavesUnicas, ~, grupoIdx] = unique(claves);

keepIdx = zeros(1, numel(clavesUnicas));
for u = 1:numel(clavesUnicas)
    idxGrupo = find(grupoIdx == u);
    [~, iMasReciente] = max([escenarios(idxGrupo).fecha_num]);
    keepIdx(u) = idxGrupo(iMasReciente);

    if numel(idxGrupo) > 1
        descartadas = idxGrupo(idxGrupo ~= keepIdx(u));
        fprintf(['Nota: %d corrida(s) duplicada(s) para alpha=%.2g, grid=%.0fm; ' ...
            'se usa la más reciente (%s). Descartadas: %s\n'], ...
            numel(descartadas), escenarios(keepIdx(u)).alpha_res, escenarios(keepIdx(u)).grid_size_m, ...
            escenarios(keepIdx(u)).carpeta, strjoin({escenarios(descartadas).carpeta}, ', '));
    end
end
escenarios = escenarios(keepIdx);

fprintf('\n%d escenarios únicos (alpha_res, grid) tras deduplicar.\n\n', numel(escenarios));

%% Tabla resumen: una fila por escenario (valores promedio sobre la trayectoria)
n_esc        = numel(escenarios);
alpha_v      = [escenarios.alpha_res]';
grid_v       = round([escenarios.grid_size_m]');
ntg_v        = [escenarios.n_targets]';
pot_mean_dBm = arrayfun(@(e) mean(e.potencia_dBm), escenarios)';
dxy_mean_cm  = arrayfun(@(e) mean(e.delta_xy_m), escenarios)' * 100;
dz_mean_cm   = arrayfun(@(e) mean(e.delta_z_m), escenarios)'  * 100;
res_prod_cm2 = arrayfun(@(e) mean(e.delta_xy_m .* e.delta_z_m), escenarios)' * 1e4;
costo_mean   = arrayfun(@(e) mean(e.costo_J), escenarios)';

resumen = table(grid_v, alpha_v, ntg_v, pot_mean_dBm, dxy_mean_cm, dz_mean_cm, res_prod_cm2, costo_mean, ...
    'VariableNames', {'GridSize_m', 'Alpha', 'N_targets', 'Potencia_dBm_prom', ...
                       'Delta_xy_cm_prom', 'Delta_z_cm_prom', 'Res_prod_cm2_prom', 'Costo_J_prom'});
resumen = sortrows(resumen, {'GridSize_m', 'Alpha'});
disp(resumen);

writetable(resumen, fullfile(outDir, 'resumen_escenarios.csv'));
fprintf('\nTabla resumen guardada en: %s\n', fullfile(outDir, 'resumen_escenarios.csv'));

%% Curvas de evolución agrupadas por tamaño de grid, comparando alpha
grid_sizes = unique(grid_v);

for gi = 1:numel(grid_sizes)
    idx = find(grid_v == grid_sizes(gi));
    [~, order] = sort(alpha_v(idx));
    idx = idx(order);

    figCurvas = figure('Position', [100, 100, 900, 750]);

    subplot(3, 1, 1); hold on;
    for ii = idx'
        plot(escenarios(ii).potencia_dBm, '-o', 'MarkerSize', 3, ...
            'DisplayName', sprintf('\\alpha = %.2g', alpha_v(ii)));
    end
    xlabel('Índice de posición Tx');
    ylabel('Potencia Rx (dBm)');
    title(sprintf('Potencia recibida — grid %d m', grid_sizes(gi)));
    legend('Location', 'best');
    grid on;

    subplot(3, 1, 2); hold on;
    for ii = idx'
        plot(escenarios(ii).delta_xy_m .* escenarios(ii).delta_z_m * 1e4, '-o', 'MarkerSize', 3, ...
            'DisplayName', sprintf('\\alpha = %.2g', alpha_v(ii)));
    end
    set(gca, 'YScale', 'log');
    xlabel('Índice de posición Tx');
    ylabel('\delta_{xy} \cdot \delta_z (cm^2)');
    title(sprintf('Resolución esperada — grid %d m', grid_sizes(gi)));
    legend('Location', 'best');
    grid on;

    subplot(3, 1, 3); hold on;
    for ii = idx'
        plot(escenarios(ii).costo_J, '-o', 'MarkerSize', 3, ...
            'DisplayName', sprintf('\\alpha = %.2g', alpha_v(ii)));
    end
    xlabel('Índice de posición Tx');
    ylabel('Costo J');
    title(sprintf('Función de costo — grid %d m  (J no es comparable entre \\alpha distintos)', grid_sizes(gi)));
    legend('Location', 'best');
    grid on;

    saveas(figCurvas, fullfile(outDir, sprintf('evolucion_grid_%dm.png', grid_sizes(gi))));
end
fprintf('Curvas de evolución por grid guardadas en: %s\n', outDir);

%% Curva de tradeoff potencia-resolución (cantidades físicas, comparables entre alpha)
figTradeoff = figure('Position', [100, 100, 750, 600]);
hold on;
markers = {'-o', '-s', '-^', '-d', '-v', '-p'};

for gi = 1:numel(grid_sizes)
    idx = find(grid_v == grid_sizes(gi));
    [~, order] = sort(alpha_v(idx));
    idx = idx(order);

    marker = markers{mod(gi - 1, numel(markers)) + 1};
    plot(res_prod_cm2(idx), pot_mean_dBm(idx), marker, 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'DisplayName', sprintf('grid %d m', grid_sizes(gi)));

    for jj = idx'
        text(res_prod_cm2(jj), pot_mean_dBm(jj), sprintf('  \\alpha=%.2g', alpha_v(jj)), 'FontSize', 8);
    end
end
set(gca, 'XScale', 'log');
xlabel('Resolución promedio \delta_{xy} \cdot \delta_z (cm^2)');
ylabel('Potencia promedio (dBm)');
title('Tradeoff potencia-resolución (cada punto = un valor de \alpha)');
legend('Location', 'best');
grid on;

saveas(figTradeoff, fullfile(outDir, 'tradeoff_potencia_resolucion.png'));
fprintf('Curva de tradeoff guardada en: %s\n', fullfile(outDir, 'tradeoff_potencia_resolucion.png'));

%% Conclusiones automáticas
fprintf('\n=== CONCLUSIONES ===\n');
for gi = 1:numel(grid_sizes)
    idx = find(grid_v == grid_sizes(gi));
    if numel(idx) >= 2
        corr_pot = corr(alpha_v(idx), pot_mean_dBm(idx));
        corr_res = corr(alpha_v(idx), res_prod_cm2(idx));
        fprintf(['Grid %d m: correlación alpha-potencia = %+.2f  |  ' ...
                 'correlación alpha-resolución(prod) = %+.2f\n'], ...
            grid_sizes(gi), corr_pot, corr_res);
    end
end

if numel(grid_sizes) >= 2
    alpha_ref = min(alpha_v);
    idx_ref   = grid_v > 0 & alpha_v == alpha_ref;
    if nnz(idx_ref) >= 2
        [grid_ref_sorted, order] = sort(grid_v(idx_ref));
        pot_ref_sorted = pot_mean_dBm(idx_ref);
        pot_ref_sorted = pot_ref_sorted(order);
        fprintf('\nEfecto del tamaño de grid en potencia (alpha = %.2g):\n', alpha_ref);
        for jj = 1:numel(grid_ref_sorted)
            fprintf('  grid %d m -> %.1f dBm\n', grid_ref_sorted(jj), pot_ref_sorted(jj));
        end
    end
end

fprintf(['\nNota: costo_J no debe usarse para comparar la calidad de la solución entre\n' ...
    'corridas con alpha_res distinto (cambia de definición con alpha). El tradeoff\n' ...
    'real de diseño se lee en la figura tradeoff_potencia_resolucion.png y en las\n' ...
    'columnas Potencia_dBm_prom / Res_prod_cm2_prom de la tabla resumen.\n']);

fprintf('\nLog de esta corrida guardado en: %s\n', logFile);
