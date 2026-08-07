clc;
clear;
close all;

% RUN_BARRIDO_ALPHA_GRID_MULTISLAB  Corre el pipeline de optimización
% multislab (run_plano_de_voo_multislab) para todas las combinaciones de
% alpha_res y tamaño de grid de targets del estudio de sensibilidad,
% pasando los valores como argumentos de función en vez de hardcodearlos
% en el script principal.
%
% Cada combinación queda guardada en su propia carpeta timestamped bajo
% io/plan_vuelo_multislab/, con su alpha_res y su grid ya registrados
% dentro del .mat (Resultados.config). analisis_plan_vuelo_multislab.m
% escanea automáticamente todas esas carpetas para consolidar resultados.

addpath(genpath('flightpath'))

%% Preparar log de esta corrida (junto a los resultados que describe)
logDir = fullfile('io', 'plan_vuelo_multislab');
if ~exist(logDir, 'dir')
    mkdir(logDir);
end
logFile = fullfile(logDir, sprintf('barrido_log_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
diary(logFile);
cleanupDiary = onCleanup(@() diary('off'));

alphas    = [0, 0.3, 0.6, 1];
gridSizes = [30, 60, 100];

n_esc = numel(alphas) * numel(gridSizes);
k = 0;
for gridSize = gridSizes
    for alpha_res = alphas
        k = k + 1;
        fprintf('\n########## Escenario %d/%d: alpha_res = %.2f, grid = %d m ##########\n', ...
            k, n_esc, alpha_res, gridSize);
        run_plano_de_voo_multislab(alpha_res, gridSize, gridSize);
    end
end

fprintf('\nBarrido completo (%d escenarios). Ejecuta analisis_plan_vuelo_multislab.m para consolidar resultados.\n', n_esc);
fprintf('Log de esta corrida guardado en: %s\n', logFile);
