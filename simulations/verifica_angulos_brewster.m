clc;
clear;
close all;
%
% VERIFICA_ANGULOS_BREWSTER  Verificación visual de los ángulos de incidencia
% en cada tramo del enlace biestático Tx-Tg-Rx contra sus respectivos ángulos
% de Brewster (aire->suelo y suelo->aire, que son distintos y complementarios:
% theta_B,aire = arctan(n2/n1), theta_B,suelo = arctan(n1/n2)).
%
% Reproduce EXACTAMENTE los mismos cálculos que hacen precomputeTxData.m
% (tramo Tx->Tg) y objective_function.m (tramo Tg->Rx) -- misma refracción de
% Fermat, mismo ángulo de incidencia, mismo calculateTMcoef -- para UN par
% Tx-Rx elegido por configuración (no todos los puntos a la vez, porque un
% gráfico con todos los pares sería ilegible).
%
% Para cada tramo dibuja:
%   - el rayo REAL refractado (fuente -> punto de refracción -> Tg/Rx), y
%   - un vector de referencia PUNTEADO anclado en el mismo punto de
%     refracción, en el mismo azimut que el rayo real, apuntando en la
%     dirección que tendría un rayo exactamente en el ángulo de Brewster de
%     ese lado (aire o suelo) -- si la punta del vector punteado casi
%     coincide con el marcador real (Tx, Rx o Tg), ese tramo está muy cerca
%     de Brewster; si hay una separación visible, esa es la desviación.

addpath(genpath('tools'))
addpath(genpath('flightpath'))

%% ==================== CONFIGURACIÓN ====================
resultFile   = '';          % ruta al .mat de run_plano_de_voo.m o run_otimiza_tx.m;
                             % dejar vacío para usar el más reciente de io/plan_vuelo
                             % o io/otimiza_tx (el que sea más nuevo de los dos)
tx_idx       = 1;           % índice del par Tx-Rx a visualizar (1..N)
targetSelect = 'centroid';  % 'centroid' (usa el centroide 3D del grid de alvos,
                             % el mismo que usan buildRxSearchGeometry y
                             % calculateBistaticResolution) o un índice entero
                             % 1..Ntargets para un alvo puntual específico
showFullGrid = true;        % mostrar el resto del grid de alvos como puntos tenues
%% =========================================================

if isempty(resultFile)
    resultFile = findLatestFlightPlan();
end

fprintf('=== ARCHIVO CARGADO ===\n%s\n\n', resultFile);
S = load(resultFile, 'Resultados');
Res = S.Resultados;

n1 = Res.parametros.n1;
n2 = Res.parametros.n2;

N = size(Res.Tx_positions, 2);
if tx_idx < 1 || tx_idx > N
    error('verifica_angulos_brewster:tx_idx', 'tx_idx=%d fuera de rango (1..%d).', tx_idx, N);
end

Tx_pos = Res.Tx_positions(:, tx_idx);
Rx_pos = Res.optimizacion.Rx_optimo(:, tx_idx);
tg     = Res.targets;

if strcmpi(targetSelect, 'centroid')
    target_pos  = mean(tg, 2);
    targetLabel = sprintf('centroide 3D del grid (%d alvos)', size(tg,2));
else
    target_pos  = tg(:, targetSelect);
    targetLabel = sprintf('alvo #%d del grid', targetSelect);
end

%% Ángulos de Brewster (los dos, distintos)
theta_B_air    = atan(n2 / n1);   % 63.43° para n1=1,n2=2 -- lado AIRE (Tx->Tg y Tg->Rx en el aire)
theta_B_ground = atan(n1 / n2);   % 26.57°                 -- lado SUELO (Tg->Rx en el suelo, Tx->Tg en el suelo)

%% ============================================================
%% Tramo 1: Tx -> Tg (mismo cálculo que precomputeTxData.m)
%% ============================================================
intersec1 = calculateRefractionPointFermat(Tx_pos, target_pos, n1, n2);
ang_inc1  = calculateIncidenceAngle(Tx_pos, intersec1, target_pos);
[~,~,~,~, T1, ~] = calculateTMcoef(ang_inc1, n1, n2);
R_T = norm(Tx_pos - intersec1) + norm(intersec1 - target_pos);

%% ============================================================
%% Tramo 2: Tg -> Rx (mismo cálculo que objective_function.m)
%% ============================================================
intersec2 = calculateRefractionPointFermat(target_pos, Rx_pos, n2, n1);
ang_inc2  = calculateIncidenceAngle(target_pos, intersec2, Rx_pos);
[~,~,~,~, T2, ~] = calculateTMcoef(ang_inc2, n2, n1);
R_R = norm(target_pos - intersec2) + norm(intersec2 - Rx_pos);

%% ==================== Resumen numérico ====================
fprintf('=== Par Tx-Rx #%d de %d  |  alvo: %s ===\n', tx_idx, N, targetLabel);
fprintf('Tx = [%7.2f, %7.2f, %6.2f]\n', Tx_pos);
fprintf('Rx = [%7.2f, %7.2f, %6.2f]\n', Rx_pos);
fprintf('Tg = [%7.2f, %7.2f, %6.2f]\n', target_pos);

fprintf('\n-- Tramo Tx->Tg (incidencia medida en el AIRE) --\n');
fprintf('  Punto de refracción      : [%7.2f, %7.2f, 0.00]\n', intersec1(1:2));
fprintf('  ang_inc (real, aire)     : %6.2f°\n', ang_inc1*180/pi);
fprintf('  theta_B,aire = atan(n2/n1): %6.2f°\n', theta_B_air*180/pi);
fprintf('  |diferencia|             : %6.2f°\n', abs(ang_inc1 - theta_B_air)*180/pi);
fprintf('  T1 = %.4f  (%.3f dB)\n', T1, 10*log10(T1));

fprintf('\n-- Tramo Tg->Rx (incidencia medida en el SUELO) --\n');
fprintf('  Punto de refracción      : [%7.2f, %7.2f, 0.00]\n', intersec2(1:2));
fprintf('  ang_inc (real, suelo)    : %6.2f°\n', ang_inc2*180/pi);
fprintf('  theta_B,suelo = atan(n1/n2): %6.2f°\n', theta_B_ground*180/pi);
fprintf('  |diferencia|             : %6.2f°\n', abs(ang_inc2 - theta_B_ground)*180/pi);
fprintf('  T2 = %.4f  (%.3f dB)\n', T2, 10*log10(T2));

fprintf('\nR_T (Tx->Tg) = %.2f m   R_R (Tg->Rx) = %.2f m\n', R_T, R_R);

%% ==================== Visualización 3D ====================
col_tx   = [0.00 0.45 0.74];   % azul   -- tramo Tx->Tg
col_rx   = [0.85 0.33 0.10];   % naranja-- tramo Tg->Rx
col_brew = [0.49 0.18 0.56];   % morado -- referencias de Brewster
col_tg   = [0.80 0.10 0.10];   % rojo   -- alvo

fig = figure('Position', [100, 100, 1100, 850]);
set(fig, 'WindowState', 'maximized');
hold on;

% Plano de la interfaz aire-suelo
maxExtent = max([norm(Tx_pos(1:2)), norm(Rx_pos(1:2)), norm(target_pos(1:2))]) * 1.15 + 5;
[Xs, Ys] = meshgrid(-maxExtent:maxExtent/15:maxExtent, -maxExtent:maxExtent/15:maxExtent);
surf(Xs, Ys, zeros(size(Xs)), 'FaceAlpha', 0.20, 'FaceColor', [0.8 0.8 0.6], ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');

if showFullGrid
    plot3(tg(1,:), tg(2,:), tg(3,:), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 8, ...
        'HandleVisibility', 'off');
end

% Rayo real Tx -> intersec1 -> Tg
plot3([Tx_pos(1) intersec1(1) target_pos(1)], [Tx_pos(2) intersec1(2) target_pos(2)], ...
      [Tx_pos(3) intersec1(3) target_pos(3)], '-', 'Color', col_tx, 'LineWidth', 2.2, ...
      'DisplayName', 'Rayo real Tx \rightarrow Tg (refractado)');

% Rayo real Tg -> intersec2 -> Rx
plot3([target_pos(1) intersec2(1) Rx_pos(1)], [target_pos(2) intersec2(2) Rx_pos(2)], ...
      [target_pos(3) intersec2(3) Rx_pos(3)], '-', 'Color', col_rx, 'LineWidth', 2.2, ...
      'DisplayName', 'Rayo real Tg \rightarrow Rx (refractado)');

% Marcadores
plot3(Tx_pos(1), Tx_pos(2), Tx_pos(3), 'o', 'MarkerSize', 10, 'MarkerFaceColor', col_tx, ...
      'MarkerEdgeColor', 'k', 'DisplayName', sprintf('Tx #%d', tx_idx));
plot3(Rx_pos(1), Rx_pos(2), Rx_pos(3), 's', 'MarkerSize', 10, 'MarkerFaceColor', col_rx, ...
      'MarkerEdgeColor', 'k', 'DisplayName', sprintf('Rx #%d', tx_idx));
plot3(target_pos(1), target_pos(2), target_pos(3), 'd', 'MarkerSize', 11, 'MarkerFaceColor', col_tg, ...
      'MarkerEdgeColor', 'k', 'DisplayName', ['Tg: ' targetLabel]);
plot3(intersec1(1), intersec1(2), 0, 'x', 'Color', col_tx, 'MarkerSize', 10, 'LineWidth', 2, ...
      'HandleVisibility', 'off');
plot3(intersec2(1), intersec2(2), 0, 'x', 'Color', col_rx, 'MarkerSize', 10, 'LineWidth', 2, ...
      'HandleVisibility', 'off');

% -- Vectores de referencia de Brewster (punteados, misma longitud que el
%    tramo real correspondiente, mismo azimut, anclados en el mismo punto de
%    refracción -- si el rayo real está en Brewster, la punta punteada cae
%    sobre el marcador real) --
drawBrewsterRef(intersec1, Tx_pos,     theta_B_air,    col_brew, 'Ref. Brewster aire (63.4°-tipo)');
drawBrewsterRef(intersec1, target_pos, theta_B_ground, col_brew, 'Ref. Brewster suelo (26.6°-tipo)');
drawBrewsterRef(intersec2, target_pos, theta_B_ground, col_brew, '');
drawBrewsterRef(intersec2, Rx_pos,     theta_B_air,    col_brew, '');

xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title(sprintf(['Verificación de ángulos de Brewster -- par Tx-Rx #%d\n' ...
    '\\theta_{inc,1}(aire)=%.2f° vs \\theta_{B,aire}=%.2f° (\\Delta=%.2f°)   |   ' ...
    '\\theta_{inc,2}(suelo)=%.2f° vs \\theta_{B,suelo}=%.2f° (\\Delta=%.2f°)'], ...
    tx_idx, ang_inc1*180/pi, theta_B_air*180/pi, abs(ang_inc1-theta_B_air)*180/pi, ...
    ang_inc2*180/pi, theta_B_ground*180/pi, abs(ang_inc2-theta_B_ground)*180/pi));
legend('Location', 'northeastoutside', 'FontSize', 9);
grid on; view(35, 18); axis tight; box on;
hold off;

outputDir = fullfile('io', 'verifica_angulos_brewster');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
fig_path = fullfile(outputDir, sprintf('angulos_brewster_tx%03d_%s.jpg', tx_idx, timestamp));
saveas(fig, fig_path);
fprintf('\nFigura guardada en: %s\n', fig_path);

%% ==================== Funciones locales ====================

function drawBrewsterRef(anchorPoint, realTargetPoint, theta_B, col, dispName)
% DRAWBREWSTERREF  Vector punteado anclado en anchorPoint (punto de
% refracción), en el mismo azimut que (realTargetPoint - anchorPoint), con
% longitud igual a esa distancia real, pero apuntando en el ángulo de
% Brewster theta_B respecto a la normal en vez del ángulo real. Si el rayo
% real está en Brewster, la punta de este vector coincide con realTargetPoint.

delta   = realTargetPoint - anchorPoint;
horiz   = norm(delta(1:2));
signZ   = sign(delta(3));
if signZ == 0, signZ = 1; end
if horiz < 1e-9
    azimut = [1; 0];
else
    azimut = delta(1:2) / horiz;
end

L = norm(delta);   % misma longitud que el tramo real, para comparación directa
dir3 = [sin(theta_B)*azimut(1); sin(theta_B)*azimut(2); signZ*cos(theta_B)];
tip  = anchorPoint + L * dir3;

quiver3(anchorPoint(1), anchorPoint(2), anchorPoint(3), ...
        tip(1)-anchorPoint(1), tip(2)-anchorPoint(2), tip(3)-anchorPoint(3), ...
        0, '--', 'Color', col, 'LineWidth', 1.4, 'MaxHeadSize', 0.6, ...
        'DisplayName', dispName, 'HandleVisibility', mat2handlevis(dispName));
end

function s = mat2handlevis(dispName)
if isempty(dispName)
    s = 'off';
else
    s = 'on';
end
end

function filepath = findLatestFlightPlan()
candidates = [dir(fullfile('io', 'plan_vuelo', 'optimizacion_*.mat')); ...
              dir(fullfile('io', 'otimiza_tx', 'otimiza_tx_*.mat'))];
if isempty(candidates)
    error(['No se encontró ningún plan de vuelo en io/plan_vuelo/ ni io/otimiza_tx/.\n' ...
           'Correr primero run_plano_de_voo.m o run_otimiza_tx.m, o fijar resultFile a mano.']);
end
[~, idx] = max([candidates.datenum]);
filepath = fullfile(candidates(idx).folder, candidates(idx).name);
end
