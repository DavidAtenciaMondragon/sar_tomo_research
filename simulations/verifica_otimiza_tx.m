clc;
clear;
close all;
%
% VERIFICA_OTIMIZA_TX  Contrasta la simulación de geometría real de
% run_otimiza_tx.m (Fermat punto-a-ponto, grid completo de alvos, hélice real)
% contra la aproximación analítica de forma cerrada de doc/otimiza_tx.tex
% (theta_i,1 ~ psi0, alvo único no centroide, misma fórmula usada en
% utils/generate_figures_tx.py).
%
% Carga por defecto el .mat MÁS RECIENTE de io/otimiza_tx/. Para comparar una
% corrida específica, asignar explícitamente resultFile más abajo.

addpath(genpath('tools'))
addpath(genpath('flightpath'))

%% Selección de archivo de entrada
resultDir  = fullfile('io', 'otimiza_tx');
resultFile = '';   % dejar vacío para usar el .mat más reciente de resultDir

if isempty(resultFile)
    resultFile = latestMatFile(resultDir, 'otimiza_tx_*.mat');
end

fprintf('=== ARCHIVO CARGADO ===\n%s\n\n', resultFile);
S = load(resultFile, 'Resultados');
Res = S.Resultados;

n1         = Res.parametros.n1;
n2         = Res.parametros.n2;
lambda     = Res.parametros.lambda;
B          = Res.parametros.B;
B_helix    = Res.parametros.B_helix;
theta_B    = Res.parametros.theta_B;
psi0_nom   = Res.parametros.psi0_nom;
R0_nom     = Res.parametros.R0_nom;
z0_nom     = Res.parametros.z0_nom;
gammaRef   = Res.parametros.gammaRef;
scenarios  = Res.parametros.scenarios;
Nscen      = numel(scenarios);

psi0_sweep_deg = Res.sweep.psi0_deg;
T1_sim         = Res.sweep.T1_sim;
delta_z_sim    = Res.sweep.delta_z_sim;

c = physconst('lightspeed');

%% Curva teórica de forma cerrada (misma aproximación de generate_figures_tx.py):
%   T1_theory(psi0)      = T_Fresnel(theta_i=psi0)         [Ec. otimiza_tx.tex 2.1]
%   delta_z_theory(psi0) = c / (2*Wz(psi0))                 [Ec. otimiza_tx.tex 2.2]
% evaluadas en un único punto de diseño (rho0,z0) por escenario -- sin
% refracción de Fermat, sin grid de alvos, sin promediar sobre la hélice.
Npsi = numel(psi0_sweep_deg);
T1_theory      = nan(Nscen, Npsi);
delta_z_theory = nan(Nscen, Npsi);

for s = 1:Nscen
    scenario = scenarios{s};
    for k = 1:Npsi
        psi0_rad = deg2rad(psi0_sweep_deg(k));

        [~,~,~,~, T1_theory(s,k), ~] = calculateTMcoef(psi0_rad, n1, n2);

        switch scenario
            case 'R0fixed'
                R0 = R0_nom;
            case 'z0fixed'
                R0 = z0_nom / cos(psi0_rad);
            otherwise
                error('Escenario desconocido: %s', scenario);
        end

        sin_t = min(sin(psi0_rad) / n2, 1.0);
        cos_t = sqrt(max(1 - sin_t^2, 1e-9));
        Wz = n2*B*cos_t + (c * B_helix * sin(psi0_rad) * cos(psi0_rad)) / (lambda * R0 * n2 * cos_t);
        delta_z_theory(s,k) = c / (2*Wz);
    end
end

%% Óptimo teórico psi0*(gamma=gammaRef) por sección áurea sobre la curva cerrada
tolDeg = 1e-4;
psi0_star_theory = nan(1, Nscen);
for s = 1:Nscen
    scenario = scenarios{s};
    costFun = @(psi0_deg) costTxTheory(psi0_deg, scenario, R0_nom, z0_nom, B_helix, n1, n2, B, lambda, gammaRef);
    psi0_star_theory(s) = goldenSectionTx(costFun, min(psi0_sweep_deg), max(psi0_sweep_deg), tolDeg);
end

psi0_star_sim = Res.optimo.psi0_star_sim;

%% Curva psi0*(gamma) teórica, para comparar con Res.gammaCurve.psi0_star_sim
gammaSweep = Res.gammaCurve.gamma;
Ngamma     = numel(gammaSweep);
psi0_star_gamma_theory = nan(Nscen, Ngamma);
for s = 1:Nscen
    scenario = scenarios{s};
    for g = 1:Ngamma
        costFun = @(psi0_deg) costTxTheory(psi0_deg, scenario, R0_nom, z0_nom, B_helix, n1, n2, B, lambda, gammaSweep(g));
        psi0_star_gamma_theory(s,g) = goldenSectionTx(costFun, min(psi0_sweep_deg), max(psi0_sweep_deg), tolDeg);
    end
end

%% Métricas de discrepancia simulación vs. teoría
fprintf('=== SIMULACIÓN (geometría real) vs. TEORÍA (forma cerrada) ===\n');
fprintf('theta_B = %.2f°   |   psi0 en uso actualmente = %.2f°\n\n', theta_B*180/pi, psi0_nom*180/pi);

err_T1  = nan(1, Nscen);
err_dz  = nan(1, Nscen);
for s = 1:Nscen
    dT1 = T1_sim(s,:) - T1_theory(s,:);
    ddz = (delta_z_sim(s,:) - delta_z_theory(s,:)) * 100;   % cm
    err_T1(s) = sqrt(mean(dT1.^2, 'omitnan'));
    err_dz(s) = sqrt(mean(ddz.^2, 'omitnan'));

    fprintf('-- %s --\n', scenarios{s});
    fprintf('  psi0*_sim    = %6.2f°   T1=%.4f   delta_z=%.2f cm\n', ...
        psi0_star_sim(s), Res.optimo.T1_star_sim(s), Res.optimo.delta_z_star_sim(s)*100);
    fprintf('  psi0*_theory = %6.2f°   diferencia = %.2f°\n', ...
        psi0_star_theory(s), psi0_star_sim(s) - psi0_star_theory(s));
    fprintf('  RMS(T1_sim - T1_theory)        = %.4f\n', err_T1(s));
    fprintf('  RMS(delta_z_sim - delta_z_theory) = %.2f cm\n\n', err_dz(s));
end

%% Guardar CSV comparativo
outputDir = fullfile('io', 'otimiza_tx_comparacion');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

csv_sweep = fullfile(outputDir, sprintf('comparacion_sweep_%s.csv', timestamp));
fid = fopen(csv_sweep, 'w');
fprintf(fid, 'escenario,psi0_deg,T1_sim,T1_theory,delta_T1,delta_z_sim_cm,delta_z_theory_cm,delta_delta_z_cm\n');
for s = 1:Nscen
    for k = 1:Npsi
        fprintf(fid, '%s,%.4f,%.6f,%.6f,%.6f,%.4f,%.4f,%.4f\n', ...
            scenarios{s}, psi0_sweep_deg(k), ...
            T1_sim(s,k), T1_theory(s,k), T1_sim(s,k)-T1_theory(s,k), ...
            delta_z_sim(s,k)*100, delta_z_theory(s,k)*100, (delta_z_sim(s,k)-delta_z_theory(s,k))*100);
    end
end
fclose(fid);
fprintf('CSV (varredura comparada) guardado en: %s\n', csv_sweep);

csv_optimo = fullfile(outputDir, sprintf('comparacion_optimo_%s.csv', timestamp));
fid = fopen(csv_optimo, 'w');
fprintf(fid, 'escenario,psi0_star_sim_deg,psi0_star_theory_deg,diferencia_deg,psi0_nom_deg,RMS_T1,RMS_delta_z_cm\n');
for s = 1:Nscen
    fprintf(fid, '%s,%.4f,%.4f,%.4f,%.4f,%.6f,%.4f\n', scenarios{s}, ...
        psi0_star_sim(s), psi0_star_theory(s), psi0_star_sim(s)-psi0_star_theory(s), ...
        psi0_nom*180/pi, err_T1(s), err_dz(s));
end
fprintf(fid, '\nescenario,gamma,psi0_star_sim_deg,psi0_star_theory_deg,diferencia_deg\n');
for s = 1:Nscen
    for g = 1:Ngamma
        fprintf(fid, '%s,%.2f,%.4f,%.4f,%.4f\n', scenarios{s}, gammaSweep(g), ...
            Res.gammaCurve.psi0_star_sim(s,g), psi0_star_gamma_theory(s,g), ...
            Res.gammaCurve.psi0_star_sim(s,g) - psi0_star_gamma_theory(s,g));
    end
end
fclose(fid);
fprintf('CSV (óptimo comparado) guardado en: %s\n', csv_optimo);

%% Figura 1: T1 y delta_z, simulación vs. teoría, por escenario (estilo fig01_psi0_tradeoff)
fig1 = figure('Position', [100, 100, 1300, 550]);
for s = 1:Nscen
    subplot(1, Nscen, s);
    yyaxis left
    hold on;
    plot(psi0_sweep_deg, T1_theory(s,:), '-',  'LineWidth', 1.8, 'Color', [0.17 0.48 0.71], 'DisplayName', 'T_1 teoría (forma cerrada)');
    plot(psi0_sweep_deg, T1_sim(s,:),    '--', 'LineWidth', 1.8, 'Color', [0.17 0.48 0.71], 'DisplayName', 'T_1 simulación (Fermat real)');
    ylabel('T_1'); ylim([0 1.05]);
    yyaxis right
    plot(psi0_sweep_deg, delta_z_theory(s,:)*100, '-',  'LineWidth', 1.8, 'Color', [0.84 0.10 0.11], 'DisplayName', '\delta_z teoría');
    plot(psi0_sweep_deg, delta_z_sim(s,:)*100,    '--', 'LineWidth', 1.8, 'Color', [0.84 0.10 0.11], 'DisplayName', '\delta_z simulación');
    ylabel('\delta_z [cm]');
    xline(theta_B*180/pi, ':', 'Color', [0.88 0.51 0.08], 'LineWidth', 1.3);
    xline(45, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.3);
    xline(psi0_nom*180/pi, '-.', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.3);
    xlabel('\psi_0 [deg]');
    title(sprintf('%s', scenarios{s}), 'Interpreter', 'none');
    legend('Location', 'southwest', 'FontSize', 7);
    grid on; hold off;
end
sgtitle('Simulación (geometría real) vs. teoría (forma cerrada): T_1 y \delta_z', 'FontWeight', 'bold');
fig1_path = fullfile(outputDir, sprintf('comparacion_tradeoff_%s.jpg', timestamp));
saveas(fig1, fig1_path);
fprintf('Figura guardada en: %s\n', fig1_path);

%% Figura 2: psi0*(gamma), simulación vs. teoría (estilo fig02_gamma_sweep b)
fig2 = figure('Position', [100, 100, 1300, 550]);
for s = 1:Nscen
    subplot(1, Nscen, s); hold on;
    plot(gammaSweep, psi0_star_gamma_theory(s,:), '-o', 'LineWidth', 1.8, 'Color', [0.17 0.48 0.71], 'DisplayName', 'teoría');
    plot(gammaSweep, Res.gammaCurve.psi0_star_sim(s,:), '--s', 'LineWidth', 1.8, 'Color', [0.84 0.10 0.11], 'DisplayName', 'simulación');
    yline(psi0_nom*180/pi, '-.', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.3, 'DisplayName', 'sistema atual');
    yline(theta_B*180/pi, ':', 'Color', [0.88 0.51 0.08], 'LineWidth', 1.3, 'DisplayName', '\theta_B');
    xlabel('\gamma'); ylabel('\psi_0^* [deg]');
    title(sprintf('%s', scenarios{s}), 'Interpreter', 'none');
    legend('Location', 'best', 'FontSize', 8);
    grid on; hold off;
end
sgtitle('\psi_0^*(\gamma): simulación vs. teoría', 'FontWeight', 'bold');
fig2_path = fullfile(outputDir, sprintf('comparacion_gamma_%s.jpg', timestamp));
saveas(fig2, fig2_path);
fprintf('Figura guardada en: %s\n', fig2_path);

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

function J = costTxTheory(psi0_deg, scenario, R0_nom, z0_nom, B_helix, n1, n2, B, lambda, gamma)
% COSTTXTHEORY  J_Tx(psi0;gamma) evaluado con la aproximación analítica de
% forma cerrada (misma fórmula que utils/generate_figures_tx.py).

c = physconst('lightspeed');
psi0_rad = deg2rad(psi0_deg);

[~,~,~,~, T1, ~] = calculateTMcoef(psi0_rad, n1, n2);

switch scenario
    case 'R0fixed'
        R0 = R0_nom;
    case 'z0fixed'
        R0 = z0_nom / cos(psi0_rad);
    otherwise
        error('costTxTheory:escenario', 'Escenario desconocido: %s', scenario);
end

sin_t = min(sin(psi0_rad) / n2, 1.0);
cos_t = sqrt(max(1 - sin_t^2, 1e-9));
Wz = n2*B*cos_t + (c * B_helix * sin(psi0_rad) * cos(psi0_rad)) / (lambda * R0 * n2 * cos_t);
delta_z = c / (2*Wz);

if ~isfinite(T1) || ~isfinite(delta_z) || T1 <= 0 || delta_z <= 0
    J = 1e6;
else
    J = -(1 - gamma) * log10(T1) + gamma * log10(delta_z);
end

end

function psi0_star_deg = goldenSectionTx(costFun, psiMin, psiMax, tolDeg)
gr = (sqrt(5) - 1) / 2;
a = psiMin;
b = psiMax;
cc = b - gr*(b - a);
dd = a + gr*(b - a);
fc = costFun(cc);
fd = costFun(dd);

while abs(b - a) > tolDeg
    if fc < fd
        b  = dd;
        dd = cc;
        fd = fc;
        cc = b - gr*(b - a);
        fc = costFun(cc);
    else
        a  = cc;
        cc = dd;
        fc = fd;
        dd = a + gr*(b - a);
        fd = costFun(dd);
    end
end

psi0_star_deg = (a + b) / 2;

end
