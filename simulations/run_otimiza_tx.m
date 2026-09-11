clc;
clear;
close all;
%
% RUN_OTIMIZA_TX  Implementa o Algoritmo Consolidado (Cap. 6) de
% doc/otimiza_tx.tex: Etapa 1 (otimização do Tx) seguida da Etapa 2
% (otimização do Rx dado o Tx ótimo) -- as duas plataformas, em um único
% plano de voo biestático.
%
% ETAPA 1 (executada uma vez, antes do voo): para cada psi0 candidato
% (ângulo de look do Tx), constrói-se a hélice cônica real que satisfaz a
% condição beta=psi0 (Seção 2.3) com B_helix fixo, e mede-se T1 e delta_z
% com geometria REAL (refração de Fermat ponto-a-ponto, grid completo de
% alvos, várias posições ao longo da hélice) -- não a aproximação analítica
% de primeira ordem (theta_i,1 ~ psi0, alvo único no centroide) usada em
% utils/generate_figures_tx.py para as figuras do documento teórico. O
% resultado é psi0* (busca de seção áurea) para el escenario de restrição
% elegido (strOpt.scenarioFinal), que fixa rho0,z0,beta,B_helix da hélice
% final do Tx.
%
% ETAPA 2 (executada a cada pulso do Tx, durante o voo): com o Tx já fixo
% pela Etapa 1, roda-se exatamente a mesma metodologia de
% run_plano_de_voo.m (objective_function + fmincon, 5 arranques guiados por
% Brewster) para otimizar a posição (x,y) do Rx em cada posição decimada do
% Tx -- sem nenhuma alteração em relação ao que já está validado ali.
%
% Saída: Resultados.mat + CSVs em io/otimiza_tx/, com a varredura completa
% de psi0 da Etapa 1, o ótimo psi0* por cenário/gamma, e a trajetória
% biestática completa {Tx_k, Rx_k*} da Etapa 2. Usar verifica_otimiza_tx.m
% para contrastar a Etapa 1 (simulação) contra a curva teórica de forma
% fechada.

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

optJSON    = json2struct(strcat('parametros',filesep,'otimiza_tx.json'));
strOpt     = optJSON.otimiza_tx; clear optJSON;

%% Grid de alvos subsuperficiales (mismo grid que run_plano_de_voo.m)
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

theta_B = atan(n2 / n1);
fprintf('Ángulo de Brewster: theta_B = %.2f°\n', theta_B*180/pi);
fprintf('n1=%.1f, n2=%.1f, f0=%.0f MHz, B=%.0f MHz, lambda=%.3f m\n', ...
    n1, n2, strRadarTx.FreqPortadora/1e6, B/1e6, lambda);

%% Geometría de referencia (hélice actualmente en uso)
rho0_nom      = (strRadarTx.RaioMenorEsp + strRadarTx.RaioMaiorEsp) / 2;
z0_nom        = (strRadarTx.AltMenorEsp  + strRadarTx.AltMaiorEsp)  / 2;
psi0_nom      = atan(rho0_nom / z0_nom);
R0_nom        = hypot(rho0_nom, z0_nom);
Delta_rho_nom = strRadarTx.RaioMaiorEsp - strRadarTx.RaioMenorEsp;
Delta_z_nom   = strRadarTx.AltMaiorEsp  - strRadarTx.AltMenorEsp;
B_helix       = hypot(Delta_rho_nom, Delta_z_nom);   % abertura fixa em toda a varredura

fprintf('\nSistema atual: psi0=%.2f°, R0=%.2f m, z0=%.2f m, B_helix=%.2f m\n', ...
    psi0_nom*180/pi, R0_nom, z0_nom, B_helix);

%% Varredura de psi0, simulando geometria REAL por cenário de restrição
psi0_sweep_deg = strOpt.psiMinDeg:strOpt.psiStepDeg:strOpt.psiMaxDeg;
Npsi      = numel(psi0_sweep_deg);
scenarios = strOpt.scenarios;
Nscen     = numel(scenarios);
Nsamples  = strOpt.azimuthSamples;

T1_sim      = nan(Nscen, Npsi);
delta_z_sim = nan(Nscen, Npsi);
rho0_sweep  = nan(Nscen, Npsi);
z0_sweep    = nan(Nscen, Npsi);

fprintf('\n=== VARREDURA DE psi0 (geometria real, %d amostras/hélice, %d alvos) ===\n', ...
    Nsamples, size(tg,2));

for s = 1:Nscen
    scenario = scenarios{s};
    fprintf('\n-- Escenario: %s --\n', scenario);
    for k = 1:Npsi
        [T1_sim(s,k), delta_z_sim(s,k), rho0_sweep(s,k), z0_sweep(s,k)] = simulateTxDesign( ...
            psi0_sweep_deg(k), scenario, R0_nom, z0_nom, B_helix, strRadarTx, ...
            tg, target_center_2d, n1, n2, B, lambda, Nsamples);
        if mod(k, 20) == 0 || k == Npsi
            fprintf('  psi0=%5.1f°  T1_sim=%.4f  delta_z_sim=%.2f cm\n', ...
                psi0_sweep_deg(k), T1_sim(s,k), delta_z_sim(s,k)*100);
        end
    end
end

%% Costo J_Tx(psi0;gamma) simulado, para el gamma de referencia
gammaRef  = strOpt.gammaRef;
J_sim_ref = -(1 - gammaRef) * log10(T1_sim) + gammaRef * log10(delta_z_sim);

%% Óptimo por escenario (busca de sección áurea sobre la simulación real)
tolDeg = strOpt.goldenSectionTolDeg;
psi0_star_sim  = nan(1, Nscen);
T1_star_sim    = nan(1, Nscen);
dz_star_sim    = nan(1, Nscen);

fprintf('\n=== ÓPTIMO psi0* POR ESCENARIO (gamma=%.2f, búsqueda de sección áurea) ===\n', gammaRef);
for s = 1:Nscen
    scenario = scenarios{s};
    costFun = @(psi0_deg) costTxSim(psi0_deg, scenario, R0_nom, z0_nom, B_helix, ...
        strRadarTx, tg, target_center_2d, n1, n2, B, lambda, Nsamples, gammaRef);
    psi0_star_sim(s) = goldenSectionTx(costFun, strOpt.psiMinDeg, strOpt.psiMaxDeg, tolDeg);
    [T1_star_sim(s), dz_star_sim(s)] = simulateTxDesign(psi0_star_sim(s), scenario, ...
        R0_nom, z0_nom, B_helix, strRadarTx, tg, target_center_2d, n1, n2, B, lambda, Nsamples);
    fprintf('  %-10s  psi0*_sim = %6.2f°   T1=%.4f (%.3f dB)   delta_z=%.2f cm\n', ...
        scenario, psi0_star_sim(s), T1_star_sim(s), 10*log10(T1_star_sim(s)), dz_star_sim(s)*100);
end

%% Barrido de gamma (psi0*(gamma) simulado), para el escenario R0fixed y z0fixed
gammaSweep    = strOpt.gammaSweep(:)';
Ngamma        = numel(gammaSweep);
psi0_star_gamma_sim = nan(Nscen, Ngamma);

fprintf('\n=== BARRIDO DE gamma -> psi0*(gamma) (geometria real) ===\n');
for s = 1:Nscen
    scenario = scenarios{s};
    for g = 1:Ngamma
        costFun = @(psi0_deg) costTxSim(psi0_deg, scenario, R0_nom, z0_nom, B_helix, ...
            strRadarTx, tg, target_center_2d, n1, n2, B, lambda, Nsamples, gammaSweep(g));
        psi0_star_gamma_sim(s,g) = goldenSectionTx(costFun, strOpt.psiMinDeg, strOpt.psiMaxDeg, tolDeg);
    end
    fprintf('  %-10s  gamma=[%s]\n', scenario, sprintf('%.1f ', gammaSweep));
    fprintf('  %-10s  psi0*  =[%s]\n', '', sprintf('%.1f ', psi0_star_gamma_sim(s,:)));
end

%% ============================================================
%% ETAPA 1 -> ETAPA 2: construir la hélice final del Tx
%% ============================================================
scenarioFinal = strOpt.scenarioFinal;
idxFinal      = find(strcmp(scenarios, scenarioFinal), 1);
if isempty(idxFinal)
    error('run_otimiza_tx:scenarioFinal', ...
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

Delta_rho_star = B_helix * cos(psi0_star_rad);   % condición beta=psi0 (Sección 2.3)
Delta_z_star   = B_helix * sin(psi0_star_rad);

RaioMenor_star = rho0_star - Delta_rho_star/2;
RaioMaior_star = rho0_star + Delta_rho_star/2;
AltMenor_star  = z0_star   - Delta_z_star/2;
AltMaior_star  = z0_star   + Delta_z_star/2;
beta_helix     = psi0_star_rad;   % beta = psi0*, ya validado (Sección 2.3)

fprintf('\n=== ETAPA 1 -> ETAPA 2: hélice final del Tx (escenario %s) ===\n', scenarioFinal);
fprintf('psi0* = %.2f°  (theta_B=%.2f°, 45° tomográfico)\n', psi0_star_final, theta_B*180/pi);
fprintf('rho0=%.2f m, z0=%.2f m, R0=%.2f m\n', rho0_star, z0_star, hypot(rho0_star, z0_star));
fprintf('RaioMenor=%.2f m, RaioMaior=%.2f m, AltMenor=%.2f m, AltMaior=%.2f m\n', ...
    RaioMenor_star, RaioMaior_star, AltMenor_star, AltMaior_star);

%% ============================================================
%% ETAPA 2: optimización del Rx dado el Tx óptimo (metodología de
%% run_plano_de_voo.m, sin alteraciones -- Sección 5.2 de otimiza_tx.tex)
%% ============================================================
[PxT_, PyT_, PzT_, t] = funcao_espiral(strRadarTx.NumVoltasEsp, RaioMaior_star, RaioMenor_star, ...
    AltMaior_star, AltMenor_star, strRadarTx.Vt, strRadarTx.PRF, strRadarTx.NorthOffset);

decimationFactor = strOpt.decimationFactor;
PxT = PxT_(1:decimationFactor:end);
PyT = PyT_(1:decimationFactor:end);
PzT = PzT_(1:decimationFactor:end);
t_  = t(1:decimationFactor:end);

Tx_pos = [PxT; PyT; PzT];

Pt    = strRadarTx.PotenciaTx;
sigma = strTarget.rcs;

radPattern = createRadiationPattern(strRadarTx.AperturaElev, strRadarTx.AperturaAzimut);

alpha_res = strOpt.alphaRes;
fprintf('\nOptimizando %d posiciones Tx (decimationFactor=%d)  |  %d targets  |  alpha_res=%.2f\n', ...
    length(PxT), decimationFactor, size(tg,2), alpha_res);

lb = [-strSystem.searchRadius, -strSystem.searchRadius];
ub = [ strSystem.searchRadius,  strSystem.searchRadius];
n_starts = 5;

options = optimoptions('fmincon', 'Display', 'none', ...
    'MaxIterations', 200, 'FunctionTolerance', 1e-10, ...
    'Algorithm', 'interior-point');

Rx_opt_all     = zeros(3, length(PxT));
max_power_all  = zeros(1, length(PxT));
delta_xy_all   = zeros(1, length(PxT));
delta_z_all    = zeros(1, length(PxT));
cost_opt_all   = zeros(1, length(PxT));

for tx_idx = 1:length(PxT)

    fprintf('Tx %d/%d: [%.1f, %.1f, %.1f]...\n', tx_idx, length(PxT), PxT(tx_idx), PyT(tx_idx), PzT(tx_idx));

    Tx_pos_current = [PxT(tx_idx); PyT(tx_idx); PzT(tx_idx)];
    Rx_z           = PzT(tx_idx);

    [~, ~, Gt_current] = calculateTxGainsForTargets(Tx_pos_current, tg, radPattern);
    [R_T_pre, T1_pre]  = precomputeTxData(Tx_pos_current, tg, n1, n2);
    search_geo         = buildRxSearchGeometry(Tx_pos_current, Rx_z, tg, theta_B);

    results = cell(n_starts, 1);
    fvals   = zeros(n_starts, 1);

    for i = 1:n_starts
        x0 = chooseCandidateRxPos(i, search_geo.target_center, search_geo.targets_distance, ...
            search_geo.opposite_dir, search_geo.perp_dir, search_geo.distancia_brewster, search_geo.azimut_opuesto);
        x0 = max(lb, min(ub, x0));

        [x_opt, fval] = fmincon(...
            @(xy) objective_function(xy, tg, Tx_pos_current, Rx_z, radPattern, ...
                Gt_current, Pt, sigma, lambda, n1, n2, ...
                alpha_res, B, B_helix, beta_helix, target_center_2d, R_T_pre, T1_pre), ...
            x0, [], [], [], [], lb, ub, [], options);

        results{i} = x_opt;
        fvals(i)   = fval;
    end

    [~, best_idx]        = min(fvals);
    best_xy              = results{best_idx};
    Rx_opt_all(:,tx_idx) = [best_xy(1); best_xy(2); Rx_z];
    cost_opt_all(tx_idx) = fvals(best_idx);

    c_power = objective_function(best_xy, tg, Tx_pos_current, Rx_z, radPattern, ...
        Gt_current, Pt, sigma, lambda, n1, n2, ...
        0, B, B_helix, beta_helix, target_center_2d, R_T_pre, T1_pre);
    max_power_all(tx_idx) = 10^(-c_power);

    [delta_xy_all(tx_idx), delta_z_all(tx_idx)] = calculateBistaticResolution(...
        Tx_pos_current, Rx_opt_all(:,tx_idx), target_center_2d, n2, lambda, B, B_helix, beta_helix);

end

fprintf('\n=== RESULTADOS ETAPA 2 (alpha_res = %.2f) ===\n', alpha_res);
fprintf('%-5s  %-20s  %-20s  %-12s  %-12s  %-12s\n', ...
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
outputDir = fullfile('io', 'otimiza_tx');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

timestamp     = datestr(now, 'yyyymmdd_HHMMSS');
filename_base = fullfile(outputDir, sprintf('otimiza_tx_%s', timestamp));

% ── Etapa 1 ──
Resultados.parametros.n1        = n1;
Resultados.parametros.n2        = n2;
Resultados.parametros.lambda    = lambda;
Resultados.parametros.B         = B;
Resultados.parametros.B_helix   = B_helix;
Resultados.parametros.theta_B   = theta_B;
Resultados.parametros.psi0_nom  = psi0_nom;
Resultados.parametros.R0_nom    = R0_nom;
Resultados.parametros.z0_nom    = z0_nom;
Resultados.parametros.gammaRef  = gammaRef;
Resultados.parametros.gammaSweep = gammaSweep;
Resultados.parametros.scenarios = {scenarios{:}};
Resultados.parametros.azimuthSamples = Nsamples;
Resultados.parametros.Ntargets  = size(tg, 2);
Resultados.parametros.Pt        = Pt;
Resultados.parametros.sigma     = sigma;
Resultados.parametros.alpha_res = alpha_res;

Resultados.sweep.psi0_deg       = psi0_sweep_deg;
Resultados.sweep.T1_sim         = T1_sim;
Resultados.sweep.delta_z_sim    = delta_z_sim;
Resultados.sweep.rho0_sim       = rho0_sweep;
Resultados.sweep.z0_sim         = z0_sweep;
Resultados.sweep.J_sim_gammaRef = J_sim_ref;

Resultados.optimo.psi0_star_sim = psi0_star_sim;
Resultados.optimo.T1_star_sim   = T1_star_sim;
Resultados.optimo.delta_z_star_sim = dz_star_sim;

Resultados.gammaCurve.gamma           = gammaSweep;
Resultados.gammaCurve.psi0_star_sim   = psi0_star_gamma_sim;

% ── Etapa 1 -> Etapa 2: hélice final ──
Resultados.etapa1_final.scenario   = scenarioFinal;
Resultados.etapa1_final.psi0_star  = psi0_star_final;
Resultados.etapa1_final.rho0_star  = rho0_star;
Resultados.etapa1_final.z0_star    = z0_star;
Resultados.etapa1_final.beta_helix = beta_helix;
Resultados.etapa1_final.RaioMenor  = RaioMenor_star;
Resultados.etapa1_final.RaioMaior  = RaioMaior_star;
Resultados.etapa1_final.AltMenor   = AltMenor_star;
Resultados.etapa1_final.AltMaior   = AltMaior_star;

% ── Etapa 2: plan de voo biestático completo {Tx_k, Rx_k*} ──
Resultados.Tx_positions   = Tx_pos;
Resultados.targets        = tg;
Resultados.target_center  = target_center_2d;
Resultados.optimizacion.Rx_optimo           = Rx_opt_all;
Resultados.optimizacion.potencia_maxima_W   = max_power_all;
Resultados.optimizacion.potencia_maxima_dBm = 10*log10(max_power_all*1000);
Resultados.optimizacion.delta_xy_m          = delta_xy_all;
Resultados.optimizacion.delta_z_m           = delta_z_all;
Resultados.optimizacion.costo_combinado     = cost_opt_all;

Resultados.fecha_calculo = datestr(now);

save([filename_base '.mat'], 'Resultados');
fprintf('Datos .mat guardados en: %s.mat\n', filename_base);

% CSV con la varredura completa de la Etapa 1 (una fila por psi0 x escenario)
fid = fopen([filename_base '_etapa1_sweep.csv'], 'w');
fprintf(fid, 'escenario,psi0_deg,rho0_m,z0_m,T1_sim,delta_z_sim_m,J_sim_gammaRef\n');
for s = 1:Nscen
    for k = 1:Npsi
        fprintf(fid, '%s,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f\n', ...
            scenarios{s}, psi0_sweep_deg(k), rho0_sweep(s,k), z0_sweep(s,k), ...
            T1_sim(s,k), delta_z_sim(s,k), J_sim_ref(s,k));
    end
end
fclose(fid);
fprintf('CSV (Etapa 1, varredura) guardado en: %s\n', [filename_base '_etapa1_sweep.csv']);

% CSV con el óptimo de la Etapa 1 por escenario y la curva psi0*(gamma)
fid = fopen([filename_base '_etapa1_optimo.csv'], 'w');
fprintf(fid, 'escenario,psi0_star_sim_deg,T1_star_sim,delta_z_star_sim_m\n');
for s = 1:Nscen
    fprintf(fid, '%s,%.4f,%.6f,%.6f\n', scenarios{s}, psi0_star_sim(s), T1_star_sim(s), dz_star_sim(s));
end
fprintf(fid, '\nescenario,gamma,psi0_star_sim_deg\n');
for s = 1:Nscen
    for g = 1:Ngamma
        fprintf(fid, '%s,%.2f,%.4f\n', scenarios{s}, gammaSweep(g), psi0_star_gamma_sim(s,g));
    end
end
fclose(fid);
fprintf('CSV (Etapa 1, óptimo) guardado en: %s\n', [filename_base '_etapa1_optimo.csv']);

% CSV con el plan de voo biestático completo de la Etapa 2 (mismo formato que
% run_plano_de_voo.m, para reusar verifica_plano_de_voo.m si se desea comparar)
fid = fopen([filename_base '_etapa2_planvoo.csv'], 'w');
fprintf(fid, 'Tiempo_s,Tx_X,Tx_Y,Tx_Z,Rx_X,Rx_Y,Rx_Z,Pot_W,Pot_dBm,delta_xy_m,delta_z_m,costo_J\n');
for tx_idx = 1:length(PxT)
    fprintf(fid, '%.6f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.6e,%.3f,%.6f,%.6f,%.6f\n', ...
        t_(tx_idx), ...
        Tx_pos(1,tx_idx), Tx_pos(2,tx_idx), Tx_pos(3,tx_idx), ...
        Rx_opt_all(1,tx_idx), Rx_opt_all(2,tx_idx), Rx_opt_all(3,tx_idx), ...
        max_power_all(tx_idx), 10*log10(max_power_all(tx_idx)*1000), ...
        delta_xy_all(tx_idx), delta_z_all(tx_idx), ...
        cost_opt_all(tx_idx));
end
fclose(fid);
fprintf('CSV (Etapa 2, plan de voo) guardado en: %s\n', [filename_base '_etapa2_planvoo.csv']);

%% Visualización: Etapa 1 (T1, delta_z simulados vs. psi0) y Etapa 2 (configuración biestática)
fig = figure('Position', [100, 100, 1300, 500]);
for s = 1:Nscen
    subplot(1, Nscen, s);
    yyaxis left
    plot(psi0_sweep_deg, T1_sim(s,:), '-', 'LineWidth', 2, 'Color', [0.17 0.48 0.71]);
    ylabel('T_1 simulado (Fresnel real, media sobre grid)');
    ylim([0 1.05]);
    yyaxis right
    plot(psi0_sweep_deg, delta_z_sim(s,:)*100, '--', 'LineWidth', 2, 'Color', [0.84 0.10 0.11]);
    ylabel('\delta_z simulado [cm]');
    xline(theta_B*180/pi, ':', 'Color', [0.88 0.51 0.08], 'LineWidth', 1.5, 'DisplayName', '\theta_B');
    xline(45, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5, 'DisplayName', '45°');
    xline(psi0_star_sim(s), '-.', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.5, 'DisplayName', '\psi_0^* sim');
    xlabel('\psi_0 [deg]');
    title(sprintf('%s — simulación con geometría real', scenarios{s}), 'Interpreter', 'none');
    grid on;
end
sgtitle('Etapa 1 (Tx): T_1 y \delta_z simulados vs. \psi_0', 'FontWeight', 'bold');
fig_path = fullfile(outputDir, sprintf('otimiza_tx_etapa1_sweep_%s.jpg', timestamp));
saveas(fig, fig_path);
fprintf('Figura (Etapa 1) guardada en: %s\n', fig_path);

plot_bistatic_configuration(PxT, PyT, PzT, tg, Rx_opt_all, n1, n2, theta_B);

fprintf('\nListo. psi0*=%.2f° (%s), %d posiciones Tx/Rx optimizadas.\n', ...
    psi0_star_final, scenarioFinal, length(PxT));
fprintf('Correr verifica_otimiza_tx.m para contrastar la Etapa 1 (simulación) contra la curva teórica cerrada.\n');

%% ==================== Funciones locales ====================

function [T1_avg, delta_z_avg, rho0, z0] = simulateTxDesign(psi0_deg, scenario, R0_nom, z0_nom, ...
    B_helix, strRadarTx, tg, target_center_2d, n1, n2, B, lambda, Nsamples)
% SIMULATETXDESIGN  Construye la hélice real para un psi0 candidato (bajo la
% restricción de escenario dada) y mide T1/delta_z reales promediados sobre
% el grid completo de alvos y varias posiciones a lo largo de la hélice.

psi0_rad = deg2rad(psi0_deg);

switch scenario
    case 'R0fixed'
        rho0 = R0_nom * sin(psi0_rad);
        z0   = R0_nom * cos(psi0_rad);
    case 'z0fixed'
        z0   = z0_nom;
        rho0 = z0_nom * tan(psi0_rad);
    otherwise
        error('simulateTxDesign:escenario', 'Escenario desconocido: %s', scenario);
end

% Condición beta = psi0 (Sección 2.3 de otimiza_tx.tex), B_helix fijo
Delta_rho = B_helix * cos(psi0_rad);
Delta_z   = B_helix * sin(psi0_rad);

RaioMenor = rho0 - Delta_rho/2;
RaioMaior = rho0 + Delta_rho/2;
AltMenor  = z0   - Delta_z/2;
AltMaior  = z0   + Delta_z/2;

% Restricción práctica (Sección 4.5): geometría degenerada si el radio o la
% altura mínima se vuelven no positivos
if RaioMenor <= 1.0 || AltMenor <= 1.0
    T1_avg = NaN;
    delta_z_avg = NaN;
    return;
end

[Px_, Py_, Pz_, t] = funcao_espiral(strRadarTx.NumVoltasEsp, RaioMaior, RaioMenor, ...
    AltMaior, AltMenor, strRadarTx.Vt, strRadarTx.PRF, strRadarTx.NorthOffset);

idx = unique(round(linspace(1, numel(t), Nsamples)));

T1_pts = zeros(1, numel(idx));
dz_pts = zeros(1, numel(idx));
for j = 1:numel(idx)
    Tx_pos = [Px_(idx(j)); Py_(idx(j)); Pz_(idx(j))];

    % T1 real: refracción de Fermat + Fresnel TM, promediado sobre TODO el grid
    [~, T1_pre] = precomputeTxData(Tx_pos, tg, n1, n2);
    T1_pts(j) = mean(T1_pre);

    % delta_z real: mismo Wz ya validado para el Rx, evaluado en el ángulo de
    % look real de este punto de la hélice del Tx (no la aproximación psi0)
    [~, dz] = calculateBistaticResolution(Tx_pos, Tx_pos, target_center_2d, ...
        n2, lambda, B, B_helix, psi0_rad);
    dz_pts(j) = dz;
end

T1_avg      = mean(T1_pts);
delta_z_avg = mean(dz_pts);

end

function J = costTxSim(psi0_deg, scenario, R0_nom, z0_nom, B_helix, strRadarTx, ...
    tg, target_center_2d, n1, n2, B, lambda, Nsamples, gamma)
% COSTTXSIM  J_Tx(psi0;gamma) evaluado con la simulación de geometría real.

[T1_avg, delta_z_avg] = simulateTxDesign(psi0_deg, scenario, R0_nom, z0_nom, ...
    B_helix, strRadarTx, tg, target_center_2d, n1, n2, B, lambda, Nsamples);

if ~isfinite(T1_avg) || ~isfinite(delta_z_avg) || T1_avg <= 0 || delta_z_avg <= 0
    J = 1e6;   % penalización: geometría degenerada
else
    J = -(1 - gamma) * log10(T1_avg) + gamma * log10(delta_z_avg);
end

end

function psi0_star_deg = goldenSectionTx(costFun, psiMin, psiMax, tolDeg)
% GOLDENSECTIONTX  Búsqueda de sección áurea escalar (Algoritmo 1 de
% otimiza_tx.tex), válida porque J_Tx es unimodal en psi0 (Proposición 4.1).

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
