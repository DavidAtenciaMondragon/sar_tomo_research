clc;
clear;
close all;
%
% VERIFICA_ENERGIA_SNELL_MONO  Comparacion fisicamente objetiva de la
% potencia biestatica recibida (ecuacion de radar completa: ganancias de
% antena, transmitancia de Fresnel/TM en ambos tramos via refraccion de
% Fermat, y perdida de propagacion 1/R^2 en cada tramo) entre un
% escenario monoestatico real (Tx y Rx con geometria de helice IDENTICA,
% NorthOffset=0 -> trayectorias coincidentes) y un escenario especular
% con la MISMA geometria de helice pero Rx desfasado 180 grados en
% azimut (NorthOffset=180).
%
% A diferencia de verifica_energia_snell_offsets.m (donde Tx y Rx tenian
% geometrias de helice distintas por diseno: Tx optimizado vs Rx
% original), aqui Tx y Rx comparten radios, alturas y Vt exactamente, asi
% que el unico grado de libertad entre escenarios es el NorthOffset del
% Rx. Esto aisla el efecto puro de mover el Rx a un azimut opuesto,
% manteniendo la razon de aspecto (rho, z) identica a la del Tx en todo
% momento.
%
%   Pr = (Pt * Gt * Gr * sigma * T1 * T2 * lambda^2) / ((4*pi)^3 * R_T^2 * R_R^2)

addpath(genpath('flightpath'))
addpath(genpath('gs'))
addpath(genpath('tools'))

%% ==================== CONFIGURACION ====================
scenarios = {'_espiral_mono0', '_espiral_mono180'};
scenarioLabels = {'monoestatico (offset=0, Tx=Rx)', 'especular (offset=180, misma geometria)'};
decimationFactor = 10;
%% =========================================================

outDir = fullfile('io','snell','processed');
figDir = fullfile('io','snell','figure');
if ~exist(outDir,'dir'); mkdir(outDir); end
if ~exist(figDir,'dir'); mkdir(figDir); end

Nscen = numel(scenarios);
results = cell(1, Nscen);

for s = 1:Nscen
    paramSuffix = scenarios{s};
    fprintf('\n=== Escenario: %s ===\n', paramSuffix);

    radarJSON  = json2struct(fullfile('parametros', sprintf('radarTx%s.json', paramSuffix)));
    strRadarTx = radarJSON.radar; clear radarJSON;

    radarJSON  = json2struct(fullfile('parametros', sprintf('radarRx%s.json', paramSuffix)));
    strRadarRx = radarJSON.radar; clear radarJSON;

    targetJSON = json2struct(fullfile('parametros', sprintf('target%s.json', paramSuffix)));
    strTarget  = targetJSON.target; clear targetJSON;

    systemJSON = json2struct(fullfile('parametros', sprintf('system%s.json', paramSuffix)));
    strSystem  = systemJSON.system; clear systemJSON;

    %% Trayectorias (misma convencion de argumentos que GS_snell_test_script.m)
    [PxT_, PyT_, PzT_] = funcao_espiral(strRadarTx.NumVoltasEsp, strRadarTx.RaioMenorEsp, strRadarTx.RaioMaiorEsp, ...
        strRadarTx.AltMaiorEsp, strRadarTx.AltMenorEsp, strRadarTx.Vt, strRadarTx.PRF, strRadarTx.NorthOffset);
    [PxR_, PyR_, PzR_] = funcao_espiral(strRadarRx.NumVoltasEsp, strRadarRx.RaioMenorEsp, strRadarRx.RaioMaiorEsp, ...
        strRadarRx.AltMaiorEsp, strRadarRx.AltMenorEsp, strRadarRx.Vt, strRadarRx.PRF, strRadarRx.NorthOffset);

    N = min(numel(PxT_), numel(PxR_));
    idx = 1:decimationFactor:N;

    PxT = PxT_(idx); PyT = PyT_(idx); PzT = PzT_(idx);
    PxR = PxR_(idx); PyR = PyR_(idx); PzR = PzR_(idx);
    Nk  = numel(idx);
    t_k = (idx-1) / strRadarTx.PRF;

    %% Parametros fisicos
    n1     = 1;
    n2     = strSystem.IndiceRefracaoSolo;
    c      = strSystem.VelocidadeLuz;
    lambda = c / strRadarTx.FreqPortadora;
    B      = strRadarTx.FreqMayor - strRadarTx.FreqMenor;
    Pt     = strRadarTx.PotenciaTx;
    sigma  = strTarget.rcs;
    target_pos = strTarget.pos(:);

    radPatternTx = createRadiationPattern(strRadarTx.AperturaElev, strRadarTx.AperturaAzimut);
    radPatternRx = createRadiationPattern(strRadarRx.AperturaElev, strRadarRx.AperturaAzimut);

    fprintf('n1=%.1f, n2=%.1f, lambda=%.4f m, B=%.1f MHz, Pt=%.1f kW, sigma=%.2f m^2\n', ...
        n1, n2, lambda, B/1e6, Pt/1e3, sigma);
    fprintf('Muestras evaluadas: %d (decimacion=%d de %d totales)\n', Nk, decimationFactor, N);

    %% Recorrido pulso a pulso: ecuacion de radar biestatica completa
    Pr        = zeros(1, Nk);
    T1_series = zeros(1, Nk);
    T2_series = zeros(1, Nk);
    Gt_series = zeros(1, Nk);
    Gr_series = zeros(1, Nk);
    R_T_series = zeros(1, Nk);
    R_R_series = zeros(1, Nk);

    for k = 1:Nk
        Tx_pos = [PxT(k); PyT(k); PzT(k)];
        Rx_pos = [PxR(k); PyR(k); PzR(k)];

        % --- Tramo 1: Tx -> Target (refraccion de Fermat, aire->suelo) ---
        intersec1 = calculateRefractionPointFermat(Tx_pos, target_pos, n1, n2);
        ang_inc1  = calculateIncidenceAngle(Tx_pos, intersec1, target_pos);
        [~,~,~,~, T1, ~] = calculateTMcoef(ang_inc1, n1, n2);
        R_T = norm(Tx_pos - intersec1) + norm(intersec1 - target_pos);

        refApontRadar = [0;0;0];
        deltaTxRef = refApontRadar - Tx_pos;
        [azTxRef, elTxRef, ~] = cart2sph(deltaTxRef(1), deltaTxRef(2), deltaTxRef(3));
        azTxRef = mod(azTxRef + 2*pi, 2*pi);
        Gt = getAntennaGain(target_pos, Tx_pos, azTxRef, elTxRef, radPatternTx);

        % --- Tramo 2: Target -> Rx (refraccion de Fermat, suelo->aire) ---
        intersec2 = calculateRefractionPointFermat(target_pos, Rx_pos, n2, n1);
        ang_inc2  = calculateIncidenceAngle(target_pos, intersec2, Rx_pos);
        [~,~,~,~, T2, ~] = calculateTMcoef(ang_inc2, n2, n1);
        R_R = norm(target_pos - intersec2) + norm(intersec2 - Rx_pos);

        delta_Tg_Rx = target_pos - Rx_pos;
        [az_Rx, el_Rx, ~] = cart2sph(delta_Tg_Rx(1), delta_Tg_Rx(2), delta_Tg_Rx(3));
        az_Rx = mod(az_Rx + 2*pi, 2*pi);
        Gr = getAntennaGain(target_pos, Rx_pos, az_Rx, el_Rx, radPatternRx);

        % --- Ecuacion de radar biestatica completa ---
        Pr(k) = (Pt * Gt * Gr * sigma * T1 * T2 * lambda^2) / ((4*pi)^3 * R_T^2 * R_R^2);

        T1_series(k) = T1; T2_series(k) = T2;
        Gt_series(k) = Gt; Gr_series(k) = Gr;
        R_T_series(k) = R_T; R_R_series(k) = R_R;
    end

    Pr_dBm = 10*log10(Pr*1000);

    fprintf('Potencia recibida: max=%.2f dBm, media=%.2f dBm, min=%.2f dBm\n', ...
        max(Pr_dBm), mean(Pr_dBm), min(Pr_dBm));
    fprintf('T1 (Tx->Target): max=%.4f, media=%.4f\n', max(T1_series), mean(T1_series));
    fprintf('T2 (Target->Rx): max=%.4f, media=%.4f\n', max(T2_series), mean(T2_series));

    results{s} = struct('paramSuffix', paramSuffix, 'label', scenarioLabels{s}, ...
        't', t_k, 'Pr', Pr, 'Pr_dBm', Pr_dBm, ...
        'T1', T1_series, 'T2', T2_series, 'Gt', Gt_series, 'Gr', Gr_series, ...
        'R_T', R_T_series, 'R_R', R_R_series);

    %% Guardar CSV por escenario
    csvFile = fullfile(outDir, sprintf('energia_biestatica%s.csv', paramSuffix));
    tbl = table(t_k(:), Pr(:), Pr_dBm(:), T1_series(:), T2_series(:), Gt_series(:), Gr_series(:), ...
        R_T_series(:), R_R_series(:), ...
        'VariableNames', {'t_s','Pr_W','Pr_dBm','T1','T2','Gt','Gr','R_T_m','R_R_m'});
    writetable(tbl, csvFile);
    fprintf('Guardado: %s\n', csvFile);
end

%% ==================== Comparacion ====================
fprintf('\n=== COMPARACION FISICA (energia biestatica) ===\n');
fprintf('%-45s %12s %12s %12s\n', 'Escenario', 'Pr_max_dBm', 'Pr_media_dBm', 'Pr_min_dBm');
for s = 1:Nscen
    r = results{s};
    fprintf('%-45s %12.2f %12.2f %12.2f\n', r.label, max(r.Pr_dBm), mean(r.Pr_dBm), min(r.Pr_dBm));
end

delta_max_dB  = max(results{2}.Pr_dBm) - max(results{1}.Pr_dBm);
delta_mean_dB = mean(results{2}.Pr_dBm) - mean(results{1}.Pr_dBm);
fprintf('\nDelta (%s vs %s): max=%.2f dB, media=%.2f dB\n', ...
    scenarioLabels{2}, scenarioLabels{1}, delta_max_dB, delta_mean_dB);

%% Figura comparativa
figComp = figure;
hold on
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]};
for s = 1:Nscen
    plot(results{s}.t, results{s}.Pr_dBm, 'Color', colors{s}, 'LineWidth', 1.5, ...
        'DisplayName', scenarioLabels{s});
end
xlabel('t (s)')
ylabel('P_r (dBm)')
title('Potencia biestatica recibida: monoestatico vs especular (misma geometria de helice)')
legend('Location','best')
grid minor
hold off
saveas(figComp, fullfile(figDir, 'comparacion_energia_mono.png'));
fprintf('\nFigura guardada: %s\n', fullfile(figDir, 'comparacion_energia_mono.png'));

%% CSV resumen combinado
summaryFile = fullfile(outDir, 'comparacion_energia_mono_resumen.csv');
labelCol = scenarioLabels(:);
maxCol   = cellfun(@(r) max(r.Pr_dBm), results)';
meanCol  = cellfun(@(r) mean(r.Pr_dBm), results)';
minCol   = cellfun(@(r) min(r.Pr_dBm), results)';
T1meanCol = cellfun(@(r) mean(r.T1), results)';
T2meanCol = cellfun(@(r) mean(r.T2), results)';
summaryTbl = table(labelCol, maxCol, meanCol, minCol, T1meanCol, T2meanCol, ...
    'VariableNames', {'Escenario','Pr_max_dBm','Pr_media_dBm','Pr_min_dBm','T1_media','T2_media'});
writetable(summaryTbl, summaryFile);
fprintf('Resumen guardado: %s\n', summaryFile);
