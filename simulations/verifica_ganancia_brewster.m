clc;
clear;
close all;
%
% VERIFICA_GANANCIA_BREWSTER  Comparación energética monoestático vs biestático
% especular frente a objetos PEC enterrados (cubo, caja, cilindro, esfera, cono:
% los mismos sólidos de tools/rcs_core), barriendo el ángulo de observación para
% ver la ganancia que aporta el ángulo de Brewster en la transmitancia TM de la
% interfaz aire-suelo. Se ejecuta para cada objeto de la lista 'objetos' y al final
% compara todos (resumen en consola, CSV y figura).
%
% Escenario:
%   - Interfaz plana aire (n1=1) / suelo (n2=2). El Tx vuela sobre un ARCO de
%     radio constante r = 100 m centrado en el punto de la superficie encima del
%     cubo (0,0,0): Tx = r*[cos(eps), 0, sin(eps)], con eps = 20:5:90° (90° = nadir).
%   - Dos configuraciones de recepción:
%       * MONOESTÁTICO: Rx = Tx.
%       * BIESTÁTICO ESPECULAR: Rx especular del Tx respecto al origen,
%         Rx = r*[-cos(eps), 0, sin(eps)] (mismo arco, lado opuesto). En nadir
%         ambas configuraciones coinciden.
%   - El ángulo de incidencia real en la interfaz (theta_i, desde la normal) se
%     obtiene con la refracción de Fermat hacia el centro del cubo; es ~90° - eps
%     (el punto de refracción se desplaza < 1 m del origen). theta_i recorre
%     ~70°..0° y pasa por Brewster (atan(2) = 63.43°); el eps_B equivalente se
%     interpola de la curva continua.
%   - Objetos PEC con el centro a 2 m de profundidad. Geometría como en
%     rcs_core/rcs/geometry.py: centrados en el origen, cilindro y cono con eje z
%     (cono con la base en z = -L/2 y el vértice arriba). Las superficies curvas se
%     facetan con flecha < lambda_s/100 (mismo criterio que simular.py con PO).
%   - RCS por Óptica Física BIESTÁTICA, generalización de tools/rcs_core/rcs/po.py
%     (integral exacta por triángulo, misma normalización). Para Tx->Rx con
%     polarización TM->TM (E en el plano de incidencia xz, como en calculateTMcoef):
%
%       sqrt(sigma) = (2*sqrt(pi)/lambda) * sum_{caras iluminadas por Tx}
%                     [e_r . (n x (k_i x e_i))] * I_T(k*(i_hat + s_hat))
%
%     i_hat: dirección cubo->Tx, s_hat: dirección cubo->Rx, k_i = -i_hat. Con
%     s_hat = i_hat se reduce exactamente a po.py (factor n.i_hat, fase 2k).
%     Se valida al arrancar contra los .mat de tools/rcs_core/resultados y contra
%     la placa especular 4*pi*A^2*cos^2(theta)/lambda^2. La RCS se evalúa con la
%     longitud de onda del suelo, lambda_s = lambda0/n2, y en las direcciones de
%     los rayos refractados (theta_t), no las del aire.
%   - Ecuación radar con transmitancia de potencia TM en la bajada (aire->suelo)
%     y en la subida (suelo->aire), con la convención de distancias de
%     verifica_energia_snell_mono.m (R = tramo aire + tramo suelo):
%
%       Pr = Pt * Gt * Gr * lambda0^2 * sigma * T_down * T_up / ((4*pi)^3 * R_T^2 * R_R^2)
%
%   - Biestático: se estima además la reflexión especular de la propia superficie
%     (teoría de imágenes, plano infinito, |Gamma_TM|^2), que en Brewster se anula:
%
%       Pr_sup = Pt * Gt * Gr * lambda0^2 * |Gamma_TM|^2 / ((4*pi)^2 * (R_Tx0 + R_0Rx)^2)
%
%   - El resto de parámetros sale de los JSON *_espiral_mono0 (radarTx, radarRx,
%     system, target).
%
% Salidas: por objeto, figura principal, diagramas 3D de RCS (mono y biestático) y
% Monte Carlo; al final, un resumen comparativo. Todo en io/verifica_ganancia_brewster/
% <fecha>_objetos/ (una subcarpeta por objeto) + log en la Command Window.

addpath(genpath('flightpath'))
addpath(genpath('tools'))

%% ==================== CONFIGURACIÓN ====================
paramSuffix   = '_espiral_mono0';   % JSON de parámetros faltantes (radarTx, radarRx, system, target)
fcOverride    = 425e6;              % portadora [Hz]; [] = usar FreqPortadora del JSON (mono0 trae 10 GHz)
r_radar       = 100;                % radio del arco de Tx y Rx, centrado en (0,0,0) [m]
eps_deg       = 20:5:90;            % ángulos de elevación del Tx [°] (90 = nadir)
dEpsFino      = 0.05;               % paso de la curva continua de referencia [°]
objetos       = {'cubo', 'caja', 'cilindro', 'esfera', 'cono'};   % se simulan en este orden
dimObjetos    = struct( ...          % dimensiones [m] (mismo convenio que simular.py --dim)
    'cubo',     1.0, ...             %   arista
    'caja',     [1.0 0.5 0.3], ...   %   lx ly lz
    'cilindro', [0.3 1.0], ...       %   radio, largo (eje z)
    'esfera',   0.5, ...             %   radio
    'cono',     [0.3 1.0]);          %   radio de la base, alto (eje z, vértice arriba)
z_obj         = -2.0;               % profundidad del CENTRO del objeto [m]
yaw_deg       = 0;                  % giro del objeto alrededor de z [°] (0 = caras alineadas con x,y)
incluirPerdas = false;              % true: atenuación dieléctrica exp(-alpha*(R_suelo_T+R_suelo_R)) con TangentePerdaSolo
                                    %  (las simulaciones anteriores no la incluyen; ver log)
validarRCS    = true;               % comparar el PO de MATLAB con los .mat de tools/rcs_core
patrones360   = true;               % figuras 2 y 3: diagramas 3D de RCS (mono y biestático, uno por figura)
dPsi          = 0.05;               % resolución de los cortes 360° usados en el log [°]
dAng3D        = 0.5;                % resolución az/el de los diagramas 3D [°] (0.25 si fc es de GHz)
rangoDin      = 30;                 % rango dinámico de los diagramas 3D [dB] (radio = dBsm - (máx - rangoDin))
montecarlo    = true;               % figura 4: Monte Carlo con el objeto rotado aleatoriamente
Nmc           = 10000;              % repeticiones (orientaciones independientes) POR POSICIÓN del barrido.
                                    %  La media depende de los aciertos especulares (lóbulos ~9° a 425 MHz,
                                    %  ~0.4° a 10 GHz): el log da el IC 95% de cada media; súbelo si es ancho
Nboot         = 2000;               % remuestreos bootstrap para los IC 95%
modoRotacion  = 'xy';               % 'xy': giros aleatorios alrededor de x y de y, uniformes en [-rotMax, rotMax]
                                    % 'so3': orientación 3D uniforme (x, y, z), cuaternión aleatorio
rotMax        = 180;                % [°] con 'xy'. 180° permite cualquier orientación de cualquier objeto.
                                    %  (45° solo basta para el cubo, por su simetría de 90°: con la caja,
                                    %  el cilindro o el cono deja orientaciones fuera y sesga el resultado)
semillaMC     = 1;                  % semilla del generador aleatorio (reproducible)
guardar       = true;               % guardar CSV, PNG y log
cerrarFiguras = true;               % cerrar las figuras de cada objeto tras guardarlas (evita decenas de ventanas)
%% =========================================================

stamp   = datestr(now, 'yyyymmdd_HHMMSS');
outRoot = fullfile('io', 'verifica_ganancia_brewster', sprintf('%s_objetos', stamp));
if guardar
    if ~exist(outRoot, 'dir'); mkdir(outRoot); end
    logFile = fullfile(outRoot, sprintf('log_%s.txt', stamp));
    diary(logFile);
end

%% ==================== Parámetros ====================
radarJSON  = json2struct(fullfile('parametros', sprintf('radarTx%s.json', paramSuffix)));
strRadarTx = radarJSON.radar; clear radarJSON;
radarJSON  = json2struct(fullfile('parametros', sprintf('radarRx%s.json', paramSuffix)));
strRadarRx = radarJSON.radar; clear radarJSON;
systemJSON = json2struct(fullfile('parametros', sprintf('system%s.json', paramSuffix)));
strSystem  = systemJSON.system; clear systemJSON;
targetJSON = json2struct(fullfile('parametros', sprintf('target%s.json', paramSuffix)));
strTarget  = targetJSON.target; clear targetJSON;

c        = strSystem.VelocidadeLuz;
n1       = 1;
n2       = strSystem.IndiceRefracaoSolo;
tanDelta = strSystem.TangentePerdaSolo;
fc       = strRadarTx.FreqPortadora;
if ~isempty(fcOverride); fc = fcOverride; end
lambda0  = c / fc;
lambdaS  = lambda0 / n2;                        % longitud de onda dentro del suelo
Pt       = strRadarTx.PotenciaTx;
Gt_dBi   = strRadarTx.GananciaAnt;              % se interpreta en dBi
Gr_dBi   = strRadarRx.GananciaAnt;
Gt       = 10^(Gt_dBi/10);                      % haces apuntados a su punto de refracción -> ganancia máxima
Gr       = 10^(Gr_dBi/10);
B        = strRadarTx.FreqMayor - strRadarTx.FreqMenor;
Tsys     = strRadarTx.Temperatura;
kB       = 1.380649e-23;
Pn       = kB * Tsys * B;                       % ruido térmico (sin figura de ruido, un solo pulso)
Pn_dBm   = 10*log10(Pn*1e3);
sigmaRef = strTarget.rcs;                       % blanco puntual de las simulaciones anteriores [m^2]
alphaS   = pi * n2 * tanDelta / lambda0;        % atenuación del suelo [Np/m] (bajas pérdidas)

theta_B  = atand(n2/n1);                        % Brewster lado aire (incidencia desde la normal)
target   = [0; 0; z_obj];

fprintf('================ VERIFICA_GANANCIA_BREWSTER (%s) ================\n', stamp);
fprintf('Parámetros: %s\n', paramSuffix);
fprintf('  fc = %.3f GHz | lambda0 = %.4f m | lambda_suelo = %.4f m\n', fc/1e9, lambda0, lambdaS);
fprintf('  Pt = %.1f kW | Gt = %.1f dBi | Gr = %.1f dBi | B = %.1f MHz | T = %g K | Pn = kTB = %.2f dBm\n', ...
    Pt/1e3, Gt_dBi, Gr_dBi, B/1e6, Tsys, Pn_dBm);
fprintf('  n1 = %g | n2 = %g | tan_delta = %g (pérdidas %s)\n', n1, n2, tanDelta, ...
    ternario(incluirPerdas, 'INCLUIDAS', 'no incluidas'));
fprintf('  Tx en arco r = %g m centrado en (0,0,0), lado +x | Rx: monoestático (= Tx) y especular (lado -x)\n', r_radar);
fprintf('  Objetos PEC: %s | centro en z = %.2f m, yaw = %g° | polarización TM->TM\n', ...
    strjoin(objetos, ', '), z_obj, yaw_deg);
fprintf('  Brewster (aire): theta_B = %.2f°\n', theta_B);
fprintf('  Barrido: eps = %g:%g:%g° (%d posiciones)\n\n', eps_deg(1), eps_deg(2)-eps_deg(1), eps_deg(end), numel(eps_deg));

%% ==================== Validación del PO ====================
if validarRCS
    validarPOcontraRcsCore(fullfile('tools', 'rcs_core', 'resultados'), c);
    validarPObiestatico(1.0, lambdaS);
end

cM = [0 0.4470 0.7410]; cBi = [0.8500 0.3250 0.0980]; cB = [0.4940 0.1840 0.5560];
cS = [0.4660 0.6740 0.1880]; cGr = [0.6 0.6 0.6];
res = struct([]);                                   % resumen comparativo entre objetos

%% ==================== Bucle de objetos ====================
for iObj = 1:numel(objetos)
    nombreObj = objetos{iObj};
    dimObj    = dimObjetos.(nombreObj);
    [Vc, Fc, descObj] = mallaObjeto(nombreObj, dimObj, lambdaS);
    Rz = [cosd(yaw_deg) -sind(yaw_deg) 0; sind(yaw_deg) cosd(yaw_deg) 0; 0 0 1];
    Vc = (Rz * Vc')';
    tapaPlana = ismember(nombreObj, {'cubo', 'caja', 'cilindro'});   % cara superior plana y horizontal
    Lc  = max(max(Vc) - min(Vc));                      % mayor dimensión del objeto [m]
    outDir = fullfile(outRoot, nombreObj);
    if guardar && ~exist(outDir, 'dir'); mkdir(outDir); end

    % ---------- Figura 0: geometría simulada ----------
    fig0 = figure('Name', ['Geometría: ' nombreObj], 'Position', [80 80 1300 560]);
    tl0  = tiledlayout(fig0, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl0, sprintf('Geometría simulada: %s (PEC, %d triángulos)', descObj, size(Fc,1)));
    ax = nexttile(tl0);
    patch(ax, 'Vertices', Vc, 'Faces', Fc, 'FaceColor', [0.72 0.78 0.86], 'EdgeColor', [0.25 0.30 0.40], 'LineWidth', 0.6);
    axis(ax, 'equal'); grid(ax, 'on'); view(ax, 35, 22);
    camlight(ax, 'headlight'); lighting(ax, 'gouraud');
    xlabel(ax, 'x [m]'); ylabel(ax, 'y [m]'); zlabel(ax, 'z [m]');
    title(ax, '(a) Malla de facetas (centrada en el objeto)');
    ax = nexttile(tl0); hold(ax, 'on');
    xg = [-2 2]; zb = z_obj - 1.2;
    fill(ax, [xg fliplr(xg)], [0 0 zb zb], [0.85 0.75 0.6], 'FaceAlpha', 0.35, 'EdgeColor', 'none');
    plot(ax, xg, [0 0], 'k-', 'LineWidth', 1.5);
    kh = convhull(Vc(:,1), Vc(:,3));
    fill(ax, Vc(kh,1), Vc(kh,3) + z_obj, [0.45 0.50 0.60], 'EdgeColor', 'k');
    plot(ax, 0, z_obj, 'k+', 'MarkerSize', 10, 'LineWidth', 1.2);
    xd = max(Vc(:,1)) + 0.35;
    plot(ax, [xd xd], [0 z_obj], 'k--');
    plot(ax, xd + [-0.08 0.08], [z_obj z_obj], 'k-');
    text(ax, xd + 0.1, z_obj/2, sprintf('%g m', -z_obj), 'FontWeight', 'bold');
    text(ax, xg(1) + 0.1, 0.35, 'aire (n_1 = 1)');
    text(ax, xg(1) + 0.1, zb + 0.25, sprintf('suelo (n_2 = %g)', n2));
    axis(ax, 'equal'); xlim(ax, xg); ylim(ax, [zb 0.8]); grid(ax, 'on');
    xlabel(ax, 'x [m]'); ylabel(ax, 'z [m]');
    title(ax, '(b) Sección en el plano xz: objeto enterrado (+: centro)');
    if guardar
        exportgraphics(fig0, fullfile(outDir, sprintf('geometria_%s.png', stamp)), 'Resolution', 150);
    end

    fprintf('\n\n############################## OBJETO %d/%d: %s ##############################\n', ...
        iObj, numel(objetos), upper(descObj));
    fprintf('-- RCS del objeto en el suelo (PO, lambda_s = %.4f m, f equivalente en vacío = %.2f GHz) --\n', ...
        lambdaS, n2*fc/1e9);
    fprintf('  Malla: %d triángulos | mayor dimensión %.2f m = %.2f lambda_s\n', size(Fc,1), Lc, Lc/lambdaS);
    sigmaNadir = rcsPO(Vc, Fc, [0 0 1], [0 0 1], lambdaS);
    fprintf('  sigma(nadir) = %.2f dBsm\n', 10*log10(sigmaNadir));
    Dff = norm(max(Vc) - min(Vc));
    fprintf('  Campo lejano en el suelo: 2D^2/lambda_s = %.0f m (D = diagonal de la caja envolvente %.2f m)\n\n', ...
        2*Dff^2/lambdaS, Dff);

    %% ==================== Barrido 5° y curva continua ====================
    args = {r_radar, target, n1, n2, Vc, Fc, lambda0, lambdaS, Pt, Gt, Gr, sigmaRef, alphaS, incluirPerdas};
    S  = simularBarrido(eps_deg, args{:});
    Sf = simularBarrido(eps_deg(1):dEpsFino:eps_deg(end), args{:});
    eps_B = interp1(fliplr(Sf.th_i), fliplr(Sf.eps), theta_B);   % elevación del arco con theta_i = theta_B

    %% ==================== Log ====================
    fprintf('-- Resultados por posición (TM->TM, R = R_aire + R_suelo) --\n');
    fprintf('  sigM/PrM: monoestático | sigB/PrB: biestático especular | beta: ángulo biestático en el suelo\n');
    fprintf('  PrSup: reflexión especular de la superficie en el Rx biestático | SCR = PrB - PrSup\n');
    fprintf('%6s %7s %7s %7s %8s %9s %9s %9s %9s %8s %9s %8s %8s %8s\n', ...
        'eps°', 'th_i°', 'th_t°', 'beta°', 'T^2[dB]', 'sigM[dBsm]', 'sigB[dBsm]', 'PrM[dBm]', 'PrB[dBm]', ...
        'dPr[dB]', 'PrSup', 'SCR[dB]', 'SNRM', 'SNRB');
    for k = 1:numel(eps_deg)
        fprintf('%6.1f %7.2f %7.2f %7.2f %8.3f %9.2f %9.2f %9.2f %9.2f %8.2f %9.2f %8.2f %8.2f %8.2f\n', ...
            S.eps(k), S.th_i(k), S.th_t(k), S.beta(k), 10*log10(S.T2_M(k)), ...
            10*log10(S.sigM(k)), 10*log10(S.sigB(k)), S.PrM_dBm(k), S.PrB_dBm(k), S.PrB_dBm(k) - S.PrM_dBm(k), ...
            S.PrSup_dBm(k), S.PrB_dBm(k) - S.PrSup_dBm(k), S.PrM_dBm(k) - Pn_dBm, S.PrB_dBm(k) - Pn_dBm);
    end

    fprintf('\nChequeo geométrico: |Tx| = %.4f..%.4f m | |Rx| = %.4f..%.4f m | max |theta_i - (90-eps)| = %.3f°\n', ...
        min(S.rTx), max(S.rTx), min(S.rRx), max(S.rRx), max(abs(S.th_i - (90 - S.eps))));
    fprintf('Simetría del biestático especular: max |T_up,B - T_down| = %.2e | max |R_R - R_T| = %.2e m\n', ...
        max(abs(S.T_upB - S.T_down)), max(abs(S.R_RB - S.R_T)));
    fprintf('Brewster en el arco: eps_B = %.2f° (theta_i = %.2f°)\n', eps_B, theta_B);

    [~, kB5] = min(abs(S.th_i - theta_B));
    [~, kN5] = max(S.eps);
    [T2maxF, kT2F] = max(Sf.T2_M);
    fprintf('\n-- Transmitancia (común a ambas configuraciones por simetría) --\n');
    fprintf('  Curva continua: T^2 máx = %.4f en eps = %.2f° (theta_i = %.2f°)\n', T2maxF, Sf.eps(kT2F), Sf.th_i(kT2F));
    fprintf('  Punto de la malla de 5° más cercano a Brewster: eps = %g° (theta_i = %.2f°), T^2 = %.4f\n', ...
        S.eps(kB5), S.th_i(kB5), S.T2_M(kB5));
    fprintf('  Ganancia de transmitancia Brewster vs nadir: %+.2f dB\n', 10*log10(T2maxF / S.T2_M(kN5)));
    fprintf('  Con sigma cte = %g m^2: Pr(eps_B) - Pr(nadir) = %+.2f dB (mono y biestático coinciden aquí)\n', ...
        sigmaRef, interp1(Sf.eps, Sf.PrRef_dBm, eps_B) - S.PrRef_dBm(kN5));

    fprintf('\n-- Monoestático vs biestático especular (%s) --\n', descObj);
    fprintf('  sigma: mono %.2f..%.2f dBsm | biestático %.2f..%.2f dBsm (barrido 5°)\n', ...
        min(10*log10(S.sigM)), max(10*log10(S.sigM)), min(10*log10(S.sigB)), max(10*log10(S.sigB)));
    if tapaPlana
        fprintf('  Especular de la cara superior plana (~4*pi*A^2*cos^2(theta_t)/lambda_s^2):\n');
        fprintf('    el Rx especular la ve en TODO el barrido; el monoestático solo en nadir (lóbulo de ~%.2f°).\n', ...
            rad2deg(0.886*lambdaS/(2*Lc)));
    end
    dPr = S.PrB_dBm - S.PrM_dBm;
    fprintf('  Pr_B - Pr_M: mínimo %+.2f dB, medio %+.2f dB, máximo %+.2f dB (excluyendo nadir, donde es 0)\n', ...
        min(dPr(1:end-1)), mean(dPr(1:end-1)), max(dPr(1:end-1)));
    fprintf('  En eps = %g° (más cercano a Brewster): PrM = %.2f dBm, PrB = %.2f dBm (%+.2f dB)\n', ...
        S.eps(kB5), S.PrM_dBm(kB5), S.PrB_dBm(kB5), dPr(kB5));
    [PrBmax, kPrB] = max(Sf.PrB_dBm(1:end-1));
    fprintf('  Pr_B continuo (sin nadir) máximo: %.2f dBm en eps = %.2f° (theta_i = %.2f°)\n', ...
        PrBmax, Sf.eps(kPrB), Sf.th_i(kPrB));

    fprintf('\n-- Reflexión especular de superficie en el Rx biestático --\n');
    [SCRmax, kS] = max(Sf.PrB_dBm - Sf.PrSup_dBm);
    fprintf('  SCR máximo (continuo) = %.1f dB en eps = %.2f° (Gamma_TM -> 0 en Brewster)\n', SCRmax, Sf.eps(kS));
    fprintf('  En eps = %g°: PrSup = %.2f dBm, SCR = %.2f dB | en eps = 45°: SCR = %.2f dB\n', ...
        S.eps(kB5), S.PrSup_dBm(kB5), S.PrB_dBm(kB5) - S.PrSup_dBm(kB5), ...
        interp1(S.eps, S.PrB_dBm - S.PrSup_dBm, 45));
    fprintf('  (plano perfecto: Gamma_TM es exactamente 0 en Brewster; con rugosidad quedaría un residuo)\n');

    if incluirPerdas
        fprintf('\n  Pérdidas dieléctricas INCLUIDAS: alpha = %.2f Np/m (%.1f dB/m).\n', alphaS, 8.686*alphaS);
    else
        fprintf('\n  NOTA: pérdidas dieléctricas no incluidas. Con tan_delta = %g a %.1f GHz: alpha = %.2f Np/m\n', ...
            tanDelta, fc/1e9, alphaS);
        fprintf('        (%.1f dB/m) -> a nadir 2*%.1f m de suelo restarían %.0f dB al objeto (no a la superficie).\n', ...
            8.686*alphaS, S.R_sueT(kN5), 2*8.686*alphaS*S.R_sueT(kN5));
    end

    %% ==================== Figura ====================
    fig = figure('Name', ['verifica_ganancia_brewster: ' nombreObj], 'Position', [60 60 1500 850]);
    tl  = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Monoestático vs biestático especular | arco r = %g m, %s PEC a z = %g m, n_2 = %g, f_c = %g GHz, TM', ...
        r_radar, descObj, z_obj, n2, fc/1e9));

    % (1) Geometría
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, r_radar*cosd(0:0.5:180), r_radar*sind(0:0.5:180), ':', 'Color', cGr);
    for k = 1:numel(eps_deg)
        plot(ax, [S.Tx(1,k) S.PT(1,k) target(1)], [S.Tx(3,k) 0 target(3)], '-', 'Color', [cM 0.5]);
        plot(ax, [S.Rx(1,k) S.PR(1,k) target(1)], [S.Rx(3,k) 0 target(3)], '-', 'Color', [cBi 0.5]);
    end
    hT = plot(ax, S.Tx(1,:), S.Tx(3,:), 'v', 'Color', cM, 'MarkerFaceColor', cM);
    hR = plot(ax, S.Rx(1,:), S.Rx(3,:), '^', 'Color', cBi, 'MarkerFaceColor', cBi);
    hBw = plot(ax, [S.Tx(1,kB5) S.PT(1,kB5) target(1) S.PR(1,kB5) S.Rx(1,kB5)], ...
        [S.Tx(3,kB5) 0 target(3) 0 S.Rx(3,kB5)], '--', 'Color', cB, 'LineWidth', 2);
    xl = [-r_radar-10 r_radar+10];
    plot(ax, xl, [0 0], 'k-', 'LineWidth', 1.2);
    fill(ax, [xl fliplr(xl)], [0 0 -12 -12], [0.85 0.75 0.6], 'FaceAlpha', 0.35, 'EdgeColor', 'none');
    kh = convhull(Vc(:,1), Vc(:,3));                   % silueta del objeto en el plano xz
    fill(ax, Vc(kh,1), Vc(kh,3) + z_obj, [0.3 0.3 0.3], 'EdgeColor', 'k');
    xlim(ax, xl); ylim(ax, [-12 r_radar+10]); daspect(ax, [1 1 1]); grid(ax, 'on');
    xlabel(ax, 'x [m]'); ylabel(ax, 'z [m]');
    legend(ax, [hT hR hBw], {'Tx (= Rx mono)', 'Rx especular', sprintf('\\epsilon = %g° (~Brewster)', S.eps(kB5))}, ...
        'Location', 'southoutside', 'Orientation', 'horizontal');
    title(ax, 'Geometría (plano xz)');

    % (2) Transmitancia
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, Sf.eps, Sf.T_down, '-', 'Color', cM, 'LineWidth', 1.3, 'DisplayName', 'T_{down} (aire\rightarrowsuelo)');
    plot(ax, Sf.eps, Sf.T2_M, '-', 'Color', cBi, 'LineWidth', 1.3, 'DisplayName', 'T_{down}\cdotT_{up} (mono = biest.)');
    plot(ax, Sf.eps, Sf.Gamma2, '-', 'Color', cS, 'LineWidth', 1.3, 'DisplayName', '|\Gamma_{TM}|^2 (superficie)');
    plot(ax, S.eps, S.T_down, 'o', 'Color', cM, 'MarkerFaceColor', cM, 'HandleVisibility', 'off');
    plot(ax, S.eps, S.T2_M, 's', 'Color', cBi, 'MarkerFaceColor', cBi, 'HandleVisibility', 'off');
    xline(ax, eps_B, '--', sprintf('Brewster %.2f°', eps_B), 'Color', cB, 'LineWidth', 1.3, 'HandleVisibility', 'off');
    grid(ax, 'on'); xlabel(ax, 'Elevación en el arco \epsilon [°]  (\theta_i \approx 90° - \epsilon)');
    ylabel(ax, 'Transmitancia / reflectancia de potencia');
    title(ax, 'Coeficientes de Fresnel TM'); legend(ax, 'Location', 'east'); xlim(ax, eps_deg([1 end]));

    % (3) RCS mono vs biestática
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, Sf.eps, 10*log10(Sf.sigM), '-', 'Color', [cM 0.35], 'HandleVisibility', 'off');
    plot(ax, Sf.eps, 10*log10(Sf.sigB), '-', 'Color', [cBi 0.6], 'HandleVisibility', 'off');
    plot(ax, S.eps, 10*log10(S.sigM), 'o-', 'Color', cM, 'MarkerFaceColor', cM, 'LineWidth', 1.3, 'DisplayName', '\sigma monoestática');
    plot(ax, S.eps, 10*log10(S.sigB), 's-', 'Color', cBi, 'MarkerFaceColor', cBi, 'LineWidth', 1.3, 'DisplayName', '\sigma biestática especular');
    yline(ax, 10*log10(sigmaRef), ':', sprintf('\\sigma_{ref} = %g m^2', sigmaRef), 'HandleVisibility', 'off');
    xline(ax, eps_B, '--', 'Color', cB, 'LineWidth', 1.3, 'HandleVisibility', 'off');
    grid(ax, 'on'); xlabel(ax, '\epsilon [°]'); ylabel(ax, '\sigma [dBsm]');
    title(ax, sprintf('RCS PO: %s en el suelo (\\lambda_s = %.1f mm)', nombreObj, lambdaS*1e3));
    legend(ax, 'Location', 'southeast'); xlim(ax, eps_deg([1 end]));
    ylim(ax, [max(10*log10(Sf.sigB))-80, max(10*log10(Sf.sigB))+5]);

    % (4) Potencia recibida
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, Sf.eps, Sf.PrM_dBm, '-', 'Color', [cM 0.35], 'HandleVisibility', 'off');
    plot(ax, Sf.eps, Sf.PrB_dBm, '-', 'Color', [cBi 0.6], 'HandleVisibility', 'off');
    plot(ax, S.eps, S.PrM_dBm, 'o-', 'Color', cM, 'MarkerFaceColor', cM, 'LineWidth', 1.3, 'DisplayName', [nombreObj ', monoestático']);
    plot(ax, S.eps, S.PrB_dBm, 's-', 'Color', cBi, 'MarkerFaceColor', cBi, 'LineWidth', 1.3, 'DisplayName', [nombreObj ', biestático']);
    plot(ax, Sf.eps, Sf.PrSup_dBm, '-', 'Color', cS, 'LineWidth', 1.3, 'DisplayName', 'superficie (Rx biestático)');
    plot(ax, Sf.eps, Sf.PrRef_dBm, ':', 'Color', 'k', 'LineWidth', 1.1, 'DisplayName', sprintf('\\sigma = %g m^2', sigmaRef));
    yline(ax, Pn_dBm, 'k--', 'kTB', 'HandleVisibility', 'off');
    xline(ax, eps_B, '--', 'Color', cB, 'LineWidth', 1.3, 'HandleVisibility', 'off');
    grid(ax, 'on'); xlabel(ax, '\epsilon [°]'); ylabel(ax, 'P_r [dBm]');
    title(ax, 'Potencia recibida (ecuación radar con T_{down}T_{up})');
    legend(ax, 'Location', 'southeast'); xlim(ax, eps_deg([1 end]));
    ylim(ax, [max(Pn_dBm - 10, min(S.PrM_dBm) - 20), max([Sf.PrSup_dBm Sf.PrB_dBm]) + 10]);

    % (5) Diferencia energética biestático - monoestático
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, Sf.eps, Sf.PrB_dBm - Sf.PrM_dBm, '-', 'Color', cGr, 'DisplayName', 'continuo');
    bar(ax, S.eps, dPr, 0.6, 'FaceColor', cBi, 'EdgeColor', 'none', 'DisplayName', 'barrido 5°');
    xline(ax, eps_B, '--', 'Color', cB, 'LineWidth', 1.3, 'HandleVisibility', 'off');
    yline(ax, 0, 'k-', 'HandleVisibility', 'off');
    grid(ax, 'on'); xlabel(ax, '\epsilon [°]'); ylabel(ax, 'P_{r,B} - P_{r,M} [dB]');
    title(ax, 'Ganancia energética del biestático especular'); legend(ax, 'Location', 'northeast');
    xlim(ax, [eps_deg(1)-3 eps_deg(end)+3]);

    % (6) Relación blanco / reflexión de superficie en el Rx biestático
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, Sf.eps, Sf.PrB_dBm - Sf.PrSup_dBm, '-', 'Color', cS, 'LineWidth', 1.3, 'DisplayName', 'continuo');
    plot(ax, S.eps, S.PrB_dBm - S.PrSup_dBm, 'o', 'Color', cS, 'MarkerFaceColor', cS, 'DisplayName', 'barrido 5°');
    xline(ax, eps_B, '--', 'Color', cB, 'LineWidth', 1.3, 'HandleVisibility', 'off');
    yline(ax, 0, 'k-', 'HandleVisibility', 'off');
    grid(ax, 'on'); xlabel(ax, '\epsilon [°]'); ylabel(ax, 'P_{r,B} - P_{r,sup} [dB]');
    title(ax, ['Biestático: ' nombreObj ' vs reflexión especular de superficie']); legend(ax, 'Location', 'northeast');
    xlim(ax, eps_deg([1 end]));
    ylim(ax, [min(S.PrB_dBm - S.PrSup_dBm) - 10, min(SCRmax, 80) + 5]);

    %% ==================== Patrones de RCS en 360° ====================
    % Dentro del suelo, con lambda_s. Tx fijo en la dirección refractada de Brewster
    % (theta_t,B desde la vertical, lado +x). Dos cortes:
    %   - Plano xz: psi = ángulo desde el cenit (+z), positivo hacia +x, 0..360°.
    %       mono: Tx = Rx = psi | biestático: Tx fijo en +theta_t,B, Rx = psi
    %       (especular en psi = -theta_t,B).
    %   - Azimut: elevación fija 90° - theta_t,B, phi = 0..360° desde +x.
    %       mono: Tx = Rx = phi | biestático: Tx fijo en phi = 0, Rx = phi (especular en 180°).
    % Desde el aire solo se alcanzan direcciones dentro del cono crítico |psi| <= asin(n1/n2).
    if patrones360
        thtB  = asind(n1/n2*sind(theta_B));          % dirección refractada de Brewster en el suelo
        thCr  = asind(n1/n2);                        % cono de direcciones accesibles desde el aire
        psi   = 0:dPsi:360-dPsi;
        dXZ   = [sind(psi(:)) zeros(numel(psi),1) cosd(psi(:))];
        iTxXZ = repmat([sind(thtB) 0 cosd(thtB)], numel(psi), 1);
        sMonoXZ = 10*log10(max(rcsPO(Vc, Fc, dXZ, dXZ, lambdaS), 1e-30));
        sBiXZ   = 10*log10(max(rcsPO(Vc, Fc, iTxXZ, dXZ, lambdaS), 1e-30));

        elA   = 90 - thtB;
        dAZ   = [cosd(elA)*cosd(psi(:)) cosd(elA)*sind(psi(:)) sind(elA)*ones(numel(psi),1)];
        iTxAZ = repmat([cosd(elA) 0 sind(elA)], numel(psi), 1);
        sMonoAZ = 10*log10(max(rcsPO(Vc, Fc, dAZ, dAZ, lambdaS), 1e-30));
        sBiAZ   = 10*log10(max(rcsPO(Vc, Fc, iTxAZ, dAZ, lambdaS), 1e-30));

        % Tolerancia del lóbulo especular biestático (corte xz): ancho a -3 dB alrededor de -theta_t,B
        psiS = wrapTo180loc(psi);
        [sPk, kPk] = max(sBiXZ .* (abs(psiS(:) + thtB) < 2) - 1e3*(abs(psiS(:) + thtB) >= 2));
        k1 = kPk; while k1 > 1 && sBiXZ(k1-1) >= sPk - 3; k1 = k1 - 1; end
        k2 = kPk; while k2 < numel(psi) && sBiXZ(k2+1) >= sPk - 3; k2 = k2 + 1; end
        dPsi3 = (psi(k2) - psi(k1)) / 2;             % semiancho en el suelo [°]
        % Equivalente en el aire. Lado hacia nadir: siempre válido. Lado hacia rasante: si el lóbulo
        % pasa del ángulo crítico, cubre todas las direcciones accesibles hasta la incidencia rasante.
        dThI  = theta_B - asind(n2/n1*sind(max(thtB - dPsi3, 0)));   % [°] tolerancia hacia nadir
        haciaRasanteCompleto = (thtB + dPsi3) >= thCr;
        if haciaRasanteCompleto
            dThIr = 90 - theta_B;
        else
            dThIr = asind(n2/n1*sind(thtB + dPsi3)) - theta_B;      % [°] tolerancia hacia rasante
        end

        fprintf('\n-- Patrones de RCS en 360° (Tx fijo en la dirección refractada de Brewster, theta_t,B = %.2f°) --\n', thtB);
        fprintf('  Cono accesible desde el aire: |psi| <= asin(n1/n2) = %.1f° (fuera de él la RCS no es observable)\n', thCr);
        fprintf('  Corte xz  | mono: máx %.2f dBsm en psi = %.2f° | biestático: máx %.2f dBsm en psi = %.2f°\n', ...
            max(sMonoXZ), psiS(find(sMonoXZ == max(sMonoXZ), 1)), max(sBiXZ), psiS(find(sBiXZ == max(sBiXZ), 1)));
        fprintf('            | mono en psi = +theta_t,B: %.2f dBsm | biestático en el especular psi = -theta_t,B: %.2f dBsm\n', ...
            interp1(psiS, sMonoXZ, thtB, 'nearest'), sPk);
        fprintf('  Corte az  | mono: máx %.2f dBsm | biestático: máx %.2f dBsm en phi = %.2f°\n', ...
            max(sMonoAZ), max(sBiAZ), psi(find(sBiAZ == max(sBiAZ), 1)));
        fprintf('  El máximo biestático absoluto está en psi = %.2f°: dirección de avance (forward scatter,\n', ...
            psiS(find(sBiXZ == max(sBiXZ), 1)));
        fprintf('  sombra del objeto), apunta hacia abajo y no es observable desde el aire.\n');
        fprintf('  Lóbulo especular biestático: semiancho -3 dB = %.2f° en el suelo\n', dPsi3);
        fprintf('    -> en el aire: -%.2f° (hacia nadir) / +%.2f° (hacia rasante%s) de theta_B\n', dThI, dThIr, ...
            ternario(haciaRasanteCompleto, ', el lóbulo alcanza el ángulo crítico', ''));
        fprintf('    -> a r = %g m el Rx puede alejarse ~%.1f m hacia nadir y ~%.1f m hacia rasante del punto especular\n', ...
            r_radar, r_radar*deg2rad(dThI), r_radar*deg2rad(dThIr));
        if tapaPlana
            fprintf('  Por qué la sigma biestática del barrido es casi constante: con Rx especular la bisectriz\n');
            fprintf('  Tx-Rx es SIEMPRE vertical, así que la cara superior se ve siempre en su lóbulo especular\n');
            fprintf('  (teorema de equivalencia mono/biestático: sigma_B(beta) ~ sigma_M(bisectriz) = nadir).\n');
            fprintf('  Solo cae el factor cos^2(theta_t): de 0 a %.1f°, %.2f dB.\n', max(S.th_t), 20*log10(cosd(max(S.th_t))));
        end

        % Diagramas 3D de RCS sobre toda la esfera (convención de rcs_core: azimut desde +x en el
        % plano xy, elevación desde el plano xy). Mono: Tx = Rx en (az, el). Biestático: Tx fijo en
        % la dirección refractada de Brewster (az = 0°, el = 90° - theta_t,B) y Rx en (az, el).
        az3 = 0:dAng3D:360;
        el3 = -90:dAng3D:90;
        [AZ3, EL3] = meshgrid(az3, el3);                         % Nel x Naz
        d3  = [cosd(EL3(:)).*cosd(AZ3(:)), cosd(EL3(:)).*sind(AZ3(:)), sind(EL3(:))];
        iTx3 = [sind(thtB) 0 cosd(thtB)];
        sMono3 = zeros(size(d3,1), 1); sBi3 = sMono3;
        bloque = max(1000, floor(2e6 / size(Fc,1)));          % direcciones por bloque (memoria)
        for j0 = 1:bloque:size(d3,1)
            jj = j0:min(j0+bloque-1, size(d3,1));
            sMono3(jj) = rcsPO(Vc, Fc, d3(jj,:), d3(jj,:), lambdaS);
            sBi3(jj)   = rcsPO(Vc, Fc, repmat(iTx3, numel(jj), 1), d3(jj,:), lambdaS);
        end
        sMono3 = reshape(10*log10(max(sMono3, 1e-30)), size(AZ3));
        sBi3   = reshape(10*log10(max(sBi3,   1e-30)), size(AZ3));
        [~, kM3] = max(sMono3(:)); [~, kB3] = max(sBi3(:));
        fprintf('  Esfera 3D (paso %.2f°, %d direcciones): mono máx %.2f dBsm en az = %.1f°, el = %.1f°\n', ...
            dAng3D, numel(AZ3), sMono3(kM3), AZ3(kM3), EL3(kM3));
        fprintf('                                           biestático máx %.2f dBsm en az = %.1f°, el = %.1f°\n', ...
            sBi3(kB3), AZ3(kB3), EL3(kB3));

        sMax3 = max([sMono3(:); sBi3(:)]);
        rMin3 = sMax3 - rangoDin;
        fig3M = figure('Name', ['RCS 3D monoestática: ' nombreObj], 'Position', [80 80 1000 850]);
        diagramaRCS3D(fig3M, AZ3, EL3, sMono3, rMin3, sMax3, Vc, Fc, iTx3, [], thCr, ...
            sprintf('RCS monoestática 3D: %s en el suelo (PO, \\lambda_s = %.1f mm)', descObj, lambdaS*1e3));
        fig3B = figure('Name', ['RCS 3D biestática: ' nombreObj], 'Position', [120 80 1000 850]);
        diagramaRCS3D(fig3B, AZ3, EL3, sBi3, rMin3, sMax3, Vc, Fc, iTx3, [-sind(thtB) 0 cosd(thtB)], thCr, ...
            sprintf('RCS biestática 3D: %s (Tx fijo en \\theta_{t,B} = %.2f°, Rx en toda la esfera)', ...
            descObj, thtB));
    end

    %% ==================== Monte Carlo: cubo con orientación aleatoria ====================
    % Para CADA posición del barrido (eps = 20:5:90°) se generan Nmc orientaciones aleatorias e
    % independientes del cubo (giro alrededor de su centro). En cada orientación se calcula la
    % potencia recibida en monoestático y en biestático especular con EL MISMO cubo, así que la
    % comparación es pareada. Geometría, transmitancias y factores K de la ecuación radar son los
    % del barrido sin rotar: Pr = K * sigma. Polarización: TM en Tx y TM en Rx (co-polar).
    if montecarlo
        rng(semillaMC);
        Ne   = numel(S.eps);
        Ntot = Nmc * Ne;                                        % filas ordenadas por posición: r = (k-1)*Nmc + m
        [Rall, angMC] = orientacionesAleatorias(Ntot, modoRotacion, rotMax);

        % Vectores globales por posición: direcciones y polarizaciones TM (las mismas que usa rcsPO)
        yv  = repmat([0 1 0], Ne, 1);
        eIg = normalizaFilas(cross(yv, -S.iHat, 2));             % TM incidente
        eMg = normalizaFilas(cross(yv,  S.iHat, 2));             % TM recibida mono (s = i)
        eBg = normalizaFilas(cross(yv,  S.sHat, 2));             % TM recibida biestática
        % Cubo rotado R con geometría fija  <=>  cubo sin rotar con vectores R'*v (más rápido)
        rotT = @(Rb, v) [Rb(:,1).*v(:,1) + Rb(:,4).*v(:,2) + Rb(:,7).*v(:,3), ...
                         Rb(:,2).*v(:,1) + Rb(:,5).*v(:,2) + Rb(:,8).*v(:,3), ...
                         Rb(:,3).*v(:,1) + Rb(:,6).*v(:,2) + Rb(:,9).*v(:,3)];
        kRow = repelem((1:Ne)', Nmc, 1);                        % posición de cada fila
        [sM, sB] = deal(zeros(Ntot, 1));
        bloqueMC = max(1000, floor(3e6 / size(Fc,1)));        % filas por bloque (memoria)
        for r0 = 1:bloqueMC:Ntot
            rr = (r0:min(r0 + bloqueMC - 1, Ntot))';
            Rb = Rall(rr,:); kk = kRow(rr);
            iR = rotT(Rb, S.iHat(kk,:)); sR = rotT(Rb, S.sHat(kk,:));
            eIr = rotT(Rb, eIg(kk,:));
            sM(rr) = rcsPO(Vc, Fc, iR, iR, lambdaS, eIr, rotT(Rb, eMg(kk,:)));
            sB(rr) = rcsPO(Vc, Fc, iR, sR, lambdaS, eIr, rotT(Rb, eBg(kk,:)));
        end
        sigMmc = reshape(sM, Nmc, Ne);                          % [Nmc x Ne] m^2
        sigBmc = reshape(sB, Nmc, Ne);
        angMC  = reshape(angMC, Nmc, Ne, 2);

        % Comprobación: vectores rotados == malla rotada (giro de 17° en x y 23° en y)
        Rt  = rotY(23) * rotX(17);
        Rtb = repmat(reshape(Rt', 1, 9), Ne, 1);
        errRot = max(abs(10*log10(rcsPO(Vc, Fc, rotT(Rtb, S.iHat), rotT(Rtb, S.sHat), lambdaS, rotT(Rtb, eIg), rotT(Rtb, eBg))) ...
                         - 10*log10(rcsPO((Rt * Vc')', Fc, S.iHat, S.sHat, lambdaS))));

        % ---------- Estadísticas por posición ----------
        PrMmc = S.KM .* sigMmc;                                 % [Nmc x Ne] W
        PrBmc = S.KB .* sigBmc;
        dBm   = @(P) 10*log10(max(P, 1e-40)*1e3);
        meanM = dBm(mean(PrMmc, 1));  meanB = dBm(mean(PrBmc, 1));   % media energética (promedio en W)
        dMean = meanB - meanM;
        [icM, icB, icD] = deal(zeros(2, Ne));                   % IC 95% bootstrap (pareado)
        for k = 1:Ne
            ib = randi(Nmc, Nmc, Nboot);
            mMb = mean(PrMmc(ib + (k-1)*Nmc), 1);               % medias remuestreadas [1 x Nboot]
            mBb = mean(PrBmc(ib + (k-1)*Nmc), 1);
            icM(:,k) = pctl(dBm(mMb), [2.5 97.5]);
            icB(:,k) = pctl(dBm(mBb), [2.5 97.5]);
            icD(:,k) = pctl(10*log10(mBb ./ mMb), [2.5 97.5]);
        end
        pM = pctl(dBm(PrMmc), [10 50 90]);  pB = pctl(dBm(PrBmc), [10 50 90]);   % [3 x Ne] dBm
        pGana = mean(PrBmc > PrMmc, 1);
        % Acierto = reflexión especular: sigma a menos de 10 dB del máximo monoestático del objeto
        if patrones360
            sigMaxObj = max(sMono3(:));
        else
            [aa, ee] = meshgrid(0:1:360, -90:1:90);
            dd = [cosd(ee(:)).*cosd(aa(:)), cosd(ee(:)).*sind(aa(:)), sind(ee(:))];
            sigMaxObj = 10*log10(max(rcsPO(Vc, Fc, dd, dd, lambdaS)));
        end
        umbralAc = sigMaxObj - 10;
        hitM  = sum(10*log10(sigMmc) > umbralAc, 1);  hitB = sum(10*log10(sigBmc) > umbralAc, 1);

        fprintf('\n================ MONTE CARLO: %s CON ORIENTACIÓN ALEATORIA ================\n', upper(descObj));
        fprintf('  Polarización: TM en todo el cálculo (Fresnel TM en la interfaz y RCS TM->TM co-polar;\n');
        fprintf('                no se calcula TE ni polarización cruzada).\n');
        fprintf('  %d orientaciones independientes POR POSICIÓN (%d posiciones, %d en total), modo ''%s''', ...
            Nmc, Ne, Ntot, modoRotacion);
        if strcmp(modoRotacion, 'xy'); fprintf(': giro en x y en y uniformes en [-%g°, %g°]', rotMax, rotMax); end
        fprintf(', semilla %d.\n', semillaMC);
        fprintf('  Mono y especular usan el mismo objeto en cada repetición (comparación pareada).\n');
        fprintf('  Chequeo vectores rotados vs malla rotada: max error = %.1e dB\n\n', errRot);
        fprintf('  Cómo leer la tabla:\n');
        fprintf('    Media M / Media B : potencia media recibida (promedio en W, pasado a dBm) con su IC 95%%.\n');
        fprintf('                        Es la ENERGÍA que llega en promedio: la dominan las pocas orientaciones\n');
        fprintf('                        en las que una cara del objeto queda especular.\n');
        fprintf('    B-M               : cuánto más recibe en media el especular que el mono [dB]. Si su IC\n');
        fprintf('                        contiene el 0, no hay diferencia estadísticamente clara.\n');
        fprintf('    Mediana M / B     : el caso "típico" (la mitad de las repeticiones recibe menos que esto).\n');
        fprintf('    P(B>M)            : fracción de repeticiones en que el especular recibe más que el mono.\n');
        fprintf('    Aciertos M/B      : repeticiones con sigma > %.1f dBsm (reflexión especular: a menos de 10 dB\n', umbralAc);
        fprintf('                        del máximo monoestático del objeto); de ellas depende la media.\n\n');
        fprintf('%6s %7s | %-24s | %-24s | %-22s | %9s %9s | %6s | %9s\n', 'eps°', 'th_i°', ...
            'Media M [dBm] (IC95)', 'Media B [dBm] (IC95)', 'B-M [dB] (IC95)', 'Med. M', 'Med. B', 'P(B>M)', 'Aciert M/B');
        for k = 1:Ne
            fprintf('%6.1f %7.2f | %7.2f [%6.1f,%6.1f] | %7.2f [%6.1f,%6.1f] | %+6.2f [%+5.1f,%+5.1f] | %9.2f %9.2f | %6.2f | %4d/%4d\n', ...
                S.eps(k), S.th_i(k), meanM(k), icM(1,k), icM(2,k), meanB(k), icB(1,k), icB(2,k), ...
                dMean(k), icD(1,k), icD(2,k), pM(2,k), pB(2,k), pGana(k), hitM(k), hitB(k));
        end
        kOff  = find(S.eps < 90);                               % en nadir mono = especular
        nMas  = sum(icD(1,kOff) > 0); nMenos = sum(icD(2,kOff) < 0); nIgual = numel(kOff) - nMas - nMenos;
        fprintf('\n  CONCLUSIÓN (posiciones fuera del nadir, %d):\n', numel(kOff));
        fprintf('    - Especular significativamente MEJOR en media: %d posiciones\n', nMas);
        fprintf('    - Especular significativamente PEOR en media : %d posiciones\n', nMenos);
        fprintf('    - Sin diferencia clara (IC contiene 0)       : %d posiciones\n', nIgual);
        fprintf('    - Media de todas las posiciones: B-M = %+.2f dB | P(B>M) media = %.2f\n', ...
            10*log10(mean(PrBmc(:,kOff), 'all') / mean(PrMmc(:,kOff), 'all')), mean(pGana(kOff)));
        fprintf('    - Referencia con el objeto SIN rotar: B-M medio = %+.2f dB\n', ...
            10*log10(sum(S.PrB(kOff)) / sum(S.PrM(kOff))));
        fprintf('    - Aciertos especulares por posición: mono %d..%d, especular %d..%d de %d. Con pocos\n', ...
            min(hitM(kOff)), max(hitM(kOff)), min(hitB(kOff)), max(hitB(kOff)), Nmc);
        fprintf('      aciertos la media es ruidosa: mirar el IC antes de concluir (sube Nmc para estrecharlo).\n');
        fprintf('    - La mediana (caso típico) está en la zona no especular, donde PO puede errar 5-20 dB.\n');

        % ---------- Figura 4 ----------
        fig4 = figure('Name', ['Monte Carlo: ' nombreObj], 'Position', [60 60 1500 560]);
        tl4  = tiledlayout(fig4, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        titMC = sprintf('Monte Carlo: %d orientaciones aleatorias de OBJ POR POSICIÓN (modo ''%s''', Nmc, modoRotacion);
        if strcmp(modoRotacion, 'xy'); titMC = [titMC sprintf(', \\pm%g° en x e y', rotMax)]; end
        titMC = strrep(titMC, 'aleatorias de OBJ', ['aleatorias de ' descObj]);
        title(tl4, {[titMC ')'], sprintf('Polarización TM\\rightarrowTM | arco r = %g m, n_2 = %g, f_c = %g GHz', r_radar, n2, fc/1e9)});
        xlimMC = [eps_deg(1)-3 eps_deg(end)+3];

        % (a) Media energética con IC 95%
        ax = nexttile(tl4); hold(ax, 'on');
        errorbar(ax, S.eps - 0.6, meanM, meanM - icM(1,:), icM(2,:) - meanM, 'o-', 'Color', cM, ...
            'MarkerFaceColor', cM, 'LineWidth', 1.3, 'CapSize', 4, 'DisplayName', 'monoestático');
        errorbar(ax, S.eps + 0.6, meanB, meanB - icB(1,:), icB(2,:) - meanB, 's-', 'Color', cBi, ...
            'MarkerFaceColor', cBi, 'LineWidth', 1.3, 'CapSize', 4, 'DisplayName', 'especular');
        plot(ax, S.eps, S.PrRef_dBm, 'k:', 'LineWidth', 1.1, 'DisplayName', sprintf('\\sigma = %g m^2 (referencia)', sigmaRef));
        xline(ax, eps_B, '--', 'Brewster', 'Color', cB, 'LineWidth', 1.2, 'HandleVisibility', 'off');
        grid(ax, 'on'); xlim(ax, xlimMC); xlabel(ax, '\epsilon [°]'); ylabel(ax, 'P_r media [dBm]');
        title(ax, '(a) Potencia media recibida (energía), IC 95%'); legend(ax, 'Location', 'southwest');

        % (b) Diferencia de medias especular - mono con IC 95%
        ax = nexttile(tl4); hold(ax, 'on');
        fill(ax, [xlimMC fliplr(xlimMC)], [0 0 1 1]*max(5, max(icD(2,:))+2), cBi, 'FaceAlpha', 0.06, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        fill(ax, [xlimMC fliplr(xlimMC)], [0 0 1 1]*min(-5, min(icD(1,:))-2), cM, 'FaceAlpha', 0.06, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        errorbar(ax, S.eps, dMean, dMean - icD(1,:), icD(2,:) - dMean, 'kd', 'MarkerFaceColor', 'k', 'LineWidth', 1.3, 'CapSize', 5);
        yline(ax, 0, 'k-');
        text(ax, xlimMC(1)+1, max(5, max(icD(2,:))+2)*0.85, 'especular mejor', 'Color', cBi, 'FontWeight', 'bold');
        text(ax, xlimMC(1)+1, min(-5, min(icD(1,:))-2)*0.85, 'mono mejor', 'Color', cM, 'FontWeight', 'bold');
        xline(ax, eps_B, '--', 'Color', cB, 'LineWidth', 1.2);
        grid(ax, 'on'); xlim(ax, xlimMC); ylim(ax, [min(-5, min(icD(1,:))-2) max(5, max(icD(2,:))+2)]);
        xlabel(ax, '\epsilon [°]'); ylabel(ax, 'media especular - media mono [dB]');
        title(ax, sprintf('(b) Diferencia de potencia media (IC 95%%): mejor %d, peor %d, igual %d', nMas, nMenos, nIgual));

    end

    %% ==================== Guardar ====================
    if guardar
        tbl = table(S.eps(:), S.th_i(:), S.th_t(:), S.beta(:), S.Tx(1,:)', S.Tx(3,:)', S.Rx(1,:)', S.Rx(3,:)', ...
            S.R_T(:), S.R_RB(:), S.T_down(:), S.T_upB(:), S.Gamma2(:), ...
            10*log10(S.sigM(:)), 10*log10(S.sigB(:)), S.PrM_dBm(:), S.PrB_dBm(:), S.PrSup_dBm(:), S.PrRef_dBm(:), ...
            'VariableNames', {'eps_deg','theta_i_deg','theta_t_deg','beta_suelo_deg','x_Tx_m','z_Tx_m','x_Rx_m','z_Rx_m', ...
            'R_T_m','R_R_m','T_down','T_up','Gamma_TM2', ...
            'sigma_mono_dBsm','sigma_bi_dBsm','Pr_mono_dBm','Pr_bi_dBm','Pr_superficie_dBm','Pr_sigmaRef_dBm'});
        csvFile = fullfile(outDir, sprintf('verifica_ganancia_brewster_%s.csv', stamp));
        pngFile = fullfile(outDir, sprintf('verifica_ganancia_brewster_%s.png', stamp));
        writetable(tbl, csvFile);
        exportgraphics(fig, pngFile, 'Resolution', 150);
        fprintf('\nGuardado: %s\n          %s\n          %s\n', csvFile, pngFile, logFile);
        if patrones360
            pngMFile = fullfile(outDir, sprintf('rcs3d_mono_%s.png', stamp));
            pngBFile = fullfile(outDir, sprintf('rcs3d_biestatico_%s.png', stamp));
            csv2File = fullfile(outDir, sprintf('patrones_rcs_360_%s.csv', stamp));
            mat3File = fullfile(outDir, sprintf('rcs3d_%s.mat', stamp));
            exportgraphics(fig3M, pngMFile, 'Resolution', 150);
            exportgraphics(fig3B, pngBFile, 'Resolution', 150);
            writetable(table(psi(:), sMonoXZ, sBiXZ, sMonoAZ, sBiAZ, 'VariableNames', ...
                {'angulo_deg','mono_xz_dBsm','bi_xz_dBsm','mono_az_dBsm','bi_az_dBsm'}), csv2File);
            % Mismo formato que rcs_core (Nel x Naz): imagesc(az_deg, el_deg, rcs_mono_dbsm)
            az_deg = az3; el_deg = el3; rcs_mono_dbsm = sMono3; rcs_bi_dbsm = sBi3;
            tx_bi_az_el_deg = [0, 90 - thtB]; lambda_suelo_m = lambdaS;
            save(mat3File, 'az_deg', 'el_deg', 'rcs_mono_dbsm', 'rcs_bi_dbsm', 'tx_bi_az_el_deg', 'lambda_suelo_m');
            fprintf('          %s\n          %s\n          %s\n          %s\n', pngMFile, pngBFile, csv2File, mat3File);
        end
        if montecarlo
            png4File = fullfile(outDir, sprintf('montecarlo_%s.png', stamp));
            csv4File = fullfile(outDir, sprintf('montecarlo_resumen_%s.csv', stamp));
            mat4File = fullfile(outDir, sprintf('montecarlo_%s.mat', stamp));
            exportgraphics(fig4, png4File, 'Resolution', 150);
            writetable(table(S.eps(:), S.th_i(:), meanM(:), icM(1,:)', icM(2,:)', meanB(:), icB(1,:)', icB(2,:)', ...
                dMean(:), icD(1,:)', icD(2,:)', pM(2,:)', pB(2,:)', pGana(:), hitM(:), hitB(:), ...
                'VariableNames', {'eps_deg','theta_i_deg','Pr_mono_media_dBm','IC95_mono_inf','IC95_mono_sup', ...
                'Pr_esp_media_dBm','IC95_esp_inf','IC95_esp_sup','dif_media_dB','IC95_dif_inf','IC95_dif_sup', ...
                'Pr_mono_mediana_dBm','Pr_esp_mediana_dBm','P_esp_mayor','aciertos_mono','aciertos_esp'}), csv4File);
            % Datos crudos: sigma por repetición y posición (Pr = K .* sigma)
            eps_mc = S.eps; sigma_mono_m2 = sigMmc; sigma_esp_m2 = sigBmc; K_mono = S.KM; K_esp = S.KB;
            rot_xy_deg = angMC;                                   % [Nmc x Ne x 2] (NaN en modo 'so3')
            save(mat4File, 'eps_mc', 'sigma_mono_m2', 'sigma_esp_m2', 'K_mono', 'K_esp', 'rot_xy_deg', ...
                'modoRotacion', 'rotMax', 'semillaMC', 'Nmc');
            fprintf('          %s\n          %s\n          %s\n', png4File, csv4File, mat4File);
        end
    end

    % ---------- Resumen de este objeto ----------
    kOffS = find(S.eps < 90);
    res(iObj).nombre     = nombreObj;
    res(iObj).desc       = descObj;
    res(iObj).Ntri       = size(Fc, 1);
    res(iObj).sigNadir   = 10*log10(sigmaNadir);
    res(iObj).dSinRot    = 10*log10(sum(S.PrB(kOffS)) / sum(S.PrM(kOffS)));
    res(iObj).dSinRotEpB = S.PrB_dBm(kB5) - S.PrM_dBm(kB5);
    [res(iObj).dMC, res(iObj).icMC1, res(iObj).icMC2, res(iObj).pGana, res(iObj).nMas, ...
        res(iObj).nMenos, res(iObj).nIgual, res(iObj).hitMin] = deal(NaN);
    if montecarlo
        rB = zeros(Nboot, 1);
        for b = 1:Nboot                                   % IC 95% bootstrap de la diferencia global
            ib = randi(Nmc, Nmc, 1);
            rB(b) = 10*log10(sum(PrBmc(ib,kOff), 'all') / sum(PrMmc(ib,kOff), 'all'));
        end
        res(iObj).dMC    = 10*log10(sum(PrBmc(:,kOff), 'all') / sum(PrMmc(:,kOff), 'all'));
        res(iObj).icMC1  = pctl(rB, 2.5);  res(iObj).icMC2 = pctl(rB, 97.5);
        res(iObj).pGana  = mean(pGana(kOff));
        res(iObj).nMas   = nMas; res(iObj).nMenos = nMenos; res(iObj).nIgual = nIgual;
        res(iObj).hitMin = min([hitM(kOff) hitB(kOff)]);
    end
    if cerrarFiguras; close all; end
end   % bucle de objetos

%% ==================== Resumen comparativo entre objetos ====================
fprintf('\n\n==================== RESUMEN: MONOESTÁTICO vs ESPECULAR POR OBJETO ====================\n');
fprintf('  B-M = 10log10(energía especular / energía mono) sumando las posiciones fuera del nadir.\n');
fprintf('  "Sin rotar": objeto en su orientación nominal. "Monte Carlo": %d orientaciones por posición (''%s'').\n', ...
    Nmc, modoRotacion);
fprintf('%-26s %6s %10s | %12s %13s | %22s %7s %14s %8s\n', 'Objeto', 'Ntri', 'sig nadir', ...
    'Sin rotar', 'Sin rot. eB', 'Monte Carlo B-M (IC95)', 'P(B>M)', 'mejor/peor/=', 'aciertos');
for i = 1:numel(res)
    fprintf('%-26s %6d %8.2f dBsm | %+9.2f dB %+10.2f dB | %+7.2f [%+5.2f,%+5.2f] dB %7.2f %6d/%d/%d %8d\n', ...
        res(i).desc, res(i).Ntri, res(i).sigNadir, res(i).dSinRot, res(i).dSinRotEpB, ...
        res(i).dMC, res(i).icMC1, res(i).icMC2, res(i).pGana, res(i).nMas, res(i).nMenos, res(i).nIgual, res(i).hitMin);
end
fprintf('  (Sin rot. eB: diferencia en la posición del barrido más cercana a Brewster; aciertos: mínimo por posición)\n');

figR = figure('Name', 'Resumen por objeto', 'Position', [100 100 1400 520]);
tlR  = tiledlayout(figR, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tlR, sprintf('Monoestático vs biestático especular por objeto | f_c = %g GHz, TM\\rightarrowTM, arco r = %g m, z = %g m', ...
    fc/1e9, r_radar, z_obj));
xo = 1:numel(res); etq = {res.nombre};
ax = nexttile(tlR); hold(ax, 'on');
bar(ax, xo, [res.dSinRot], 0.55, 'FaceColor', cBi, 'EdgeColor', 'none');
yline(ax, 0, 'k-');
text(ax, xo, [res.dSinRot], compose('%+.1f dB', [res.dSinRot]), 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom');
set(ax, 'XTick', xo, 'XTickLabel', etq); grid(ax, 'on'); xlim(ax, [0.4 numel(res)+0.6]);
ylabel(ax, 'energía especular / mono [dB]');
title(ax, '(a) Objeto en orientación nominal (sin rotar)');
ax = nexttile(tlR); hold(ax, 'on');
dmc = [res.dMC];
errorbar(ax, xo, dmc, dmc - [res.icMC1], [res.icMC2] - dmc, 'kd', 'MarkerFaceColor', 'k', 'LineWidth', 1.3, 'CapSize', 8);
yline(ax, 0, 'k-');
text(ax, xo + 0.12, dmc, compose('%+.2f dB', dmc), 'VerticalAlignment', 'middle');
set(ax, 'XTick', xo, 'XTickLabel', etq); grid(ax, 'on'); xlim(ax, [0.4 numel(res)+0.6]);
ylabel(ax, 'energía especular / mono [dB]');
title(ax, sprintf('(b) Monte Carlo: orientación aleatoria (%d por posición), IC 95%%', Nmc));

if guardar
    resFile = fullfile(outRoot, sprintf('resumen_objetos_%s.csv', stamp));
    resPng  = fullfile(outRoot, sprintf('resumen_objetos_%s.png', stamp));
    writetable(struct2table(res), resFile);
    exportgraphics(figR, resPng, 'Resolution', 150);
    fprintf('\nResumen guardado: %s\n                  %s\n', resFile, resPng);
    diary off;
end

%% =====================================================================
%% Funciones locales
%% =====================================================================
function S = simularBarrido(eps_deg, r, target, n1, n2, V, F, lambda0, lambdaS, Pt, Gt, Gr, sigmaRef, alphaS, incluirPerdas)
% Tx en el arco de radio r (centro (0,0,0), lado +x) y Rx especular (lado -x). Traza los
% rayos con Fermat/Snell hasta el centro del cubo y evalúa la ecuación radar para el
% caso monoestático (Rx = Tx) y el biestático especular.
    N = numel(eps_deg);
    S.eps = eps_deg;
    S.Tx  = [ r*cosd(eps_deg); zeros(1,N); r*sind(eps_deg)];
    S.Rx  = [-r*cosd(eps_deg); zeros(1,N); r*sind(eps_deg)];
    S.rTx = sqrt(sum(S.Tx.^2, 1));
    S.rRx = sqrt(sum(S.Rx.^2, 1));
    [S.PT, S.PR] = deal(zeros(3,N));
    [S.iHat, S.sHat] = deal(zeros(N,3));
    [S.th_i, S.th_t, S.th_iR, S.beta, S.R_T, S.R_RB, S.R_sueT, S.R_sueR, S.T_down, S.T_upM, S.T_upB, ...
        S.Gamma2, S.sigM, S.sigB] = deal(zeros(1,N));
    for k = 1:N
        Tx = S.Tx(:,k); Rx = S.Rx(:,k);
        % --- Tramo Tx -> cubo (mismas funciones que el resto del pipeline) ---
        PT  = calculateRefractionPointFermat(Tx, target, n1, n2);
        thi = calculateIncidenceAngle(Tx, PT, target);
        [~, Gam, ~, ~, Td] = calculateTMcoef(thi, n1, n2);
        tht = asin(n1/n2*sin(thi));
        [~,~,~,~, TuM] = calculateTMcoef(tht, n2, n1);          % subida monoestática (mismo rayo)
        % --- Tramo cubo -> Rx especular ---
        PR   = calculateRefractionPointFermat(target, Rx, n2, n1);
        thtR = calculateIncidenceAngle(target, PR, Rx);          % ángulo en el suelo
        [~,~,~,~, TuB] = calculateTMcoef(thtR, n2, n1);
        S.PT(:,k) = PT; S.PR(:,k) = PR;
        S.th_i(k) = rad2deg(thi);
        S.th_t(k) = rad2deg(tht);
        S.th_iR(k) = rad2deg(asin(n2/n1*sin(thtR)));
        S.R_sueT(k) = norm(PT - target);
        S.R_sueR(k) = norm(PR - target);
        S.R_T(k)  = norm(Tx - PT) + S.R_sueT(k);
        S.R_RB(k) = norm(PR - Rx) + S.R_sueR(k);
        S.T_down(k) = Td; S.T_upM(k) = TuM; S.T_upB(k) = TuB;
        S.Gamma2(k) = Gam^2;
        % Direcciones de observación desde el cubo hacia cada punto de refracción (en el suelo)
        iHat = (PT - target)' / S.R_sueT(k);
        sHat = (PR - target)' / S.R_sueR(k);
        S.iHat(k,:) = iHat; S.sHat(k,:) = sHat;
        S.beta(k) = acosd(max(min(iHat*sHat', 1), -1));
        S.sigM(k) = rcsPO(V, F, iHat, iHat, lambdaS);
        S.sigB(k) = rcsPO(V, F, iHat, sHat, lambdaS);
    end
    S.T2_M = S.T_down .* S.T_upM;
    S.T2_B = S.T_down .* S.T_upB;
    [LM, LB] = deal(ones(1,N));
    if incluirPerdas
        LM = exp(-2*alphaS*S.R_sueT);
        LB = exp(-alphaS*(S.R_sueT + S.R_sueR));
    end
    KM = Pt * Gt * Gr * lambda0^2 .* S.T2_M .* LM ./ ((4*pi)^3 * S.R_T.^4);
    KB = Pt * Gt * Gr * lambda0^2 .* S.T2_B .* LB ./ ((4*pi)^3 * S.R_T.^2 .* S.R_RB.^2);
    S.KM = KM; S.KB = KB;                           % Pr = K * sigma (se reutiliza en el Monte Carlo)
    S.PrM   = KM .* S.sigM;
    S.PrB   = KB .* S.sigB;
    S.PrRef = KM .* sigmaRef;
    % Reflexión especular de la superficie (imagen del Tx), camino Tx -> (0,0,0) -> Rx
    S.PrSup = Pt * Gt * Gr * lambda0^2 .* S.Gamma2 ./ ((4*pi)^2 * (S.rTx + S.rRx).^2);
    S.PrM_dBm   = 10*log10(S.PrM*1e3);
    S.PrB_dBm   = 10*log10(S.PrB*1e3);
    S.PrRef_dBm = 10*log10(S.PrRef*1e3);
    S.PrSup_dBm = 10*log10(max(S.PrSup, 1e-30)*1e3);
end

function sigma = rcsPO(V, F, iHat, sHat, lambda, eI, eR)
% RCS biestática PO (TM->TM) de un PEC cerrado. Generaliza rcs_core/rcs/po.py:
%   sqrt(sigma) = (2*sqrt(pi)/lambda) * sum_{caras con n.i > 0} [e_r.(n x (k_i x e_i))] * I_T(k*(i+s))
% iHat, sHat: [Nd x 3] direcciones unitarias objeto->Tx y objeto->Rx. Con sHat = iHat
% el factor de polarización vale n.i y la fase es 2k*i, idéntico a po.py.
% TM: E en el plano que contiene z y la dirección de propagación (V de rcs_core);
% se construye como y x k para que sea continuo en nadir dentro del plano xz.
% eI, eR (opcionales, [Nd x 3]): polarizaciones explícitas; permiten evaluar un cubo rotado
% R sin rotar la malla, pasando R'*i, R'*s, R'*e_i, R'*e_r (todo es invariante a rotaciones).
    k  = 2*pi/lambda;
    p0 = V(F(:,1),:); p1 = V(F(:,2),:); p2 = V(F(:,3),:);
    cr = cross(p1 - p0, p2 - p0, 2);
    area2 = sqrt(sum(cr.^2, 2));
    n = cr ./ area2;                                % [M x 3]
    Nd = size(iHat, 1);
    kI = -iHat;                                     % propagación incidente
    if nargin < 6
        y  = repmat([0 1 0], Nd, 1);
        eI = normalizaFilas(cross(y, kI, 2));
        eR = normalizaFilas(cross(y, sHat, 2));
    end
    hI = cross(kI, eI, 2);                          % dirección de H incidente
    % factor(d,m) = e_r . (n_m x h_i) = n_m . (h_i x e_r)
    q  = cross(hI, eR, 2);                          % [Nd x 3]
    fac = q * n';                                   % [Nd x M]
    lit = (iHat * n') > 0;
    w = k*(iHat + sHat);                            % [Nd x 3]
    a = w * p0'; b = w * p1'; cc = w * p2';
    I = area2' .* exp(1i*a) .* triIntegral(b - a, cc - a);
    s = (2*sqrt(pi)/lambda) * sum(lit .* fac .* I, 2);
    sigma = abs(s).^2;
end

function u = normalizaFilas(u)
    u = u ./ sqrt(sum(u.^2, 2));
end

function out = triIntegral(x, y)
% G(x,y) = int_0^1 int_0^{1-u} exp(j(xu + yv)) dv du  (forma cerrada estable, como en po.py)
    swap = abs(x) > abs(y);
    xt = x; x(swap) = y(swap); y(swap) = xt(swap);  % ahora |y| >= |x|
    out = complex(zeros(size(x)));
    sm = abs(y) < 1e-3;
    xs = x(sm); ys = y(sm);
    out(sm) = 0.5 + 1i*(xs + ys)/6 - (xs.^2 + xs.*ys + ys.^2)/24 ...
              - 1i*(xs.^3 + xs.^2.*ys + xs.*ys.^2 + ys.^3)/120;
    xl = x(~sm); yl = y(~sm);
    out(~sm) = (exp(1i*yl) .* phi1(xl - yl) - phi1(xl)) ./ (1i*yl);
end

function out = phi1(t)
% (exp(jt) - 1)/(jt), estable cerca de 0
    out = complex(ones(size(t)));
    sm = abs(t) < 1e-4;
    out(sm) = 1 + 1i*t(sm)/2 - t(sm).^2/6;
    out(~sm) = (exp(1i*t(~sm)) - 1) ./ (1i*t(~sm));
end

function [V, F] = mallaCubo(L)
% Cubo centrado en el origen, 12 triángulos con normales salientes (PO integra la fase
% exactamente por triángulo, no hace falta refinar la malla).
    V = L/2 * [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1];
    F = [1 3 2; 1 4 3;    % z = -L/2
         5 6 7; 5 7 8;    % z = +L/2
         1 2 6; 1 6 5;    % y = -L/2
         4 8 7; 4 7 3;    % y = +L/2
         1 5 8; 1 8 4;    % x = -L/2
         2 3 7; 2 7 6];   % x = +L/2
end

function [V, F, desc] = mallaObjeto(nombre, dim, lam)
% Malla triangular cerrada (normales salientes) de los sólidos de rcs_core/rcs/geometry.py,
% centrados en el origen. Las caras planas no se subdividen (PO integra la fase exacta por
% triángulo); las curvas se facetan con flecha < lam/100, como simular.py con PO.
    nSeg = @(r) max(12, ceil(pi / acos(max(1 - lam/(100*r), -1))));   % segmentos por vuelta
    switch nombre
        case {'cubo', 'caja'}
            if isscalar(dim); dim = dim*[1 1 1]; end
            [V, F] = mallaCubo(1);
            V = V .* dim(:)';
            if strcmp(nombre, 'cubo')
                desc = sprintf('cubo %g m', dim(1));
            else
                desc = sprintf('caja %gx%gx%g m', dim);
            end
        case 'cilindro'
            r = dim(1); L = dim(2); N = nSeg(r);
            ph = 2*pi*(0:N-1)'/N;
            V = [r*cos(ph) r*sin(ph) -L/2*ones(N,1); r*cos(ph) r*sin(ph) L/2*ones(N,1); 0 0 -L/2; 0 0 L/2];
            i1 = (1:N)'; i2 = mod(i1, N) + 1;
            F = [i1 i2 N+i2; i1 N+i2 N+i1; ...                      % lateral
                 (2*N+1)*ones(N,1) i2 i1; (2*N+2)*ones(N,1) N+i1 N+i2];   % tapas
            desc = sprintf('cilindro r=%g m L=%g m', r, L);
        case 'cono'
            r = dim(1); L = dim(2); N = nSeg(r);
            ph = 2*pi*(0:N-1)'/N;
            V = [r*cos(ph) r*sin(ph) -L/2*ones(N,1); 0 0 L/2; 0 0 -L/2];
            i1 = (1:N)'; i2 = mod(i1, N) + 1;
            F = [i1 i2 (N+1)*ones(N,1); (N+2)*ones(N,1) i2 i1];
            desc = sprintf('cono r=%g m L=%g m', r, L);
        case 'esfera'
            r = dim(1); N = nSeg(r); M = ceil(N/2);
            ph = 2*pi*(0:N-1)/N; th = pi*(1:M-1)'/M;            % anillos sin los polos
            X = r*sin(th)*cos(ph); Y = r*sin(th)*sin(ph); Z = r*cos(th)*ones(1,N);
            V = [0 0 r; X(:) Y(:) Z(:); 0 0 -r];                   % índice (anillo a, meridiano b)
            idx = @(a, b) 1 + (mod(b-1, N))*(M-1) + a;
            F = zeros(0, 3);
            for b = 1:N
                F = [F; 1 idx(1,b) idx(1,b+1)];                    %#ok<AGROW> casquete norte
                for a = 1:M-2
                    F = [F; idx(a,b) idx(a+1,b) idx(a+1,b+1); idx(a,b) idx(a+1,b+1) idx(a,b+1)]; %#ok<AGROW>
                end
                F = [F; size(V,1) idx(M-1,b+1) idx(M-1,b)];        %#ok<AGROW> casquete sur
            end
            desc = sprintf('esfera r=%g m', r);
        otherwise
            error('verifica_ganancia_brewster:objeto', 'Objeto desconocido: %s', nombre);
    end
    F = orientarFuera(V, F);
end

function F = orientarFuera(V, F)
% Normales salientes (válido para cuerpos estrellados respecto al centroide), como _orient_outward.
    c  = mean(V, 1);
    p0 = V(F(:,1),:); p1 = V(F(:,2),:); p2 = V(F(:,3),:);
    n  = cross(p1 - p0, p2 - p0, 2);
    flip = sum(n .* ((p0 + p1 + p2)/3 - c), 2) < 0;
    F(flip,:) = F(flip, [1 3 2]);
end

function validarPOcontraRcsCore(dirRes, c)
% Compara el PO de MATLAB (con sHat = iHat) con los .mat generados por tools/rcs_core/simular.py
% (cualquier sólido que sepa mallar mallaObjeto: cubo, caja, cilindro, esfera, cono).
    fprintf('-- Validación del PO de MATLAB contra tools/rcs_core/resultados (monoestático) --\n');
    archivos = dir(fullfile(dirRes, '*.mat'));
    if isempty(archivos)
        fprintf('  No hay .mat en %s; se omite la validación.\n\n', dirRes);
        return;
    end
    for i = 1:numel(archivos)
        r = load(fullfile(dirRes, archivos(i).name));
        if ~isfield(r, 'forma') || ~isfield(r, 'rcs_dbsm'); continue; end
        lam = c / r.freq_Hz;
        [V, F] = mallaObjeto(strtrim(r.forma), r.dim_m, lam);
        % Corte az = 0° (mismo plano de incidencia que la simulación), el de 0 a 90°
        el = r.el_deg(r.el_deg >= 0);
        [~, i0] = min(abs(r.az_deg));
        ref = r.rcs_dbsm(r.el_deg >= 0, i0)';
        rh = [cosd(el(:)) zeros(numel(el),1) sind(el(:))];
        po = 10*log10(max(rcsPO(V, F, rh, rh, lam), 1e-30))';
        m = ref > max(ref) - 20;                    % lóbulos relevantes (20 dB bajo el máximo)
        [~, iN] = max(el);
        fprintf('  %-38s (%s, %.2f lambda): nadir %.2f vs %.2f dBsm | error medio |dB| (lóbulos) = %.2f dB\n', ...
            archivos(i).name, strtrim(r.metodo), max(max(V) - min(V))/lam, po(iN), ref(iN), mean(abs(po(m) - ref(m))));
    end
    fprintf('  (con PO el error debe ser ~0 dB; frente a MoM la diferencia refleja las limitaciones de PO)\n\n');
end

function validarPObiestatico(L, lam)
% Placa especular (cara superior del cubo): sigma = 4*pi*A^2*cos^2(theta)/lambda^2 (TM).
    fprintf('-- Validación del PO biestático: especular de la cara superior --\n');
    [V, F] = mallaCubo(L);
    for th = [0 15 30 45]
        iHat = [ sind(th) 0 cosd(th)];
        sHat = [-sind(th) 0 cosd(th)];
        sPO  = rcsPO(V, F, iHat, sHat, lam);
        sTeo = 4*pi*L^4*cosd(th)^2/lam^2;
        fprintf('  theta = %4.1f°: PO %.3f dBsm | teórico %.3f dBsm\n', th, 10*log10(sPO), 10*log10(sTeo));
    end
    fprintf('\n');
end

function diagramaRCS3D(fig, AZ, EL, SdB, rMin, sMax, Vc, Fc, iTx, sEsp, thCr, titulo)
% Diagrama 3D de RCS como en rcs_core/leer_rcs.m: radio = dBsm - rMin (recortado a 0),
% color = dBsm. Añade el objeto (a escala del diagrama), la dirección del Tx (morado), la
% especular (naranja, si se da) y el borde del cono accesible desde el aire (punteado).
    ax = axes(fig); hold(ax, 'on');
    R = max(SdB - rMin, 0);
    surf(ax, R.*cosd(EL).*cosd(AZ), R.*cosd(EL).*sind(AZ), R.*sind(EL), SdB, 'EdgeColor', 'none');
    colormap(ax, turbo); caxis(ax, [rMin sMax]);
    cb = colorbar(ax); cb.Label.String = '\sigma [dBsm]';
    Rmax = sMax - rMin;
    L = 1.25*Rmax;
    esc = 0.04*Rmax / max(abs(Vc(:)));
    patch(ax, 'Vertices', Vc*esc, 'Faces', Fc, ...
        'FaceColor', [0.55 0.55 0.6], 'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.35);
    quiver3(ax, 0, 0, 0, L*iTx(1), L*iTx(2), L*iTx(3), 0, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2.5, 'MaxHeadSize', 0.3);
    text(ax, 1.08*L*iTx(1), 1.08*L*iTx(2), 1.08*L*iTx(3), 'Tx (Brewster)', 'Color', [0.4940 0.1840 0.5560], 'FontWeight', 'bold');
    if ~isempty(sEsp)
        quiver3(ax, 0, 0, 0, L*sEsp(1), L*sEsp(2), L*sEsp(3), 0, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2.5, 'MaxHeadSize', 0.3);
        text(ax, 1.08*L*sEsp(1), 1.08*L*sEsp(2), 1.08*L*sEsp(3), 'Rx especular', 'Color', [0.8500 0.3250 0.0980], 'FontWeight', 'bold');
    end
    ph = linspace(0, 360, 181);
    plot3(ax, L*sind(thCr)*cosd(ph), L*sind(thCr)*sind(ph), L*cosd(thCr)*ones(size(ph)), 'k:', 'LineWidth', 1.2);
    for p = 0:90:270
        plot3(ax, [0 L*sind(thCr)*cosd(p)], [0 L*sind(thCr)*sind(p)], [0 L*cosd(thCr)], 'k:');
    end
    axis(ax, 'equal'); grid(ax, 'on'); view(ax, -35, 22);
    camlight(ax, 'headlight'); lighting(ax, 'gouraud'); material(ax, 'dull');
    % El radio es dB sobre el piso, no metros: se ocultan los números de los ejes
    set(ax, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', []);
    xlabel(ax, 'x'); ylabel(ax, 'y'); zlabel(ax, 'z');
    title(ax, {titulo, sprintf('radio = \\sigma - (%.1f dBsm) | punteado: cono accesible desde el aire (%.0f°)', rMin, thCr)});
end

function P = pctl(X, p)
% Percentiles por columnas con interpolación lineal (sin Statistics Toolbox). X: [n x m].
    if isrow(X); X = X(:); end
    Xs = sort(X, 1);
    n  = size(Xs, 1);
    P  = zeros(numel(p), size(Xs, 2));
    for i = 1:numel(p)
        h  = 1 + (n - 1)*p(i)/100;
        lo = floor(h); hi = min(lo + 1, n);
        P(i,:) = Xs(lo,:) + (h - lo)*(Xs(hi,:) - Xs(lo,:));
    end
end

function [Rall, ang] = orientacionesAleatorias(N, modo, rotMax)
% N matrices de rotación aleatorias, por filas [R11 R12 R13 R21 R22 R23 R31 R32 R33].
%   'xy' : R = Ry(ry)*Rx(rx), rx y ry uniformes en [-rotMax, rotMax]  (ang = [rx ry])
%   'so3': orientación uniforme en 3D (cuaternión gaussiano normalizado) (ang = NaN)
    switch modo
        case 'xy'
            ang = rotMax*(2*rand(N, 2) - 1);
            cx = cosd(ang(:,1)); sx = sind(ang(:,1));
            cy = cosd(ang(:,2)); sy = sind(ang(:,2));
            Rall = [cy, sy.*sx, sy.*cx, zeros(N,1), cx, -sx, -sy, cy.*sx, cy.*cx];
        case 'so3'
            ang = nan(N, 2);
            q = randn(N, 4); q = q ./ sqrt(sum(q.^2, 2));
            w = q(:,1); x = q(:,2); y = q(:,3); z = q(:,4);
            Rall = [1-2*(y.^2+z.^2), 2*(x.*y-z.*w), 2*(x.*z+y.*w), ...
                    2*(x.*y+z.*w), 1-2*(x.^2+z.^2), 2*(y.*z-x.*w), ...
                    2*(x.*z-y.*w), 2*(y.*z+x.*w), 1-2*(x.^2+y.^2)];
        otherwise
            error('verifica_ganancia_brewster:modoRotacion', 'modoRotacion desconocido: %s', modo);
    end
end

function R = rotX(a)
    R = [1 0 0; 0 cosd(a) -sind(a); 0 sind(a) cosd(a)];
end

function R = rotY(a)
    R = [cosd(a) 0 sind(a); 0 1 0; -sind(a) 0 cosd(a)];
end

function a = wrapTo180loc(a)
    a = mod(a + 180, 360) - 180;
end

function s = ternario(cond, a, b)
    if cond; s = a; else; s = b; end
end
