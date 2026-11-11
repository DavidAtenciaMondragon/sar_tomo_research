% FIG_INFORME  Figuras adicionales para informe_sincronismo.tex.
%
% Requiere haber ejecutado run_all (lee los CSV de resultados/). Genera:
%   inf_escenario     vuelo circular 3D, planta con iso-rangos bistaticos,
%                     soporte espectral e IRF ideal
%   inf_senal         patrones de escucha, ventana grabada, dato comprimido
%                     y llegadas del directo antes/despues de asociar
%   inf_reloj         residuo de fase de un modelo de reloj por oscilador
%   inf_cortes        cortes de la IRF del blanco central (E0..E3, ideal)
%   inf_conclusiones  resumen grafico de los experimentos 3-6
% (~5 min)

clear; inicializa;
P = params_sinc();
c = P.c; fc = P.fc; lam = c/fc;
col = lines(7);
tol_phi = 0.47;                                   % rad, 1 dB (exp1)

% ======================= 1. Escenario ===================================
% Adquisicion continua (TCXO, SNR_d 30 dB): se usa tambien en la seccion 2
tic; [D, T, TX] = simula_adquisicion(P, 21);
S = sincroniza(D, TX, P);
fprintf('inf_escenario/senal: simulacion %.0f s, J = %d\n', toc, numel(D.m));
tg = P.targets; s = P.escena_semi;
Tv = 2*pi*P.radio/P.v;                            % periodo de la vuelta
th = linspace(0, 2*pi, 361).';
circ = [P.radio*cos(th), P.radio*sin(th), P.alt*ones(size(th))];
t_snap = P.rx.t_start + [0 0.125 0.25]*Tv;        % tres instantes
[tx_s, vtx_s] = trayectoria(P, 'tx', t_snap);
rx_s = trayectoria(P, 'rx', t_snap);
Rd = norm(tx_s(1,:) - rx_s(1,:));

f = figure('Position', [50 50 1150 900]);
subplot(2,2,1); hold on; grid on; box on;
plot3(circ(:,1), circ(:,2), circ(:,3), '-', 'Color', [0.5 0.5 0.5]);
plot3(circ(:,1), circ(:,2), 0*circ(:,3), ':', 'Color', [0.7 0.7 0.7]);
plot3([-s s s -s -s], [-s -s s s -s], zeros(1,5), 'k-', 'LineWidth', 0.5);
plot3(tg(:,1), tg(:,2), tg(:,3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
for i = numel(t_snap):-1:1                        % t0 encima
    a = 1 - 0.6*(i > 1);                          % instantes posteriores atenuados
    ct = a*col(1,:) + (1-a); cr = a*col(2,:) + (1-a);
    plot3([tx_s(i,1) rx_s(i,1)], [tx_s(i,2) rx_s(i,2)], [tx_s(i,3) rx_s(i,3)], '--', 'Color', cr, 'LineWidth', 1.2);
    plot3(tx_s(i,1), tx_s(i,2), tx_s(i,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', ct, 'Color', ct);
    plot3(rx_s(i,1), rx_s(i,2), rx_s(i,3), 'v', 'MarkerSize', 8, 'MarkerFaceColor', cr, 'Color', cr);
end
plot3([tx_s(1,1) 0 rx_s(1,1)], [tx_s(1,2) 0 rx_s(1,2)], [tx_s(1,3) 0 rx_s(1,3)], '-', 'Color', col(5,:));
quiver3(tx_s(1,1), tx_s(1,2), tx_s(1,3), 2*vtx_s(1,1), 2*vtx_s(1,2), 0, 0, 'Color', col(1,:), 'LineWidth', 1.5);
text(tx_s(1,1) - 5, tx_s(1,2), tx_s(1,3) - 14, 'Radar 1 (Tx/Rx)', 'Color', col(1,:), 'HorizontalAlignment', 'right');
text(rx_s(1,1) + 5, rx_s(1,2), rx_s(1,3) + 12, 'Radar 2 (Rx pasivo)', 'Color', col(2,:));
text(mean([tx_s(1,1) rx_s(1,1)]) + 10, mean([tx_s(1,2) rx_s(1,2)]) - 10, P.alt + 4, ...
    sprintf('directo %.1f m', Rd), 'Color', col(2,:), 'FontSize', 8);
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal;
view(-35, 25); xlim([-60 60]); ylim([-60 60]); zlim([0 110]);
title(sprintf('(a) Vuelo circular: R = %g m, h = %g m, v = %g m/s', P.radio, P.alt, P.v));

subplot(2,2,2); hold on; box on;
u = -60:0.5:60; [X, Y] = meshgrid(u, u);
tx = tx_s(1,:); rx = rx_s(1,:);
Rb = hypot(hypot(X - tx(1), Y - tx(2)), tx(3)) + hypot(hypot(X - rx(1), Y - rx(2)), rx(3));
[C, h] = contour(X, Y, Rb - Rd, 0:10:250, 'Color', [0.65 0.65 0.65]);
clabel(C, h, 100:20:250, 'FontSize', 7, 'Color', [0.45 0.45 0.45]);
plot(circ(:,1), circ(:,2), '-', 'Color', [0.5 0.5 0.5]);
plot([-s s s -s -s], [-s -s s s -s], 'k-', 'LineWidth', 0.5);
plot(tg(:,1), tg(:,2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
for i = numel(t_snap):-1:1
    a = 1 - 0.6*(i > 1); ct = a*col(1,:) + (1-a); cr = a*col(2,:) + (1-a);
    plot([tx_s(i,1) rx_s(i,1)], [tx_s(i,2) rx_s(i,2)], '--', 'Color', cr);
    plot(tx_s(i,1), tx_s(i,2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', ct, 'Color', ct);
    plot(rx_s(i,1), rx_s(i,2), 'v', 'MarkerSize', 8, 'MarkerFaceColor', cr, 'Color', cr);
end
plot([0 tx(1)], [0 tx(2)], 'k:'); plot([0 rx(1)], [0 rx(2)], 'k:');
text(6, -14, '90^\circ', 'FontSize', 9);
gx = (0 - tx(1))/norm(-tx) + (0 - rx(1))/norm(-rx); gy = (0 - tx(2))/norm(-tx) + (0 - rx(2))/norm(-rx);
quiver(0, 0, 12*gx/hypot(gx, gy), 12*gy/hypot(gx, gy), 0, 'Color', col(5,:), 'LineWidth', 1.5, 'MaxHeadSize', 0.8);
text(-34, 16, sprintf('\\nabla R_b (|\\nabla R_b| = %.2f)', hypot(gx, gy)), 'Color', col(5,:), 'FontSize', 8);
axis equal; xlim([-62 62]); ylim([-62 62]); grid on;
xlabel('x (m)'); ylabel('y (m)');
title('(b) Planta: R_b - R_d (m) en t_0; instantes t_0 + T/8, T/4 atenuados');

% (c) soporte espectral sobre el suelo del blanco central (toda la vuelta)
subplot(2,2,3); hold on; box on; grid on;
nn = T.n(1:20:end);
txv = TX.pos(nn,:); rxv = TX.rx(nn,:); p0 = tg(1,:);
gR = (p0 - txv)./vecnorm(p0 - txv, 2, 2) + (p0 - rxv)./vecnorm(p0 - rxv, 2, 2);
for ff = [fc - P.B/2, fc, fc + P.B/2]
    kk = 2*pi*ff/c*gR(:,1:2);
    plot(kk(:,1), kk(:,2), '.', 'MarkerSize', 4, 'Color', col(1,:)*(ff == fc) + [0.6 0.6 0.6]*(ff ~= fc));
end
kr = 2*pi*fc/c*mean(vecnorm(gR(:,1:2), 2, 2));
axis equal; xlim(8*[-1 1]); ylim(8*[-1 1]);
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');
title(sprintf('(c) Soporte espectral (blanco central): anillo k_r = %.2f rad/m', kr));
text(-7.5, 7, sprintf('f_c \\pm B/2 (gris); resolucion ~ 2.405/k_r = %.2f m (nulo)', 2.405/kr), 'FontSize', 8);

% (d) IRF ideal del blanco central
subplot(2,2,4);
uu = -3:0.05:3; [Xp, Yp] = meshgrid(uu, uu);
Ip = forma_imagenes(D, S, T, TX, P, [Xp(:), Yp(:), zeros(numel(Xp),1)]);
A = abs(reshape(Ip(:,5), size(Xp))); A = A/max(A(:));
imagesc(uu, uu, 20*log10(A), [-35 0]); axis xy image; colormap(gca, parula);
cb = colorbar; cb.Label.String = 'dB';
xlabel('x (m)'); ylabel('y (m)');
title('(d) Respuesta al impulso BiSAR (reloj ideal, vuelta completa)');
guarda_fig(f, DIR_FIG, 'inf_escenario');

% ======================= 2. Senal recibida ==============================

f = figure('Position', [50 50 1150 700]);
% (a) patrones de escucha
subplot(2,2,1);
rs = RandStream('mt19937ar', 'Seed', 5);
pats = {'continuo', [], 'continuo'; 'periodico', 8, 'periodico 1/8'; ...
        'aleatorio', 0.5, 'aleatorio 50 %'; 'aleatorio', 0.05, 'aleatorio 5 %'; ...
        'rafagas', [20 180], 'rafagas 20/180'};
Mv = 400; Z = zeros(size(pats,1), Mv);
for i = 1:size(pats,1)
    o = patron_escucha(P.rx.M_slots, pats{i,1}, pats{i,2}, rs);
    Z(i,:) = o(1:Mv);
end
imagesc((0:Mv-1)*P.PRT0, 1:size(Z,1), Z); colormap(gca, [1 1 1; col(1,:)]);
set(gca, 'YTick', 1:size(Z,1), 'YTickLabel', pats(:,3), 'TickLength', [0 0]);
xlabel('tiempo de escucha (s)  [cada columna = una PRT]');
title(sprintf('(a) Patrones de escucha discontinua (1 s de %.0f s)', P.rx.M_slots*P.PRT0));

% (b) una ventana grabada: filtro adaptado (directo + 5 ecos)
subplot(2,2,2);
Ts = 1/P.fs; Np = round(P.Tp*P.fs); K = P.B/P.Tp;
pref = chirp_bb((0:Np-1).'*Ts, P.Tp, K);
n = T.n(round(end/2)); txp = TX.pos(n,:); rx = TX.rx(n,:);
tau_d = norm(txp - rx)/c;
tau_e = (vecnorm(tg - txp, 2, 2) + vecnorm(tg - rx, 2, 2))/c;
Ad = sqrt(10^(P.SNR_d_dB/10)/Np); Ae = sqrt(10^(P.SNR_e_dB/10)/Np);
rs = RandStream('mt19937ar', 'Seed', 3);
t0 = -0.5e-6; Nw = round(2e-6/Ts);
y = (randn(rs, Nw, 1) + 1j*randn(rs, Nw, 1))/sqrt(2);
dl = [0; tau_e - tau_d]; am = [Ad; Ae*ones(size(tau_e))].*exp(-2j*pi*fc*dl);
for p = 1:numel(dl)
    ii = (1:Nw).'; tl = t0 + (ii-1)*Ts - dl(p);
    on = tl >= 0 & tl < P.Tp;
    y(on) = y(on) + am(p)*chirp_bb(tl(on), P.Tp, K);
end
Nf = 4096; r = ifft(fft(y, Nf).*conj(fft(pref, Nf)));
r = r(1:Nw)/sqrt(sum(abs(pref).^2));
tt = (t0 + (0:Nw-1)*Ts)*1e9;
plot(tt, 20*log10(abs(r)), 'Color', col(1,:)); hold on; grid on;
yline(P.det_umbral_dB, 'r:', 'umbral de deteccion');
plot(0, 10*log10(10^(P.SNR_d_dB/10)*sum(abs(pref).^2)/Np), 'v', 'Color', col(2,:), 'MarkerFaceColor', col(2,:));
text(15, P.SNR_d_dB + 3, 'directo: referencia de tiempo y fase', 'Color', col(2,:));
plot(dl(2:end)*1e9, P.SNR_e_dB*ones(numel(tau_e),1) + 3, 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
text(680, P.SNR_e_dB + 4, 'ecos de los 5 blancos');
xlim([t0*1e9, 1300]); ylim([-10 40]);
xlabel('tiempo local relativo a la llegada del directo (ns)'); ylabel('|filtro adaptado| sobre ruido (dB)');
title('(b) Una ventana grabada por el Rx2');

% (c) dato comprimido tras cancelar el directo (todas las ventanas)
subplot(2,2,3);
L = size(D.seg, 2);
rexc = ((0:L-1)*D.Ts_up - P.seg_pre*Ts)*c;           % exceso de camino respecto al directo (m)
tn = TX.t(T.n);
A = 20*log10(abs(double(D.seg)) + 1e-6); A = A - max(A(:));
imagesc(rexc, tn, A, [-40 0]); colormap(gca, flipud(gray)); axis xy;
xlim([-10 250]); xlabel('c (t - a_{directo}) (m)  =  R_b - R_d'); ylabel('t (s, reloj Tx)');
title('(c) Dato comprimido alineado con el directo (cancelado)');
cb = colorbar; cb.Label.String = 'dB';

% (d) llegadas del directo: lo que ve el Rx2 vs tras asociar
subplot(2,2,4);
crudo = D.a_hat - D.a_hat(1) - D.m*P.PRT0;            % solo reloj local y PRT nominal
yyaxis left;
plot(tn, crudo*1e6, '.', 'MarkerSize', 3); ylabel('a_j - a_1 - m_j PRT (\mus)');
yyaxis right;
dh = S.dhat(:,3); plot(tn, (dh - dh(1))*1e6, '.', 'MarkerSize', 3);
ylabel('\delta_j - \delta_1 estimado (\mus)');
grid on; xlabel('t (s)');
legend('crudo (stagger + reloj)', 'tras asociar con el log del Tx', 'Location', 'southwest');
title(sprintf('(d) Llegadas del directo: asociacion k = %d, \\epsilon = %.2f ppm', S.k_hat, T.clk.eps*1e6));
guarda_fig(f, DIR_FIG, 'inf_senal');

% ======================= 3. Reloj por oscilador =========================
clases = {'OCXO', 'TCXO', 'XO', 'XO_BARATO'};
tr = (0:P.rx.M_slots-1).'*P.PRT0;
f = figure('Position', [50 50 1150 560]);
for ic = 1:numel(clases)
    o = osc_clase(clases{ic});
    rs = RandStream('mt19937ar', 'Seed', 40 + ic);
    w = [0; cumsum(o.sigma_y1*sqrt(P.PRT0)*randn(rs, numel(tr)-1, 1))];
    del = 0.6*o.eps_max*tr + 0.25*o.gamma_max*tr.^2 + w;
    V = [ones(size(tr)), tr, tr.^2];
    res = del - V*(V\del);                         % lo que no captura el modelo E1
    ph = 2*pi*fc*res;
    subplot(2, 4, ic);
    plot(tr, del*1e6, 'Color', col(1,:)); grid on;
    xlabel('t (s)'); ylabel('\delta(t) - \delta_0 (\mus)');
    title(sprintf('%s: \\epsilon = %.2g', strrep(clases{ic}, '_', ' '), 0.6*o.eps_max));
    subplot(2, 4, 4 + ic); hold on; box on;
    patch([0 tr(end) tr(end) 0], tol_phi*[-1 -1 1 1], [0.85 0.95 0.85], 'EdgeColor', 'none');
    plot(tr, ph, 'Color', col(2,:)); grid on;
    lim = max(1, 1.15*max(abs(ph))); ylim(lim*[-1 1]);
    xlabel('t (s)'); ylabel('2\pi f_c (\delta - modelo) (rad)');
    title(sprintf('residuo E1: %.2f rad rms', rms(ph)));
end
guarda_fig(f, DIR_FIG, 'inf_reloj');

% ======================= 4. Cortes de la IRF ============================
% Misma realizacion que exp2 (TCXO, escucha aleatoria 50 %, semilla 21)
P2 = P;
rs = RandStream('mt19937ar', 'Seed', 21);
P2.rx.obs = patron_escucha(P2.rx.M_slots, 'aleatorio', 0.5, rs);
[D2, T2, TX2] = simula_adquisicion(P2, 21);
S2 = sincroniza(D2, TX2, P2);
xs = (-3:0.01:3).';
Ga = [xs, zeros(size(xs)), zeros(size(xs))];
Gr = [zeros(size(xs)), xs, zeros(size(xs))];
Ia = forma_imagenes(D2, S2, T2, TX2, P2, Ga);
Ir = forma_imagenes(D2, S2, T2, TX2, P2, Gr);
ref = max(abs(Ia(:,5)));
nom = {'E0 offset', 'E1 modelo', 'E2 tiempo', 'E3 fase', 'reloj ideal'};
est = {'-', '-', '-', '-', 'k--'};
f = figure('Position', [50 50 1100 380]);
for p = 1:2
    subplot(1,2,p); hold on; grid on; box on;
    if p == 1, I = Ia; else, I = Ir; end
    for e = [1 2 3 4 5]
        if e < 5
            plot(xs, 20*log10(abs(I(:,e))/ref + 1e-6), est{e}, 'Color', col(e,:));
        else
            plot(xs, 20*log10(abs(I(:,e))/ref + 1e-6), est{e}, 'LineWidth', 1.5);
        end
    end
    yline(-1, 'r:'); ylim([-45 2]);
    if p == 1
        xlabel('x (m)  (corte por y = 0)'); title('(a) Corte en x del blanco central');
    else
        xlabel('y (m)  (corte por x = 0)'); title('(b) Corte en y del blanco central');
    end
    ylabel('amplitud relativa al pico ideal (dB)');
end
legend(nom, 'Location', 'south', 'NumColumns', 3);
guarda_fig(f, DIR_FIG, 'inf_cortes');

% ======================= 5. Conclusiones ================================
% Los paneles cuyo experimento aun no tiene resultados se dejan "pendiente"
hay = @(f) exist(fullfile(DIR_RES, f), 'file') == 2;
if hay('exp3_snr.csv'), A3 = readtable(fullfile(DIR_RES, 'exp3_snr.csv')); end
if hay('exp3_osc.csv'), B3 = readtable(fullfile(DIR_RES, 'exp3_osc.csv')); end
T4 = readtable(fullfile(DIR_RES, 'exp4_patrones.csv'), 'Delimiter', ',');
if hay('exp6_navegacion.csv'), T6 = readtable(fullfile(DIR_RES, 'exp6_navegacion.csv')); end
if hay('exp5_asociacion.csv'), T5 = readtable(fullfile(DIR_RES, 'exp5_asociacion.csv')); end
cls = {'OCXO', 'TCXO', 'XO', 'XO_BARATO'};

f = figure('Position', [50 50 1200 700]);
% (a) error de fase frente a SNR: teoria, medida y region de enfoque
subplot(2,3,1); hold on; box on;
snr = 0:0.5:60;
[ct, cp] = crb_directo(P, snr);
patch([0 60 60 0], [1e-3 1e-3 tol_phi tol_phi], [0.85 0.95 0.85], 'EdgeColor', 'none');
patch([0 17 17 0], [1e-3 1e-3 10 10], [0.92 0.92 0.92], 'EdgeColor', 'none');
h1 = semilogy(snr, 2*pi*fc*ct, '--', 'Color', col(2,:));
h2 = semilogy(snr, cp, '--', 'Color', col(4,:));
if hay('exp3_snr.csv')
    s3 = unique(A3.snr);
    m2 = arrayfun(@(s) median(A3.ep2(A3.snr == s)), s3);
    m3 = arrayfun(@(s) median(A3.ep3(A3.snr == s)), s3);
else
    s3 = NaN; m2 = NaN; m3 = NaN;
    text(30, 3, 'simulacion pendiente (exp3)', 'HorizontalAlignment', 'center', 'FontSize', 9);
end
h3 = semilogy(s3, m2, 'd', 'Color', col(2,:), 'MarkerFaceColor', col(2,:));
h4 = semilogy(s3, m3, 's', 'Color', col(4,:), 'MarkerFaceColor', col(4,:));
set(gca, 'YScale', 'log'); ylim([1e-3 10]); xlim([0 60]); grid on;
text(1, 2.5e-3, 'el directo no se distingue de los ecos', 'FontSize', 7, 'Rotation', 90);
text(40, 0.3, 'enfoca (< 1 dB)', 'Color', [0 0.45 0], 'FontSize', 8);
xlabel('SNR del directo por pulso (dB)'); ylabel('\sigma_\phi residual (rad)');
legend([h1 h2 h3 h4], {'E2 teoria 2\pi f_c CRB_t', 'E3 teoria CRB_\phi', 'E2 simulado', 'E3 simulado'}, ...
    'Location', 'northeast', 'FontSize', 7);
title('(a) Fase residual vs SNR_d');

% (b) oscilador: perdida y lobulos por estrategia
subplot(2,3,2);
if ~hay('exp3_osc.csv'), pendiente('(b) Lobulos vs oscilador', 'exp3');
else
med = @(v, c) median(v(strcmp(B3.clase, c)));
Mp = zeros(4,3); Ms = zeros(4,3);
for i = 1:4
    Mp(i,:) = [med(B3.perd2, cls{i}), med(B3.perd3, cls{i}), med(B3.perd4, cls{i})];
    Ms(i,:) = [med(B3.dpslr2, cls{i}), med(B3.dpslr3, cls{i}), med(B3.dpslr4, cls{i})];
end
b = bar(Mp); grid on;
for e = 1:3, b(e).FaceColor = col(e+1,:); end
yline(-1, 'r:');
set(gca, 'XTickLabel', strrep(cls, '_', ' '));
ylabel('perdida de pico frente a reloj ideal (dB)'); ylim([-32 1]);
legend('E1 modelo', 'E2 tiempo', 'E3 fase', 'Location', 'southwest');
title('(b) Imagen vs oscilador (SNR_d 30 dB)');
end

% (c) asociacion de pulsos (TCXO)
subplot(2,3,3);
if ~hay('exp5_asociacion.csv'), pendiente('(c) P(asociacion correcta), TCXO', 'exp5');
else
t5 = T5(strcmp(T5.clase, 'TCXO'), :);
Ss = unique(t5.stagger_s); Js = unique(t5.J);
Pm = zeros(numel(Ss), numel(Js));
for i = 1:numel(Ss)
    for j = 1:numel(Js)
        Pm(i,j) = t5.p_exito(t5.stagger_s == Ss(i) & t5.J == Js(j));
    end
end
imagesc(1:numel(Js), 1:numel(Ss), Pm, [0 1]); colormap(gca, parula); axis xy;
for i = 1:numel(Ss)
    for j = 1:numel(Js)
        if Pm(i,j) > 0.6, cc = 'k'; else, cc = 'w'; end
        text(j, i, sprintf('%.2f', Pm(i,j)), 'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', cc);
    end
end
sl = arrayfun(@(s) sprintf('%g ns', s*1e9), Ss, 'UniformOutput', false); sl{1} = '0 (PRF fija)';
set(gca, 'XTick', 1:numel(Js), 'XTickLabel', Js, 'YTick', 1:numel(Ss), 'YTickLabel', sl);
xlabel('ventanas observadas J'); ylabel('stagger S');
title('(c) P(asociacion correcta), TCXO');
end

% (d) patrones de escucha: sincronismo vs muestreo
subplot(2,3,4);
lab = strrep(strrep(strrep(T4.patron, 'periodico ', '1/'), 'aleatorio ', 'p='), 'rafagas ', '');
yyaxis left;
b = bar([T4.perd_E2, T4.perd_E3]); b(1).FaceColor = col(3,:); b(2).FaceColor = col(4,:);
ylim([-4 1.5]); ylabel('perdida de pico por sincronismo (dB)');
yline(-1, 'r:');
yyaxis right;
plot(1:height(T4), T4.muestreo_dPSLR, 'ko-', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
ylabel('\DeltaPSLR por muestreo de la apertura (dB)');
yl = ylim; ylim([min(yl(1), -1), 2*max(yl(2), 1)]);
set(gca, 'XTick', 1:height(T4), 'XTickLabel', lab, 'XTickLabelRotation', 50, 'FontSize', 8);
legend('E2', 'E3', '', 'muestreo', 'Location', 'north', 'NumColumns', 3, 'FontSize', 7);
grid on; title('(d) Patrones de escucha (TCXO, 30 dB)');

% (e) puenteo de fase entre ventanas sin directo
subplot(2,3,5);
T4b = readtable(fullfile(DIR_RES, 'exp4_puenteo.csv'));
hold on; box on; hh = gobjects(4,1);
for i = 1:4
    q = T4b(strcmp(T4b.clase, cls{i}), :);
    hh(i) = loglog(q.Tg_s, q.sigma_phi_sim, 'o-', 'Color', col(i,:), 'MarkerFaceColor', col(i,:), 'MarkerSize', 4);
    loglog(q.Tg_s, q.sigma_phi_teoria, '--', 'Color', col(i,:));
end
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on;
yline(tol_phi, 'r:', 'tolerancia', 'LabelHorizontalAlignment', 'left');
xlabel('hueco T_g sin directo (s)'); ylabel('\sigma_\phi interpolada (rad)');
legend(hh, strrep(cls, '_', ' '), 'Location', 'southeast', 'FontSize', 7);
title('(e) Puenteo de fase entre huecos');

% (f) navegacion del Tx
subplot(2,3,6); hold on; box on; grid on;
if ~hay('exp6_navegacion.csv'), pendiente('(f) Error de navegacion (Tx y Rx2)', 'exp6');
else
sg = T6.sigma_nav_m*100;
plot(sg, T6.perd_E2, 'd-', 'Color', col(3,:));
plot(sg, T6.perd_E3, 's-', 'Color', col(4,:), 'LineWidth', 1.8);
plot(sg, T6.perd_ideal_reloj, 'k--');
plot(sg, T6.perd_mono, 'o-', 'Color', col(1,:));
yline(-1, 'r:');
xlabel('\sigma_{nav} por eje, en cada dron (cm)'); ylabel('perdida de pico (dB)');
legend('bi E2', 'bi E3', 'bi reloj ideal', 'monostatico', 'Location', 'southwest', 'FontSize', 7);
title('(f) Error de navegacion (Tx y Rx2)');
end
guarda_fig(f, DIR_FIG, 'inf_conclusiones');
fprintf('fig_informe: listo\n');

function pendiente(tit, ex)
axis off; box on;
text(0.5, 0.5, sprintf('%s:\nsimulacion pendiente', ex), 'HorizontalAlignment', 'center', ...
    'FontSize', 11, 'Color', [0.4 0.4 0.4]);
title(tit);
end
