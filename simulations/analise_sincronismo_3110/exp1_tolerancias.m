% EXP1_TOLERANCIAS  Precision de sincronismo NECESARIA para una imagen BiSAR enfocada.
%
% Se parte del dato bistatico real de la cadena y se forma la imagen con el
% reloj verdadero mas un error inyectado controlado. Cinco familias:
%   (a) offset de tiempo constante  Dtau           -> perdida de envolvente
%   (b) error de frecuencia residual Deps (lineal)  -> desenfoque (+ despl.)
%   (c) deriva cuadratica (QPE en el borde)         -> desenfoque
% En apertura circular completa la posicion la fija la fase (no la
% envolvente) y un error de fase de baja frecuencia solo se absorbe como
% desplazamiento en su componente sinusoidal con el azimut.
%   (d) ruido de fase blanco sigma_phi              -> perdida de pico / ISLR
%   (e) jitter de tiempo blanco sigma_t (envolvente + fase)
% y se compara con la teoria (seccion 3 de analise_teorica.md).
% Salidas: resultados/exp1_tolerancias.csv, figs/exp1_tolerancias.*

clear; inicializa;
P = params_sinc();
[D, T, TX] = simula_adquisicion(P, 11);
S = sincroniza(D, TX, P);
R = rejilla_parches(P);
J = numel(D.m); fc = P.fc; c = P.c;
tn = TX.t(T.n); tc = mean(tn); Ta = max(tn) - min(tn);
rs = RandStream('mt19937ar', 'Seed', 5);
Nrep = 4;

fam = {}; val = []; dh = zeros(J,0); ps = zeros(J,0);
add = @(fam, val, dh, ps, f, v, e, p) deal([fam, {f}], [val, v], [dh, T.delta_d + e], [ps, T.psi + p]);
for v = [0 0.5 1 2 5 10 20 50]*1e-9                    % (a)
    [fam, val, dh, ps] = add(fam, val, dh, ps, 'a_offset', v, v*ones(J,1), 2*pi*fc*v*ones(J,1));
end
for v = [1e-13 3e-13 1e-12 3e-12 1e-11 3e-11 1e-10]  % (b)
    e = v*(tn - tc);
    [fam, val, dh, ps] = add(fam, val, dh, ps, 'b_frecuencia', v, e, 2*pi*fc*e);
end
for v = [pi/8 pi/4 pi/2 3*pi/4 pi 3*pi/2 2*pi]         % (c)
    e = v/(2*pi*fc)*((tn - tc)/(Ta/2)).^2;
    [fam, val, dh, ps] = add(fam, val, dh, ps, 'c_cuadratico', v, e, 2*pi*fc*e);
end
for v = [0.05 0.1 0.2 0.3 0.5 0.7 1.0 1.4]             % (d)
    for r = 1:Nrep
        [fam, val, dh, ps] = add(fam, val, dh, ps, 'd_fase', v, zeros(J,1), v*randn(rs, J, 1));
    end
end
for v = [1 3 10 30 100 150 200 300 1000]*1e-12         % (e)
    for r = 1:Nrep
        e = v*randn(rs, J, 1);
        [fam, val, dh, ps] = add(fam, val, dh, ps, 'e_jitter', v, e, 2*pi*fc*e);
    end
end
tic;
I = forma_imagenes(D, S, T, TX, P, R.G, [], struct('dhat', dh, 'psi', ps));
fprintf('exp1: %d imagenes en %.0f s\n', size(I,2), toc);
Iid = I(:,5); I = I(:,6:end);

Q = numel(val);
met = struct('err_pos', zeros(Q,1), 'perdida', 0, 'PSLR', 0, 'ISLR', 0, 'rho', 0, 'err_c', 0);
for q = 1:Q
    Mt = metricas_imagen(I(:,q), Iid, R, P);
    met.err_pos(q) = Mt.err_pos; met.perdida(q) = Mt.perdida;
    met.PSLR(q) = Mt.PSLR; met.ISLR(q) = Mt.ISLR; met.rho(q) = Mt.rho;
    met.err_c(q) = Mt.por_blanco.err_pos(1);          % blanco central
end
Mid = metricas_imagen(Iid, Iid, R, P);

% ---------------- Teoria -------------------------------------------------
% Error de base de tiempos e_j (s) en el pulso j. Para el blanco central:
%   * envolvente: el BP muestrea el pulso comprimido (Hann) desplazado e_j
%     -> amplitud a_j = |r(e_j)|/|r(0)|
%   * fase: phi_j = 2*pi*fc*e_j. Un desplazamiento dp del pico produce
%     k*grad(Rb_j).dp, con grad(Rb_j) el gradiente sobre el suelo del rango
%     bistatico en el pulso j: se ajusta phi_j ~ c0 + k*grad(Rb_j).dp por LS
%     ponderado -> |dp| = desplazamiento previsto; el residuo da la perdida
%     |sum w_j a_j exp(j r_j)| / sum w_j.
w_ap = pesos_apertura(TX.pos(T.n,:));
pc = P.targets(1,:);
Gr = (pc - T.pos_true(T.n,:))./vecnorm(pc - T.pos_true(T.n,:), 2, 2) + ...
     (pc - T.rx_true(T.n,:))./vecnorm(pc - T.rx_true(T.n,:), 2, 2);
Gr = Gr(:,1:2);
env = envolvente(P);
teo = @(e) teoria_baja(2*pi*fc*e, env(e), w_ap, Gr, 2*pi*fc/c);
perd_fase = @(s) -10*log10(exp(1))*s.^2;               % 20log10(exp(-s^2/2))

% ---------------- Tabla --------------------------------------------------
[ufam, ~, ifam] = unique(fam, 'stable');
fid = fopen(fullfile(DIR_RES, 'exp1_tolerancias.csv'), 'w');
fprintf(fid, 'familia,valor,err_pos_m,err_pos_central_m,perdida_dB,PSLR_dB,ISLR_dB,rho,teoria,teoria_desp_m\n');
resumen = struct();
for f = 1:numel(ufam)
    ii = find(ifam == f); vv = unique(val(ii));
    for k = 1:numel(vv)
        jj = ii(val(ii) == vv(k)); th2 = NaN;
        switch ufam{f}
            case 'a_offset',     o = teo(vv(k)*ones(J,1)); th = o(2);
            case 'b_frecuencia', o = teo(vv(k)*(tn - tc)); th = o(2); th2 = o(1);
            case 'c_cuadratico', o = teo(vv(k)/(2*pi*fc)*((tn - tc)/(Ta/2)).^2); th = o(2); th2 = o(1);
            case 'd_fase',       th = perd_fase(vv(k));
            case 'e_jitter',     th = perd_fase(2*pi*fc*vv(k));
        end
        r = [vv(k), mean(met.err_pos(jj)), mean(met.err_c(jj)), mean(met.perdida(jj)), ...
             mean(met.PSLR(jj)), mean(met.ISLR(jj)), mean(met.rho(jj)), th, th2];
        fprintf(fid, '%s,%.6g,%.4f,%.4f,%.3f,%.2f,%.2f,%.4f,%.4f,%.4f\n', ufam{f}, r);
        resumen.(ufam{f})(k,:) = r;
    end
end
fclose(fid);

% ---------------- Umbrales (interpolados) --------------------------------
umb = @(x, y, lim) interp_umbral(x, y, lim);
A = resumen.a_offset; B = resumen.b_frecuencia; C = resumen.c_cuadratico;
Dd = resumen.d_fase;  E = resumen.e_jitter;
fprintf('\n=== Tolerancias (blanco central / media de blancos) ===\n');
fprintf('Imagen ideal: PSLR %.1f dB, ISLR %.1f dB\n', Mid.PSLR, Mid.ISLR);
fprintf('(a) offset de tiempo para perdida de %.1f dB: %.2f ns (desplazamiento max. %.3f m)\n', ...
    P.crit.perdida_pico_dB, 1e9*umb(A(:,1), -A(:,4), P.crit.perdida_pico_dB), max(A(:,2)));
bf = min(umb(B(:,1), -B(:,4), P.crit.perdida_pico_dB), umb(B(:,1), B(:,2), P.crit.err_pos_m));
fprintf('(b) error de frecuencia (perdida %.1f dB o despl. %.2f m): %.3g (= %.3g Hz en fc)\n', ...
    P.crit.perdida_pico_dB, P.crit.err_pos_m, bf, fc*bf);
fprintf('(c) QPE: perdida de %.1f dB con %.2f rad; desplazamiento %.2f m con %.2f rad\n', ...
    P.crit.perdida_pico_dB, umb(C(:,1), -C(:,4), P.crit.perdida_pico_dB), ...
    P.crit.err_pos_m, umb(C(:,1), C(:,2), P.crit.err_pos_m));
fprintf('(d) sigma_phi para perdida de %.1f dB: %.3f rad\n', P.crit.perdida_pico_dB, umb(Dd(:,1), -Dd(:,4), P.crit.perdida_pico_dB));
fprintf('(e) sigma_t para perdida de %.1f dB: %.0f ps\n', P.crit.perdida_pico_dB, 1e12*umb(E(:,1), -E(:,4), P.crit.perdida_pico_dB));
save(fullfile(DIR_RES, 'exp1_tolerancias.mat'), 'resumen', 'Mid', 'P');

% ---------------- Figura -------------------------------------------------
f = figure('Position', [50 50 1100 620]);
subplot(2,3,1);
semilogx(A(2:end,1)*1e9, A(2:end,4), 'o-', A(2:end,1)*1e9, A(2:end,8), 'k--'); grid on;
yline(-P.crit.perdida_pico_dB, 'r:'); xlabel('\Delta\tau constante (ns)'); ylabel('perdida de pico (dB)');
title('(a) offset de tiempo (solo envolvente)'); legend('simulado', 'teoria', 'criterio', 'Location', 'southwest');
subplot(2,3,2);
yyaxis left; semilogx(B(:,1), B(:,4), 'o-', B(:,1), B(:,8), 'k--'); ylabel('perdida de pico (dB)');
yline(-P.crit.perdida_pico_dB, 'r:');
yyaxis right; semilogx(B(:,1), B(:,2), 's-', B(:,1), B(:,9), 'k:'); ylabel('desplazamiento (m)');
grid on; xlabel('\Delta\epsilon residual'); title('(b) frecuencia residual');
legend('perdida sim.', 'teoria', '', 'despl. sim.', 'teoria', 'Location', 'west');
subplot(2,3,3);
yyaxis left; plot(C(:,1)/pi, C(:,4), 'o-', C(:,1)/pi, C(:,8), 'k--'); ylabel('perdida de pico (dB)');
yline(-P.crit.perdida_pico_dB, 'r:');
yyaxis right; plot(C(:,1)/pi, C(:,2), 's-', C(:,1)/pi, C(:,9), 'k:'); ylabel('desplazamiento (m)');
grid on; xlabel('QPE en el borde (\times\pi rad)'); title('(c) deriva cuadratica');
legend('perdida sim.', 'teoria', '', 'despl. sim.', 'teoria', 'Location', 'west');
subplot(2,3,4);
plot(Dd(:,1), Dd(:,4), 'o-', Dd(:,1), Dd(:,8), 'k--', Dd(:,1), Dd(:,6) - Mid.ISLR, 's-'); grid on;
yline(-P.crit.perdida_pico_dB, 'r:'); xlabel('\sigma_\phi (rad)'); ylabel('dB');
title('(d) ruido de fase blanco'); legend('perdida pico', 'teoria e^{-\sigma^2/2}', '\DeltaISLR', 'Location', 'southwest');
subplot(2,3,5);
semilogx(E(:,1)*1e12, E(:,4), 'o-', E(:,1)*1e12, E(:,8), 'k--'); grid on;
yline(-P.crit.perdida_pico_dB, 'r:'); xlabel('\sigma_t (ps)'); ylabel('perdida de pico (dB)');
title('(e) jitter de tiempo'); ylim([-30 1]);
subplot(2,3,6);
yyaxis left; semilogx(E(:,1)*1e12, E(:,6) - Mid.ISLR, 's-'); ylabel('\DeltaISLR (dB)');
yyaxis right; semilogx(E(:,1)*1e12, E(:,7), 'o-'); ylabel('\rho (coherencia con ideal)');
grid on; xlabel('\sigma_t (ps)'); title('(e) ISLR y coherencia');
sgtitle(sprintf('Tolerancias de sincronismo (f_c = %.0f MHz, B = %.0f MHz)', fc/1e6, P.B/1e6));
guarda_fig(f, DIR_FIG, 'exp1_tolerancias');

function out = teoria_baja(phi, amp, w, Gr, k)
% [desplazamiento (m), perdida (dB)] previstos para un error de fase phi_j
% (rad) y una atenuacion de envolvente amp_j por pulso
sw = sqrt(w);
A = [ones(size(phi)), k*Gr];
x = (sw.*A) \ (sw.*phi);
r = phi - A*x;
out = [norm(x(2:3)), 20*log10(abs(sum(w.*amp.*exp(1j*r)))/sum(w))];
end

function f = envolvente(P)
% |r(dt)|/|r(0)| del pulso comprimido con ventana Hann (como rc_segmento)
Ts = 1/P.fs; Np = round(P.Tp*P.fs); Nf = 1024; U = 64;
pref = chirp_bb((0:Np-1).'*Ts, P.Tp, P.B/P.Tp);
fr = [0:Nf/2-1, -Nf/2:-1].'/(Nf*Ts);
W = 0.5 + 0.5*cos(2*pi*fr/P.B); W(abs(fr) > P.B/2) = 0;
X = abs(fft(pref, Nf)).^2.*W;
Xu = zeros(U*Nf, 1); Xu(1:Nf/2) = X(1:Nf/2); Xu(end-Nf/2+1:end) = X(Nf/2+1:end);
r = abs(fftshift(ifft(Xu)))*U*Nf/sum(X);             % tabla con paso Ts/U
tg = ((0:U*Nf-1).' - U*Nf/2)*Ts/U;
f = @(dt) interp1(tg, r, max(min(dt(:), tg(end)), tg(1)));
end

function x = interp_umbral(v, y, lim)
% Primer cruce de y por lim (interpolacion en log del parametro)
x = NaN;
i = find(y >= lim, 1);
if isempty(i), return; end
if i == 1, x = v(1); return; end
t = (lim - y(i-1))/(y(i) - y(i-1));
x = exp(log(max(v(i-1), eps)) + t*(log(v(i)) - log(max(v(i-1), eps))));
end
