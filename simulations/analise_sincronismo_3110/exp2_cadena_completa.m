% EXP2_CADENA_COMPLETA  Demostracion de extremo a extremo (una realizacion).
%
% Vuelo circular de ambos radares (una vuelta). Rx2 con TCXO, arranque en
% instante arbitrario, reloj con offset ~10 s, escucha discontinua aleatoria
% (50 % de las PRT). Se muestran:
%   * el error de reloj verdadero y lo que estima cada estrategia
%   * la funcion de coste de la asociacion de pulsos
%   * las imagenes BiSAR (E0..E3, ideal), la monostatica y la fusion coherente
% y se tabulan las metricas de sincronismo e imagen con el veredicto.
% Salidas: resultados/exp2_*.csv, figs/exp2_*.*

clear; inicializa;
P = params_sinc();
P.mono = true;
rs = RandStream('mt19937ar', 'Seed', 21);
P.rx.obs = patron_escucha(P.rx.M_slots, 'aleatorio', 0.5, rs);
tic; [D, T, TX, Mo] = simula_adquisicion(P, 21);
S = sincroniza(D, TX, P);
fprintf('exp2: generacion + sincronismo %.1f s, %d ventanas\n', toc, numel(D.m));
wrap = @(x) angle(exp(1j*x));
fc = P.fc; J = numel(D.m);
tn = TX.t(T.n);

% ---------------- Metricas de sincronismo --------------------------------
asoc_ok = isequal(S.n, T.n);
eT = S.dhat - T.delta_d;                       % J x 4 (s)
eP = wrap(S.psi - T.psi);                      % J x 4 (rad)
eP0 = wrap(eP - angle(mean(exp(1j*eP))));      % sin la constante (no desenfoca)
[crb_t, crb_p] = crb_directo(P, P.SNR_d_dB);
fprintf('asociacion: %s (k_hat %d, margen %.0f)\n', mat2str(asoc_ok), S.k_hat, S.asoc.margen);
fprintf('eps: verdadero %.4e  estimado %.4e  (error %.2e)\n', T.clk.eps, S.eps_hat, S.eps_hat - T.clk.eps);
fprintf('CRB: sigma_t %.0f ps, sigma_phi %.4f rad\n', crb_t*1e12, crb_p);

% ---------------- Imagenes ------------------------------------------------
R = rejilla_parches(P);
I = forma_imagenes(D, S, T, TX, P, R.G);
Im = forma_mono(Mo, TX, P, R.G);
u = -P.escena_semi:P.escena_d:P.escena_semi; [X, Y] = meshgrid(u, u);
Gs = [X(:), Y(:), zeros(numel(X),1)];
Is = forma_imagenes(D, S, T, TX, P, Gs);
Ims = forma_mono(Mo, TX, P, Gs);

nom = [S.nombres, {'Ideal'}];
fid = fopen(fullfile(DIR_RES, 'exp2_metricas.csv'), 'w');
fprintf(fid, 'estrategia,err_t_rms_ps,err_phi_rms_rad,err_pos_m,perdida_dB,PSLR_dB,ISLR_dB,rho,cumple\n');
Mid = metricas_imagen(I(:,5), I(:,5), R, P);
fprintf('\n%-10s %10s %10s %8s %8s %7s %7s %6s  %s\n', 'estrategia', 'e_t(ps)', 'e_phi', 'pos(m)', ...
    'perd(dB)', 'PSLR', 'ISLR', 'rho', 'veredicto');
for e = 1:5
    Mt = metricas_imagen(I(:,e), I(:,5), R, P);
    if e <= 4, et = rms(eT(:,e))*1e12; ep = rms(eP0(:,e)); else, et = 0; ep = 0; end
    ok = Mt.perdida_peor >= -P.crit.perdida_pico_dB && Mt.err_pos_peor <= P.crit.err_pos_m && ...
         Mt.PSLR - Mid.PSLR <= P.crit.dPSLR_dB && Mt.ISLR - Mid.ISLR <= P.crit.dISLR_dB && ...
         Mt.rho_peor >= P.crit.rho_min;
    v = {'NO', 'SI'};
    fprintf('%-10s %10.1f %10.4f %8.3f %8.2f %7.1f %7.1f %6.3f  %s\n', nom{e}, et, ep, ...
        Mt.err_pos, Mt.perdida, Mt.PSLR, Mt.ISLR, Mt.rho, v{ok+1});
    fprintf(fid, '%s,%.2f,%.5f,%.4f,%.3f,%.2f,%.2f,%.4f,%d\n', nom{e}, et, ep, Mt.err_pos, ...
        Mt.perdida, Mt.PSLR, Mt.ISLR, Mt.rho, ok);
end
fclose(fid);

% ---------------- Fusion coherente mono + bi -----------------------------
% Con referencia de fase por pulso (E3) la fase absoluta de la imagen bi
% queda anclada (phi0 se cancela) y es coherente con la monostatica.
fprintf('\nCoherencia mono-bi en los picos (fase bi - fase mono, grados):\n');
fus = zeros(R.K, 5); dph = zeros(R.K, 5);
for e = [2 3 4 5]
    for k = 1:R.K
        ii = (k-1)*R.npix + (1:R.npix);
        a = Im(ii)/max(abs(Im(ii))); b = I(ii,e)/max(abs(I(ii,5)));
        [~, ip] = max(abs(Im(ii)));
        dph(k,e) = rad2deg(angle(b(ip)*conj(a(ip))));
        fus(k,e) = max(abs(a + b))/(max(abs(a)) + max(abs(b)));
    end
    fprintf('  %-10s dphi = %s   ganancia fusion/max = %.3f\n', nom{e}, mat2str(round(dph(:,e).')), mean(fus(:,e)));
end
writematrix([dph(:,2:5), fus(:,2:5)], fullfile(DIR_RES, 'exp2_fusion.csv'));

% ---------------- Figura 1: reloj -----------------------------------------
f = figure('Position', [50 50 1100 650]);
subplot(2,2,1);
plot(tn, (T.delta_d - T.delta_d(1))*1e6, '.', 'MarkerSize', 3); grid on;
xlabel('t (s, reloj Tx)'); ylabel('\delta(t) - \delta(t_1) (\mus)');
title(sprintf('Error de reloj del Rx2 (\\epsilon = %.2f ppm, \\delta_0 = %.3f s)', T.clk.eps*1e6, T.clk.d0));
subplot(2,2,2);
plot(tn, eT(:,2)*1e9, '.', tn, eT(:,3)*1e9, '.', 'MarkerSize', 4); grid on;
xlabel('t (s)'); ylabel('error de tiempo (ns)'); legend('E1 modelo', 'E2/E3 por ventana');
title('Error de la base de tiempos estimada');
subplot(2,2,3);
plot(tn, eP0(:,2), '.', tn, eP0(:,3), '.', tn, eP0(:,4), '.', 'MarkerSize', 4); grid on;
yline(0.47, 'r:'); yline(-0.47, 'r:');
xlabel('t (s)'); ylabel('error de fase (rad)'); legend('E1 modelo', 'E2 tiempo', 'E3 fase', 'tolerancia 1 dB');
title('Error de la correccion de fase (sin constante)');
subplot(2,2,4);
k = 0:numel(S.asoc.rms)-1;
semilogy(k, S.asoc.rms*1e9, '-'); hold on; semilogy(T.n(1) - 1 - D.m(1), S.asoc.rms_mejor*1e9, 'ro');
grid on; xlabel('hipotesis k (indice del 1er pulso)'); ylabel('residuo rms (ns)');
title(sprintf('Asociacion de pulsos: k = %d, margen %.0f', S.k_hat, S.asoc.margen));
guarda_fig(f, DIR_FIG, 'exp2_reloj_asociacion');

% ---------------- Figura 2: imagenes --------------------------------------
f = figure('Position', [50 50 1200 600]);
ref = max(abs(Is(:,5)));
pan = {Is(:,1), Is(:,2), Is(:,3), Is(:,4), Is(:,5), Ims, ...
       Is(:,4)/ref + Ims/max(abs(Ims))};
tit = [nom, {'Monostatica (Radar 1)', 'Fusion mono + E3'}];
for p = 1:7
    subplot(2,4,p);
    A = abs(reshape(pan{p}, size(X)));
    if p <= 5, A = A/ref; else, A = A/max(A(:)); end
    imagesc(u, u, 20*log10(A + 1e-12), [-35 0]); axis xy image; colormap(gca, 'jet');
    hold on; plot(P.targets(:,1), P.targets(:,2), 'w+', 'MarkerSize', 5);
    title(tit{p}); xlabel('x (m)'); ylabel('y (m)');
end
subplot(2,4,8); axis off;
text(0, 0.9, sprintf('Rx2: %s, SNR_d = %d dB', P.osc.nombre, P.SNR_d_dB));
text(0, 0.75, sprintf('ventanas: %d de %d PRT (%.0f %%)', J, P.rx.M_slots, 100*J/P.rx.M_slots));
text(0, 0.6, sprintf('asociacion correcta: %s', mat2str(asoc_ok)));
text(0, 0.45, sprintf('\\sigma_t (E2) = %.0f ps', rms(eT(:,3))*1e12));
text(0, 0.3, sprintf('\\sigma_\\phi E2 = %.2f, E3 = %.3f rad', rms(eP0(:,3)), rms(eP0(:,4))));
text(0, 0.15, 'escala comun (dB) respecto al pico ideal');
guarda_fig(f, DIR_FIG, 'exp2_imagenes');
save(fullfile(DIR_RES, 'exp2.mat'), 'S', 'T', 'P', 'eT', 'eP', 'dph', 'fus');
