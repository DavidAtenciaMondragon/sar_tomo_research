% EXP6_NAVEGACION  Acoplamiento entre error de navegacion y sincronismo.
%
% Ambos drones vuelan con un error de navegacion AR(1) (tau_c = 1 s) de la
% misma sigma por eje, independiente en cada uno. El camino directo mide
% reloj + distancia Tx-Rx2 a la vez: la parte del error de ambos proyectada
% sobre la linea Tx-Rx2 es indistinguible de un error de reloj y se
% "corrige" como si lo fuera. Para el eco de un blanco P el error de fase
% que queda es k*[(u_TxP - u_TxRx).eTx + (u_RxP + u_TxRx).eRx] en lugar de
% k*(u_TxP.eTx + u_RxP.eRx). Se compara la imagen BiSAR (E2, E3, reloj
% ideal) y la monostatica (solo error del Tx) frente a las mismas imagenes
% sin error de navegacion (mismo ruido).
% Salidas: resultados/exp6_navegacion.csv, figs/exp6_navegacion.*

clear; inicializa;
P0 = params_sinc(); P0.mono = true;
R = rejilla_parches(P0);
sig = [0 0.005 0.01 0.02 0.05 0.1];
Nmc = 3;
out = nan(numel(sig), 4, Nmc);                 % perdida peor: E2, E3, ideal, mono
rho = out;
tic;
for mc = 1:Nmc
    P = P0; P.nav.sigma = 0; P.nav.sigma_rx = 0;
    [D, T, TX, Mo] = simula_adquisicion(P, 6000 + mc);
    S = sincroniza(D, TX, P);
    I0 = forma_imagenes(D, S, T, TX, P, R.G); Im0 = forma_mono(Mo, TX, P, R.G);
    for is = 1:numel(sig)
        P.nav.sigma = sig(is); P.nav.sigma_rx = sig(is);
        [D, T, TX, Mo] = simula_adquisicion(P, 6000 + mc);
        S = sincroniza(D, TX, P);
        I = forma_imagenes(D, S, T, TX, P, R.G); Im = forma_mono(Mo, TX, P, R.G);
        pares = {I(:,3), I0(:,5); I(:,4), I0(:,5); I(:,5), I0(:,5); Im, Im0};
        for q = 1:4
            Mt = metricas_imagen(pares{q,1}, pares{q,2}, R, P);
            out(is, q, mc) = Mt.perdida_peor; rho(is, q, mc) = Mt.rho_peor;
        end
    end
    fprintf('exp6 realizacion %d/%d: %.0f s\n', mc, Nmc, toc);
end
om = median(out, 3); rm = median(rho, 3);
lam = P0.c/P0.fc;
fid = fopen(fullfile(DIR_RES, 'exp6_navegacion.csv'), 'w');
fprintf(fid, 'sigma_nav_m,perd_E2,perd_E3,perd_ideal_reloj,perd_mono,rho_E2,rho_E3,rho_ideal,rho_mono\n');
for is = 1:numel(sig)
    fprintf(fid, '%.3f,%.3f,%.3f,%.3f,%.3f,%.4f,%.4f,%.4f,%.4f\n', sig(is), om(is,:), rm(is,:));
end
fclose(fid);
fprintf('\nPerdida de pico (peor blanco, dB) frente a navegacion perfecta (lambda = %.2f m)\n', lam);
fprintf('%9s %8s %8s %10s %8s\n', 'sigma(m)', 'E2', 'E3', 'reloj id.', 'mono');
for is = 1:numel(sig)
    fprintf('%9.3f %8.2f %8.2f %10.2f %8.2f\n', sig(is), om(is,:));
end
save(fullfile(DIR_RES, 'exp6.mat'), 'out', 'rho', 'sig');

f = figure('Position', [50 50 800 330]);
subplot(1,2,1);
plot(sig*100, om, 'o-'); grid on; yline(-1, 'r:');
xlabel('\sigma_{nav} por eje y por dron (cm)'); ylabel('perdida de pico (dB)');
legend('Bi E2', 'Bi E3 (fase directo)', 'Bi reloj ideal', 'Mono', 'Location', 'southwest');
title('Error de navegacion de Tx y Rx2 (AR(1), \tau_c = 1 s)');
subplot(1,2,2);
plot(sig*100, rm, 'o-'); grid on; xlabel('\sigma_{nav} (cm)'); ylabel('\rho peor blanco');
title('Coherencia con la imagen sin error');
guarda_fig(f, DIR_FIG, 'exp6_navegacion');
