% EXP4_PULSOS_PERDIDOS  Efecto de la escucha discontinua del Rx2.
%
% (A) Patrones de escucha (periodico 1/G, aleatorio p, rafagas) sobre la misma
%     adquisicion continua. Se separan dos perdidas:
%       * perdida de SINCRONISMO : estrategia vs reloj ideal con el MISMO patron
%       * perdida de MUESTREO    : ideal con patron vs ideal continuo
%     (la primera es la que mide la calidad del sincronismo; la segunda es
%     inherente a tener menos pulsos y no se arregla sincronizando mejor).
% (B) Puenteo de fase: si una ventana NO tiene camino directo utilizable
%     (bloqueo, saturacion) su fase debe predecirse desde ventanas vecinas
%     separadas Tg. Error de fase vs Tg para cada oscilador y comparacion con
%     el puente browniano 2*pi*fc*sigma_y*sqrt(Tg/6).
% Salidas: resultados/exp4_patrones.csv, exp4_puenteo.csv, figs/exp4_*.*

clear; inicializa;
P0 = params_sinc();
R = rejilla_parches(P0);
wrap = @(x) angle(exp(1j*x));

% ======================= (A) patrones de escucha =========================
pat = {'continuo', [];  'periodico', 2; 'periodico', 4; 'periodico', 8; 'periodico', 16; ...
       'periodico', 32; 'periodico', 64; 'periodico', 128; 'aleatorio', 0.5; 'aleatorio', 0.2; ...
       'aleatorio', 0.05; 'aleatorio', 0.01; 'rafagas', [20 180]; 'rafagas', [5 495]};
nP = size(pat, 1); Nmc = 3;
etiq = cell(nP, 1);
for i = 1:nP
    if isempty(pat{i,2}), etiq{i} = pat{i,1}; else, etiq{i} = sprintf('%s %s', pat{i,1}, mat2str(pat{i,2})); end
end
cols = {'J', 'asoc_ok', 'perd_E1', 'perd_E2', 'perd_E3', 'pos_E3', 'dislr_E3', 'rho_E3', 'ok_E1', 'ok_E2', 'ok_E3', ...
        'muestreo_perd', 'muestreo_dPSLR', 'muestreo_dISLR', 'muestreo_rho', 'total_perd_E3'};
res = nan(nP, numel(cols), Nmc);
tic;
for mc = 1:Nmc
    P = P0;
    [D, T, TX] = simula_adquisicion(P, 4000 + mc);
    rs = RandStream('mt19937ar', 'Seed', 40 + mc);
    Ifull = [];
    for i = 1:nP
        obs = patron_escucha(P.rx.M_slots, pat{i,1}, pat{i,2}, rs);
        sel = find(obs(D.m + 1));
        [Ds, Ts] = subconjunto(D, T, sel);
        S = sincroniza(Ds, TX, P);
        I = forma_imagenes(Ds, S, Ts, TX, P, R.G);
        if i == 1, Ifull = I(:,5); end
        Mid = metricas_imagen(I(:,5), I(:,5), R, P);
        Mfull = metricas_imagen(Ifull, Ifull, R, P);
        r = zeros(1, numel(cols)); r(1) = numel(sel); r(2) = isequal(S.n, Ts.n);
        for e = 2:4
            Mt = metricas_imagen(I(:,e), I(:,5), R, P);
            ok = Mt.perdida_peor >= -P.crit.perdida_pico_dB && Mt.err_pos_peor <= P.crit.err_pos_m && ...
                 Mt.PSLR - Mid.PSLR <= P.crit.dPSLR_dB && Mt.ISLR - Mid.ISLR <= P.crit.dISLR_dB && ...
                 Mt.rho_peor >= P.crit.rho_min;
            r(1 + e) = Mt.perdida_peor; r(7 + e) = ok;
            if e == 4, r(6) = Mt.err_pos_peor; r(7) = Mt.ISLR - Mid.ISLR; r(8) = Mt.rho_peor; end
        end
        Ms = metricas_imagen(I(:,5), Ifull, R, P);
        r(12:15) = [Ms.perdida_peor, Mid.PSLR - Mfull.PSLR, Mid.ISLR - Mfull.ISLR, Ms.rho_peor];
        Mt = metricas_imagen(I(:,4), Ifull, R, P); r(16) = Mt.perdida_peor;
        res(i, :, mc) = r;
    end
    fprintf('exp4(A) realizacion %d/%d: %.0f s\n', mc, Nmc, toc);
end
med = median(res, 3);
Ta = array2table(med, 'VariableNames', cols); Ta.patron = etiq;
Ta = Ta(:, [end, 1:end-1]);
writetable(Ta, fullfile(DIR_RES, 'exp4_patrones.csv'));
fprintf('\n(A) Patrones de escucha (medianas de %d realizaciones; TCXO, SNR_d = %d dB)\n', Nmc, P0.SNR_d_dB);
disp(Ta(:, {'patron', 'J', 'asoc_ok', 'perd_E1', 'perd_E2', 'perd_E3', 'ok_E3', ...
            'muestreo_dPSLR', 'muestreo_dISLR', 'muestreo_rho'}));

% ======================= (B) puenteo de fase ==============================
clases = {'OCXO', 'TCXO', 'XO', 'XO_BARATO'};
Tg = 2.^(1:9);                                   % huecos en PRT (>= 4 anclas)
NmcB = 3;
sB = nan(numel(clases), numel(Tg), NmcB); slip = sB; teo = nan(numel(clases), numel(Tg));
[~, crb_p] = crb_directo(P0, P0.SNR_d_dB);
for ic = 1:numel(clases)
    P = P0; P.osc = osc_clase(clases{ic});
    sw = sqrt((2*pi*P.fc*P.osc.jit_rms)^2 + P.osc.phi_lo_rms^2 + crb_p^2);
    teo(ic,:) = sqrt((2*pi*P.fc*P.osc.sigma_y1)^2*Tg*P.PRT0/6 + sw^2*(1 + 2/3));
    for mc = 1:NmcB
        [D, T, TX] = simula_adquisicion(P, 5000 + 10*ic + mc);
        S = sincroniza(D, TX, P);
        tn = S.tn; J = numel(tn);
        for ig = 1:numel(Tg)
            a = 1:Tg(ig):J;
            if numel(a) < 4, continue; end
            x = tn - mean(tn(a));
            cf = polyfit(x(a), S.dm(a) - S.dm(1), 2);
            phm = 2*pi*(mod(P.fc*S.dm(1), 1) + mod(P.fc*polyval(cf, x), 1));
            ra = unwrap(wrap(S.psim(a) - phm(a)));
            psih = phm + interp1(tn(a), ra, tn, 'linear', 'extrap');
            nb = setdiff(1:J, a);
            e = wrap(psih(nb) - T.psi(nb));
            e = wrap(e - angle(mean(exp(1j*e))));
            sB(ic, ig, mc) = rms(e);
            slip(ic, ig, mc) = mean(abs(e) > pi/2);
        end
    end
    fprintf('exp4(B) %s: %.0f s\n', clases{ic}, toc);
end
sBm = median(sB, 3, 'omitnan'); slm = median(slip, 3, 'omitnan');
Tgs = Tg*P0.PRT0;
fid = fopen(fullfile(DIR_RES, 'exp4_puenteo.csv'), 'w');
fprintf(fid, 'clase,Tg_PRT,Tg_s,sigma_phi_sim,sigma_phi_teoria,frac_saltos\n');
for ic = 1:numel(clases)
    for ig = 1:numel(Tg)
        fprintf(fid, '%s,%d,%.5f,%.4f,%.4f,%.4f\n', clases{ic}, Tg(ig), Tgs(ig), sBm(ic,ig), teo(ic,ig), slm(ic,ig));
    end
end
fclose(fid);
fprintf('\n(B) Hueco maximo con sigma_phi <= 0.47 rad (perdida <= 1 dB):\n');
for ic = 1:numel(clases)
    o = osc_clase(clases{ic});
    tmax_teo = 6*max(0.47^2 - 0.03^2, 0)/(2*pi*P0.fc*o.sigma_y1)^2;
    v = sBm(ic,:); i = find(v > 0.47, 1);
    if isempty(i), s = sprintf('> %.3g s (limite del ensayo)', Tgs(find(~isnan(v), 1, 'last')));
    elseif i <= 2, s = '< 2 PRT';
    else, s = sprintf('%.3g s (%d PRT)', Tgs(i-1), Tg(i-1)); end
    fprintf('  %-10s simulado: %-20s teoria: %.3g s\n', clases{ic}, s, tmax_teo);
end
save(fullfile(DIR_RES, 'exp4.mat'), 'res', 'Ta', 'sB', 'slip', 'teo', 'Tg', 'clases');

% ---------------- Figuras -------------------------------------------------
f = figure('Position', [50 50 1100 380]);
ip = [1 2:8]; G = [1 2 4 8 16 32 64 128];
subplot(1,3,1);
semilogx(G, med(ip, 3), 'o-', G, med(ip, 4), 'd-', G, med(ip, 5), 's-', G, med(ip, 12), 'k^--'); grid on;
yline(-1, 'r:'); ylim([-15 1]); xlabel('escucha 1 de cada G PRT'); ylabel('perdida de pico (dB)');
legend('E1 vs ideal', 'E2 vs ideal', 'E3 vs ideal', 'ideal(patron) vs continuo', 'Location', 'southwest');
title('(A) Sincronismo vs muestreo');
subplot(1,3,2);
semilogx(G, med(ip, 13), 'o-', G, med(ip, 14), 's-', G, med(ip, 7), 'd-'); grid on;
xlabel('escucha 1 de cada G PRT'); ylabel('dB');
legend('\DeltaPSLR muestreo', '\DeltaISLR muestreo', '\DeltaISLR E3 vs ideal', 'Location', 'northwest');
title('(A) Lobulos secundarios');
subplot(1,3,3);
loglog(Tgs, sBm.', 'o-'); hold on; set(gca, 'ColorOrderIndex', 1); loglog(Tgs, teo.', '--');
yline(0.47, 'r:'); grid on; ylim([5e-3 5]);
xlabel('hueco entre ventanas con directo T_g (s)'); ylabel('\sigma_\phi interpolada (rad)');
legend([strrep(clases, '_', ' '), {'teoria'}], 'Location', 'northwest');
title('(B) Puenteo de fase entre huecos');
guarda_fig(f, DIR_FIG, 'exp4_pulsos_perdidos');
