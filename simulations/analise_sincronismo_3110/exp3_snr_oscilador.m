% EXP3_SNR_OSCILADOR  Precision ALCANZADA frente a SNR del directo y calidad del reloj.
%
% Monte Carlo de la cadena completa (escucha continua):
%   (A) SNR_d en {5..50} dB con TCXO          -> errores vs CRB, imagen
%   (B) clase de oscilador OCXO/TCXO/XO/XO_barato a SNR_d = 30 dB
% Para cada estrategia: error rms de tiempo/fase frente a la verdad, exito de
% la asociacion y de la deteccion, y metricas de imagen frente al reloj ideal.
% Salidas: resultados/exp3_snr.csv, exp3_osc.csv, figs/exp3_*.*

clear; inicializa;
Nmc = 4;
snrs = [5 10 15 20 25 30 40 50];
clases = {'OCXO', 'TCXO', 'XO', 'XO_BARATO'};
P0 = params_sinc();
R = rejilla_parches(P0);

casos = [num2cell(snrs(:)), repmat({'TCXO'}, numel(snrs), 1); ...
         repmat({30}, numel(clases), 1), clases(:)];
nC = size(casos, 1);
cols = {'snr','mc','det_ok','asoc_ok','et1','et2','ep1','ep2','ep3','deps', ...
        'pos1','pos2','pos3','pos4','perd1','perd2','perd3','perd4', ...
        'dpslr1','dpslr2','dpslr3','dpslr4','dislr1','dislr2','dislr3','dislr4', ...
        'rho1','rho2','rho3','rho4','ok1','ok2','ok3','ok4'};
res = nan(nC*Nmc, numel(cols)); clase = cell(nC*Nmc, 1);
wrap = @(x) angle(exp(1j*x));
fila = 0; tic;
for ic = 1:nC
    P = P0; P.SNR_d_dB = casos{ic,1}; P.osc = osc_clase(casos{ic,2});
    for mc = 1:Nmc
        fila = fila + 1;
        [D, T, TX] = simula_adquisicion(P, 1000*ic + mc);
        S = sincroniza(D, TX, P);
        I = forma_imagenes(D, S, T, TX, P, R.G);
        Mid = metricas_imagen(I(:,5), I(:,5), R, P);
        eT = S.dhat - T.delta_d; eP = wrap(S.psi - T.psi);
        eP = wrap(eP - angle(mean(exp(1j*eP))));
        r = [P.SNR_d_dB, mc, mean(D.det_ok), isequal(S.n, T.n), ...
             rms(eT(:,2)), rms(eT(:,3)), rms(eP(:,2)), rms(eP(:,3)), rms(eP(:,4)), ...
             S.eps_hat - T.clk.eps];
        m = zeros(1, 24);
        for e = 1:4
            Mt = metricas_imagen(I(:,e), I(:,5), R, P);
            ok = Mt.perdida_peor >= -P.crit.perdida_pico_dB && Mt.err_pos_peor <= P.crit.err_pos_m && ...
                 Mt.PSLR - Mid.PSLR <= P.crit.dPSLR_dB && Mt.ISLR - Mid.ISLR <= P.crit.dISLR_dB && ...
                 Mt.rho_peor >= P.crit.rho_min;
            m(e + [0 4 8 12 16 20]) = [Mt.err_pos_peor, Mt.perdida_peor, Mt.PSLR - Mid.PSLR, ...
                                       Mt.ISLR - Mid.ISLR, Mt.rho_peor, ok];
        end
        res(fila, :) = [r, m]; clase{fila} = casos{ic,2};
    end
    fprintf('exp3 caso %d/%d (%s, %d dB): %.0f s\n', ic, nC, casos{ic,2}, casos{ic,1}, toc);
end
Tb = array2table(res, 'VariableNames', cols); Tb.clase = clase;
nA = numel(snrs)*Nmc;
writetable(Tb(1:nA,:), fullfile(DIR_RES, 'exp3_snr.csv'));
writetable(Tb(nA+1:end,:), fullfile(DIR_RES, 'exp3_osc.csv'));
save(fullfile(DIR_RES, 'exp3.mat'), 'Tb', 'snrs', 'clases', 'Nmc');

% ---------------- Resumen por caso (medianas Monte Carlo) ----------------
agr = @(v) reshape(v, Nmc, []);
A = Tb(1:nA,:); B = Tb(nA+1:end,:);
[crb_t, crb_p] = crb_directo(P0, snrs);
fprintf('\n(A) TCXO, barrido de SNR del camino directo (medianas)\n');
fprintf('%5s %6s %5s %9s %8s %7s %7s | %7s %7s %7s | %s\n', 'SNR', 'det', 'asoc', 'et2(ps)', 'CRB_t', ...
    'ep2', 'ep3', 'perd1', 'perd2', 'perd3', 'cumple E1/E2/E3 (fraccion)');
for i = 1:numel(snrs)
    ii = (i-1)*Nmc + (1:Nmc);
    fprintf('%5d %6.3f %5.2f %9.0f %8.0f %7.3f %7.3f | %7.2f %7.2f %7.2f | %.2f %.2f %.2f\n', snrs(i), ...
        mean(A.det_ok(ii)), mean(A.asoc_ok(ii)), median(A.et2(ii))*1e12, crb_t(i)*1e12, ...
        median(A.ep2(ii)), median(A.ep3(ii)), median(A.perd2(ii)), median(A.perd3(ii)), ...
        median(A.perd4(ii)), mean(A.ok2(ii)), mean(A.ok3(ii)), mean(A.ok4(ii)));
end
fprintf('\n(B) SNR_d = 30 dB, clase de oscilador (medianas)\n');
fprintf('%10s %9s %8s %7s %7s | %7s %7s %7s | %s\n', 'clase', 'et1(ps)', 'deps', 'ep1', 'ep3', ...
    'perd1', 'perd2', 'perd3', 'cumple E1/E2/E3');
for i = 1:numel(clases)
    ii = (i-1)*Nmc + (1:Nmc);
    fprintf('%10s %9.0f %8.1e %7.3f %7.3f | %7.2f %7.2f %7.2f | %.2f %.2f %.2f\n', clases{i}, ...
        median(B.et1(ii))*1e12, median(abs(B.deps(ii))), median(B.ep1(ii)), median(B.ep3(ii)), ...
        median(B.perd2(ii)), median(B.perd3(ii)), median(B.perd4(ii)), ...
        mean(B.ok2(ii)), mean(B.ok3(ii)), mean(B.ok4(ii)));
end

% ---------------- Figuras -------------------------------------------------
f = figure('Position', [50 50 1100 380]);
subplot(1,3,1);
semilogy(snrs, median(agr(A.et2))*1e12, 'o-', snrs, crb_t*1e12, 'k--'); grid on;
yline(177, 'r:'); xlabel('SNR_d (dB)'); ylabel('\sigma_t (ps)');
legend('medido (E2/E3)', 'CRB', 'tolerancia 1 dB (exp1)'); title('Retardo del camino directo');
subplot(1,3,2);
semilogy(snrs, median(agr(A.ep2)), 'o-', snrs, median(agr(A.ep3)), 's-', snrs, crb_p, 'k--', ...
    snrs, 2*pi*P0.fc*crb_t, 'k:'); grid on;
yline(0.47, 'r:'); xlabel('SNR_d (dB)'); ylabel('\sigma_\phi (rad)');
legend('E2 (fase del retardo)', 'E3 (fase portadora)', 'CRB fase', '2\pi f_c CRB_t', 'tolerancia 1 dB');
title('Error de la correccion de fase');
subplot(1,3,3);
plot(snrs, median(agr(A.perd1)), 'o-', snrs, median(agr(A.perd2)), 'd-', snrs, median(agr(A.perd3)), ...
    's-', snrs, median(agr(A.perd4)), '^-'); grid on;
yline(-1, 'r:'); ylim([-25 1]); xlabel('SNR_d (dB)'); ylabel('perdida de pico peor blanco (dB)');
legend('E0', 'E1', 'E2', 'E3', 'Location', 'southeast'); title('Imagen BiSAR vs reloj ideal');
guarda_fig(f, DIR_FIG, 'exp3_snr');

f = figure('Position', [50 50 800 330]);
x = 1:numel(clases);
subplot(1,2,1);
semilogy(x, median(agr(B.ep1)), 'o-', x, median(agr(B.ep2)), 'd-', x, median(agr(B.ep3)), 's-'); grid on;
yline(0.47, 'r:'); set(gca, 'XTick', x, 'XTickLabel', strrep(clases, '_', ' '));
ylabel('\sigma_\phi (rad)'); legend('E1 modelo', 'E2 tiempo', 'E3 fase', 'Location', 'northwest');
title('Error de fase vs oscilador (SNR_d = 30 dB)');
subplot(1,2,2);
plot(x, median(agr(B.perd2)), 'o-', x, median(agr(B.perd3)), 'd-', x, median(agr(B.perd4)), 's-'); grid on;
yline(-1, 'r:'); set(gca, 'XTick', x, 'XTickLabel', strrep(clases, '_', ' '));
ylabel('perdida de pico (dB)'); title('Imagen vs oscilador');
guarda_fig(f, DIR_FIG, 'exp3_oscilador');
