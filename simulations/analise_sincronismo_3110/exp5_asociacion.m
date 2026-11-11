% EXP5_ASOCIACION  Observabilidad del indice de pulso (asociacion Tx <-> Rx2).
%
% Modelo ligero (sin formar senales): llegadas del directo en el reloj del Rx2
%   a_j = t_n + tau_d(n) + delta(t) + ruido(sigma_t)
% con el log real del Tx (stagger + trayectoria) y el reloj del Rx2 completo
% (offset, frecuencia, deriva y paseo aleatorio). Se barre:
%   * stagger maximo S (0 = PRF uniforme; en el vuelo circular la distancia
%     Tx-Rx2 es constante y la geometria no aporta ninguna firma)
%   * numero de ventanas observadas J (aleatorias dentro de la vuelta)
%   * oscilador (TCXO, XO_barato)
% y se mide la probabilidad de asociacion correcta. sigma_t se toma de la
% CRB a la SNR del directo (validada en exp3).
% Salidas: resultados/exp5_asociacion.csv, figs/exp5_asociacion.*

clear; inicializa;
P0 = params_sinc();
Ss = [0 1e-9 1e-8 1e-7 1e-6 5e-6];
Js = [3 4 5 8 16 64 256];
clases = {'TCXO', 'XO_BARATO'};
Ntr = 200;
[sig_t, ~] = crb_directo(P0, P0.SNR_d_dB);
pexito = zeros(numel(Ss), numel(Js), numel(clases));
margen = pexito;
tic;
for ic = 1:numel(clases)
    o = osc_clase(clases{ic});
    for is = 1:numel(Ss)
        P = P0; P.stagger_max = Ss(is);
        for ij = 1:numel(Js)
            ok = 0; mg = zeros(Ntr,1);
            for tr = 1:Ntr
                rs = RandStream('mt19937ar', 'Seed', tr + 1e4*ij + 1e6*is + 1e8*ic);
                TX = genera_log_tx(P, rs, rs);
                tpred = TX.t + vecnorm(TX.pos - TX.rx, 2, 2)/P.c;
                n1 = 1 + randi(rs, P.N_tx - P.rx.M_slots);
                m = sort(randperm(rs, P.rx.M_slots - 1, Js(ij) - 1)).'; m = [0; m];
                n = n1 + m; t = tpred(n);
                epsl = o.eps_max*(2*rand(rs) - 1); gam = o.gamma_max*(2*rand(rs) - 1);
                w = [0; cumsum(o.sigma_y1*sqrt(diff(t)).*randn(rs, numel(t) - 1, 1))];
                a = t + 12.3 + epsl*t + 0.5*gam*t.^2 + w + sig_t*randn(rs, size(t));
                [kh, info] = asocia_pulsos(m, a, tpred);
                ok = ok + (kh == n1 - 1);
                mg(tr) = info.margen;
            end
            pexito(is, ij, ic) = ok/Ntr; margen(is, ij, ic) = median(mg);
        end
        fprintf('exp5 %s S = %g: %.0f s\n', clases{ic}, Ss(is), toc);
    end
end
fid = fopen(fullfile(DIR_RES, 'exp5_asociacion.csv'), 'w');
fprintf(fid, 'clase,stagger_s,J,p_exito,margen_mediano\n');
for ic = 1:numel(clases)
    for is = 1:numel(Ss)
        for ij = 1:numel(Js)
            fprintf(fid, '%s,%g,%d,%.3f,%.3g\n', clases{ic}, Ss(is), Js(ij), pexito(is,ij,ic), margen(is,ij,ic));
        end
    end
end
fclose(fid);
for ic = 1:numel(clases)
    fprintf('\nP(asociacion correcta), %s (filas: stagger; columnas: J = %s)\n', clases{ic}, mat2str(Js));
    for is = 1:numel(Ss)
        fprintf('  S = %8.0e s: %s\n', Ss(is), sprintf('%6.2f', pexito(is,:,ic)));
    end
end
save(fullfile(DIR_RES, 'exp5.mat'), 'pexito', 'margen', 'Ss', 'Js', 'clases');

f = figure('Position', [50 50 900 340]);
for ic = 1:numel(clases)
    subplot(1, 2, ic);
    semilogx(Js, pexito(:,:,ic).', 'o-'); grid on; ylim([0 1.05]);
    xlabel('ventanas observadas J'); ylabel('P(asociacion correcta)');
    legend(arrayfun(@(s) sprintf('S = %g ns', s*1e9), Ss, 'UniformOutput', false), 'Location', 'southeast');
    title(sprintf('Rx2 %s, \\sigma_t = %.0f ps', strrep(clases{ic}, '_', ' '), sig_t*1e12));
end
guarda_fig(f, DIR_FIG, 'exp5_asociacion');
