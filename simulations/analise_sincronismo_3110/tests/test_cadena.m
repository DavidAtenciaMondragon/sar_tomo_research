% TEST_CADENA  Prueba de humo de la cadena completa (una realizacion).
%   Comprueba: asociacion correcta, errores de tiempo/fase coherentes con la
%   CRB, imagen ideal focalizada en la posicion verdadera, E3 ~ ideal.
clear; addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'lib'));
P = params_sinc();
P.rx.M_slots = 800; P.N_tx = 1400;            % arco de ~23 grados
tic; [D, T, TX] = simula_adquisicion(P, 1); fprintf('generacion: %.1f s\n', toc);
S = sincroniza(D, TX, P);
fprintf('k_hat = %d, k_true = %d, margen = %.1f\n', S.k_hat, T.n(1)-1, S.asoc.margen);
err_t = (S.dm - T.delta_d)*1e12;
wrap = @(x) angle(exp(1j*x));
err_p = wrap(S.psim - T.psi);
Np = round(P.Tp*P.fs); snr = 10^(P.SNR_d_dB/10);
crb_t = 1/(2*pi*P.B/sqrt(12)*sqrt(2*snr))*1e12; crb_p = 1/sqrt(2*snr);
fprintf('error tiempo: media %.2f ps, std %.2f ps (CRB %.2f ps)\n', mean(err_t), std(err_t), crb_t);
fprintf('error fase:   media %.4f rad, std %.4f rad (CRB %.4f rad)\n', mean(err_p), std(err_p), crb_p);
fprintf('eps verdadero %.3e, estimado %.3e\n', T.clk.eps, S.eps_hat);
R = rejilla_parches(P);
tic; I = forma_imagenes(D, S, T, TX, P, R.G); fprintf('BP: %.1f s\n', toc);
for e = 1:5
    Mt = metricas_imagen(I(:,e), I(:,5), R, P);
    fprintf('%-10s pos %.3f m  perdida %6.2f dB  PSLR %6.1f  ISLR %6.1f  rho %.3f\n', ...
        S.nombres{min(e,4)}, Mt.err_pos, Mt.perdida, Mt.PSLR, Mt.ISLR, Mt.rho);
end
n = R.n; A = reshape(abs(I(1:R.npix,5)).^2, n, n);
f = figure('Visible','off'); imagesc(R.u, R.u, 10*log10(A/max(A(:))), [-40 0]); axis xy image; colorbar;
exportgraphics(f, fullfile(tempdir, 'parche_ideal.png'));
fprintf('figura en %s\n', fullfile(tempdir, 'parche_ideal.png'));
