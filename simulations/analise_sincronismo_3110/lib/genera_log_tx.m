function [TX, pos_true, rx_true] = genera_log_tx(P, rs, rn)
% GENERA_LOG_TX  Log del Radar 1 (emision con stagger + navegacion) y
% navegacion del Radar 2 en los mismos instantes.
%   t_n = n*PRT0 + s_n, s_n ~ U(-S, S)   (S = P.stagger_max)
%   TX.pos / TX.rx son las navegaciones nominales de Tx y Rx2 (las que conoce
%   el procesador; la del Rx2 se registra con su reloj local y se lleva al
%   tiempo del Tx con la base de tiempos sincronizada: con delta conocido a
%   ns el error de posicion es despreciable);
%   pos_true / rx_true = nominal + error AR(1) (P.nav.sigma, P.nav.sigma_rx).
N  = P.N_tx;
s  = P.stagger_max*(2*rand(rs, N, 1) - 1);
tt = (0:N-1).'*P.PRT0 + s;
pos_nom = trayectoria(P, 'tx', tt);
rx_nom  = trayectoria(P, 'rx', tt);
pos_true = pos_nom + ar1(P.nav.sigma, N, P, rn);
rx_true  = rx_nom + ar1(P.nav.sigma_rx, N, P, rn);
TX.t = tt; TX.pos = pos_nom; TX.rx = rx_nom; TX.stagger = s;
end

function e = ar1(sigma, N, P, rn)
e = zeros(N,3);
if sigma <= 0, return; end
a = exp(-P.PRT0/P.nav.tcorr);
w = randn(rn, N, 3);
e(1,:) = sigma*w(1,:);
for k = 2:N
    e(k,:) = a*e(k-1,:) + sqrt(1-a^2)*sigma*w(k,:);
end
end
