function [a_hat, r_hat, m0, ok] = estima_directo(y, g0, Ts, pref, umbral_dB)
% ESTIMA_DIRECTO  Deteccion y estimacion ML (retardo + fase) del camino directo.
%
% 1) Filtro adaptado (sin ponderar: estimador ML para un eco en ruido blanco).
% 2) Deteccion: primera muestra que supera umbral_dB sobre el ruido y es
%    maximo local en +-20 muestras (el directo es el primer camino en llegar;
%    los lobulos de rango del propio directo no son maximos en +-20 muestras).
%    Si nada supera el umbral se toma el maximo global (ok = false).
% 3) Refinamiento: Newton sobre |r(tau)|^2 con r(tau) evaluado por
%    interpolacion de Fourier exacta en una sub-puerta de 512 muestras.
%
% a_hat : llegada del directo en tiempo local (s); r_hat : r(a_hat) complejo
% m0    : indice (muestra) de la llegada gruesa dentro de y

Np = numel(pref);
Nf = 2^nextpow2(numel(y) + Np);
c = ifft(fft(y, Nf).*conj(fft(pref, Nf)));
c = c(1:numel(y));
pw = abs(c).^2;
sig2 = median(pw)/log(2);                     % potencia de ruido (robusta)
cand = find(pw >= sig2*10^(umbral_dB/10));
ok = false; m0 = [];
for i = cand(:).'
    lo = max(1, i-20); hi = min(numel(pw), i+20);
    if pw(i) >= max(pw(lo:hi)), m0 = i; ok = true; break; end
end
if isempty(m0), [~, m0] = max(pw); end

% --- refinamiento ML en sub-puerta -------------------------------------
Lsub = 512; Nfs = 1024;
i0 = max(1, m0 - 64); i0 = min(i0, max(1, numel(y) - Lsub + 1));
ys = y(i0:min(end, i0 + Lsub - 1));
X  = fft(ys, Nfs).*conj(fft(pref, Nfs));
wk = 2*pi*[0:Nfs/2-1, -Nfs/2:-1].'/Nfs;
tau = m0 - i0;                                 % en muestras (lag 0 = i0)
if m0 > 1 && m0 < numel(pw)                    % arranque parabolico
    a1 = abs(c(m0-1)); a2 = abs(c(m0)); a3 = abs(c(m0+1));
    den = a1 - 2*a2 + a3;
    if den < 0, tau = tau + 0.5*(a1 - a3)/den; end
end
for it = 1:8
    e  = exp(1j*wk*tau);
    r0 = sum(X.*e)/Nfs;
    r1 = sum(X.*(1j*wk).*e)/Nfs;
    r2 = sum(X.*(-(wk.^2)).*e)/Nfs;
    g  = real(conj(r0)*r1);
    gp = abs(r1)^2 + real(conj(r0)*r2);
    if gp >= 0, break; end                     % fuera de la cuenca del maximo
    step = -g/gp;
    step = max(min(step, 0.5), -0.5);
    tau = tau + step;
    if abs(step) < 1e-7, break; end
end
r_hat = sum(X.*exp(1j*wk*tau))/Nfs;
a_hat = g0 + (i0 - 1 + tau)*Ts;
end
