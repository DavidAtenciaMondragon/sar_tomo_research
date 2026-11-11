function I = bp_bistatico(seg, seg_t0, Ts_up, tx, rx, t_emit, dhat, psi, w, G, fc, c)
% BP_BISTATICO  Backprojection en el dominio del tiempo con base de tiempos corregida.
%
%   Para cada ventana j y pixel g:
%     tau = (|tx_j - g| + |g - rx_j|)/c
%     t_loc = t_emit_j + tau + dhat_j            (hora local del Rx2)
%     I(g) += w_j * seg_j(t_loc) * exp(j*(2*pi*fc*tau + psi_j))
%   dhat_j y psi_j son el error de reloj y la correccion de fase estimados.
%   dhat y psi pueden tener S columnas (S imagenes con la misma geometria).
%   Monostatico: rx = tx, dhat = 0, psi = 0.
%   I se normaliza por sum(w) (imagenes comparables entre patrones de escucha).

J = size(seg, 1); L = size(seg, 2); Sn = size(dhat, 2);
I = complex(zeros(size(G,1), Sn));
for j = 1:J
    tau = (vecnorm(G - tx(j,:), 2, 2) + vecnorm(G - rx(j,:), 2, 2))/c;
    ph = exp(1j*2*pi*fc*tau);
    sj = double(seg(j,:));
    for s = 1:Sn
        if s > 1 && dhat(j,s) == dhat(j,s-1)
            I(:,s) = I(:,s) + (w(j)*exp(1j*psi(j,s)))*v;   % mismo retardo: reutiliza
            continue;
        end
        idx = (t_emit(j) + dhat(j,s) - seg_t0(j))/Ts_up + tau/Ts_up + 1;
        i0 = floor(idx); fr = idx - i0;
        ok = i0 >= 1 & i0 < L;
        v = complex(zeros(size(idx)));
        v(ok) = ((1 - fr(ok)).*sj(i0(ok)).' + fr(ok).*sj(i0(ok)+1).').*ph(ok);
        I(:,s) = I(:,s) + (w(j)*exp(1j*psi(j,s)))*v;
    end
end
I = I/sum(w);
end
