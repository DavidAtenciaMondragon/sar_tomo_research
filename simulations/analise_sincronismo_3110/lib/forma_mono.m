function I = forma_mono(M, TX, P, G)
% FORMA_MONO  BP monostatico del Radar 1 (su reloj es la referencia).
n = M.n; J = numel(n);
w = pesos_apertura(TX.pos(n,:));
I = bp_bistatico(M.seg, M.seg_t0, M.Ts_up, TX.pos(n,:), TX.pos(n,:), TX.t(n), ...
    zeros(J,1), zeros(J,1), w, G, P.fc, P.c);
end
