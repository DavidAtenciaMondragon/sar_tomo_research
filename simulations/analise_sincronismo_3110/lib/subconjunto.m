function [Ds, Ts] = subconjunto(D, T, sel)
% SUBCONJUNTO  Ventanas sel de una adquisicion continua (simula un patron de
% escucha discontinuo: las ventanas son independientes entre si y el paseo
% aleatorio del reloj ya esta muestreado de forma consistente).
sel = sel(:);
Ds = D;
for f = {'m', 'a_hat', 'r_hat', 'seg_t0', 'det_ok'}
    Ds.(f{1}) = D.(f{1})(sel);
end
Ds.seg = D.seg(sel, :);
Ts = T;
for f = {'n', 'delta_d', 'psi'}
    Ts.(f{1}) = T.(f{1})(sel);
end
end
