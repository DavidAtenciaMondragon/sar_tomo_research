function [I, nombres] = forma_imagenes(D, S, T, TX, P, G, sel, extra)
% FORMA_IMAGENES  BP bistatico para las estrategias E0..E3 y el reloj ideal.
%
%   sel   : indices de ventanas de D a usar (por defecto todas)
%   extra : struct opcional con .dhat/.psi (J x Q) -> Q imagenes adicionales
%           con el indice de pulso verdadero (inyeccion de errores, exp1)
%   I     : Npix x (5+Q), columnas {E0, E1, E2, E3, ideal, extra...}
%   El "ideal" usa el error de reloj y la fase verdaderos (y el indice de
%   pulso verdadero): aisla la perdida debida al sincronismo de la debida al
%   ruido, al patron de escucha o a la navegacion.
if nargin < 7 || isempty(sel), sel = (1:numel(D.m)).'; end
if nargin < 8, extra = []; end
nombres = [S.nombres, {'Ideal'}];
Dt = [T.delta_d(sel), zeros(numel(sel), 0)];
Pt = T.psi(sel);
if ~isempty(extra), Dt = [Dt, extra.dhat(sel,:)]; Pt = [Pt, extra.psi(sel,:)]; end
if isequal(S.n(sel), T.n(sel))
    I = bp(D, sel, TX, P, G, T.n(sel), [S.dhat(sel,:), Dt], [S.psi(sel,:), Pt]);
else                                           % asociacion erronea
    I = [bp(D, sel, TX, P, G, S.n(sel), S.dhat(sel,:), S.psi(sel,:)), ...
         bp(D, sel, TX, P, G, T.n(sel), Dt, Pt)];
end
end

function I = bp(D, sel, TX, P, G, n, dh, ps)
I = bp_bistatico(D.seg(sel,:), D.seg_t0(sel), D.Ts_up, TX.pos(n,:), ...
    TX.rx(n,:), TX.t(n), dh, ps, pesos_apertura(TX.pos(n,:)), G, P.fc, P.c);
end
