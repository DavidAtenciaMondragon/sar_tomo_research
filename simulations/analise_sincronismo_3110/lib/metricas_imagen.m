function Mt = metricas_imagen(I, Iref, R, P)
% METRICAS_IMAGEN  Calidad de focalizacion de I frente a la imagen de referencia Iref.
%
% I, Iref : vectores de pixeles de los parches (rejilla_parches)
% Por blanco k:
%   err_pos  : distancia (m) entre el pico de I y el de Iref (interpolacion
%              parabolica): desplazamiento debido al sincronismo
%   err_abs  : distancia (m) entre el pico de I y el blanco verdadero
%   perdida  : 20log10(|pico I| / |pico Iref|) (dB, negativo = perdida)
%   PSLR     : maximo fuera del lobulo principal / pico (dB)
%   ISLR     : energia fuera / energia dentro del lobulo principal (dB)
%   rho      : |<I, Iref>| / (||I|| ||Iref||) en el parche (coherencia)
% El lobulo principal es la region conexa de Iref por encima de -20 dB que
% contiene su pico (con la IRF J0 de la apertura circular llega casi al
% primer nulo; los anillos quedan fuera), limitada a un radio P.ml_radio
% (si el fondo supera -20 dB la region no invade el parche), y se traslada al pico de I para no
% mezclar desplazamiento con desenfoque.
% Se devuelven medias sobre blancos (campos *) y peor caso (*_peor).

K = R.K; n = R.n; u = R.u;
ep = zeros(K,1); ea = ep; pl = ep; ps = ep; is = ep; rh = ep;
for k = 1:K
    ii = (k-1)*R.npix + (1:R.npix);
    A  = reshape(I(ii), n, n);  B = reshape(Iref(ii), n, n);
    pa = abs(A).^2; pb = abs(B).^2;
    [xa, ya, pka] = pico(pa, u);
    [xb, yb, pkb] = pico(pb, u);
    ep(k) = hypot(xa - xb, ya - yb);
    ea(k) = hypot(xa, ya);                       % blanco en el centro del parche
    pl(k) = 10*log10(pka/pkb);
    [X, Y] = meshgrid(u, u);
    Bi = interp2(X, Y, pb, xb + (X - xa), yb + (Y - ya), 'linear', 0);
    ml = region_conexa(Bi >= 0.01*pkb, X, Y, xa, ya) & hypot(X - xa, Y - ya) <= P.ml_radio;
    if all(ml(:)) || ~any(ml(:)), ml = hypot(X - xa, Y - ya) <= 0.4; end
    ps(k) = 10*log10(max(pa(~ml))/pka);
    is(k) = 10*log10(sum(pa(~ml))/sum(pa(ml)));
    rh(k) = abs(sum(A(:).*conj(B(:))))/sqrt(sum(pa(:))*sum(pb(:)));
end
Mt.err_pos = mean(ep);  Mt.err_pos_peor = max(ep);  Mt.err_abs = mean(ea);
Mt.perdida = mean(pl);  Mt.perdida_peor = min(pl);
Mt.PSLR = mean(ps);     Mt.ISLR = mean(is);
Mt.rho = mean(rh);      Mt.rho_peor = min(rh);
Mt.por_blanco = struct('err_pos', ep, 'err_abs', ea, 'perdida', pl, 'PSLR', ps, 'ISLR', is, 'rho', rh);
end

function ml = region_conexa(m, X, Y, x0, y0)
% Componente 4-conexa de la mascara m que contiene el pixel mas cercano a (x0, y0)
[~, i] = min(hypot(X(:) - x0, Y(:) - y0));
ml = false(size(m));
if ~m(i), return; end
ml(i) = true;
while true
    g = ml;
    g(2:end,:) = g(2:end,:) | ml(1:end-1,:); g(1:end-1,:) = g(1:end-1,:) | ml(2:end,:);
    g(:,2:end) = g(:,2:end) | ml(:,1:end-1); g(:,1:end-1) = g(:,1:end-1) | ml(:,2:end);
    g = g & m;
    if isequal(g, ml), break; end
    ml = g;
end
end

function [x0, y0, pk] = pico(p, u)
[pk, i] = max(p(:)); [iy, ix] = ind2sub(size(p), i);
d = u(2) - u(1);
x0 = u(ix) + d*parab(p, iy, ix, 2);
y0 = u(iy) + d*parab(p, iy, ix, 1);
end

function dlt = parab(p, iy, ix, dim)
dlt = 0;
if dim == 2, if ix == 1 || ix == size(p,2), return; end
    a = p(iy, ix-1); b = p(iy, ix); c = p(iy, ix+1);
else,        if iy == 1 || iy == size(p,1), return; end
    a = p(iy-1, ix); b = p(iy, ix); c = p(iy+1, ix);
end
a = log(a + eps); b = log(b + eps); c = log(c + eps);
den = a - 2*b + c;
if den < 0, dlt = 0.5*(a - c)/den; end
end
