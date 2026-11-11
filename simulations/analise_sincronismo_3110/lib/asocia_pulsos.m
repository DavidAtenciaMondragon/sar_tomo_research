function [k_hat, info] = asocia_pulsos(m, a_hat, t_pred)
% ASOCIA_PULSOS  Asocia las ventanas del Rx2 con los pulsos del log del Tx.
%
%   Hipotesis k: la ventana de la ranura m_j contiene el pulso n_j = k + m_j.
%   Para cada k se ajusta por LS el modelo de reloj
%       a_j = alpha + beta*x_j + g*x_j^2,   x_j = t_pred(n_j) - t_pred(n_1)
%   (alpha ~ offset, beta ~ 1 + eps, g ~ deriva gamma/2; el termino
%   cuadratico solo con J >= 4, para que quede residuo) y se toma el k de
%   menor residuo. t_pred(n) = t_n + tau_d(n) es la llegada del directo
%   predicha con el log (instantes de emision + navegacion de Tx y Rx2).
%
%   Observabilidad: un desplazamiento de k pulsos cambia t_pred en k*PRT0
%   (absorbido por alpha) mas la diferencia de stagger s_{n+k}-s_n y de
%   geometria tau_d(n+k)-tau_d(n). Sin stagger y con tau_d suave (constante
%   en el vuelo circular comun), el residuo de las hipotesis erroneas es
%   casi identico al correcto: k es inobservable.
%
% info.rms    residuo rms de cada hipotesis (s)
% info.margen cociente entre el segundo mejor residuo y el mejor

m = m(:); a = a_hat(:) - a_hat(1); J = numel(m);
N = numel(t_pred);
kmax = N - 1 - max(m);
rms = inf(kmax+1, 1);
if J < 3
    % con 2 observaciones el ajuste afin es exacto: no hay informacion
    k_hat = NaN; info = struct('rms', rms, 'margen', 1, 'J', J); return;
end
ord = 1 + (J >= 4);
for k = 0:kmax
    x = t_pred(k + m + 1) - t_pred(k + 1);
    xs = x/max(abs(x));                         % escalado (condicionamiento)
    A = xs.^(0:ord);
    r = a - A*(A\a);
    rms(k+1) = sqrt(mean(r.^2));
end
[rs, io] = sort(rms);
k_hat = io(1) - 1;
info.rms = rms; info.margen = rs(2)/rs(1); info.J = J;
info.rms_mejor = rs(1); info.rms_segundo = rs(2);
end
