function [p, v] = trayectoria(P, quien, t)
% TRAYECTORIA  Posicion y velocidad nominales en el vuelo circular comun.
%   quien = 'tx' (Radar 1) o 'rx' (Radar 2, adelantado P.dphi_rx en azimut)
%   t     : instantes (reloj del Tx), vector columna
t = t(:);
th = P.th0 + (P.v/P.radio)*t;
if strcmpi(quien, 'rx'), th = th + P.dphi_rx; end
p = [P.radio*cos(th), P.radio*sin(th), P.alt*ones(size(t))];
v = P.v*[-sin(th), cos(th), zeros(size(t))];
end
