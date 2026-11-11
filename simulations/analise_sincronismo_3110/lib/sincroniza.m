function S = sincroniza(D, TX, P)
% SINCRONIZA  Base de tiempos comun a partir del camino directo (post-proceso).
%
% Entradas: medidas del Rx2 (D) y log del Tx (TX). No usa verdad terreno.
% Salida S con las cuatro estrategias comparadas (columnas de dhat / psi):
%   E0 'offset'  : offset de tiempo y fase medidos en la 1a ventana, constantes
%   E1 'modelo'  : modelo global delta(t) = d0 + eps*t + gamma*t^2/2 (LS sobre
%                  los retardos del directo); fase = 2*pi*fc*delta(t) + phi0
%   E2 'tiempo'  : retardo del directo medido en cada ventana; fase derivada
%                  del retardo (2*pi*fc*delta_j + phi0), sin fase portadora
%   E3 'fase'    : retardo medido + fase portadora del directo en cada ventana
%                  (referenciado de fase por pulso)
%
%   delta_j medido = a_hat_j - (t_n + tau_d(n))
%   psi_j   medido = -(arg r_hat_j + 2*pi*fc*tau_d(n))     (mod 2*pi)
%
% Entradas: medidas del Rx2 (D) y logs de navegacion (TX.pos, TX.rx).

c = P.c; fc = P.fc;
tau_d = vecnorm(TX.pos - TX.rx, 2, 2)/c;        % navegacion nominal de ambos
t_pred = TX.t + tau_d;

[k_hat, ainfo] = asocia_pulsos(D.m, D.a_hat, t_pred);
S.k_hat = k_hat; S.asoc = ainfo;
if isnan(k_hat), k_hat = 0; end
n = k_hat + D.m + 1;
S.n = n;
tn = TX.t(n);
dm = D.a_hat - t_pred(n);                     % retardo de reloj medido
psim = wrap(-angle(D.r_hat) - 2*pi*mod(fc*tau_d(n), 1));

% --- E1: modelo polinomico (centrado, sobre dm - dm(1) por precision) ---
tc = mean(tn); x = tn - tc;
ord = min(2, numel(dm) - 1);
cf = polyfit(x, dm - dm(1), ord);
dpoly = polyval(cf, x);
d1 = dm(1) + dpoly;
ph1 = 2*pi*(mod(fc*dm(1), 1) + mod(fc*dpoly, 1));
phi0_1 = angle(sum(exp(1j*(psim - ph1))));
S.eps_hat = cf(end-1);                         % d(delta)/dt en tc ~ eps
S.cf = cf; S.tc = tc;

% --- E2: retardo por ventana, fase derivada del retardo ----------------
ph2 = 2*pi*mod(fc*dm, 1);
phi0_2 = angle(sum(exp(1j*(psim - ph2))));

J = numel(dm);
S.nombres = {'E0 offset', 'E1 modelo', 'E2 tiempo', 'E3 fase'};
S.dhat = [dm(1)*ones(J,1), d1, dm, dm];
S.psi  = [psim(1)*ones(J,1), wrap(ph1 + phi0_1), wrap(ph2 + phi0_2), psim];
S.dm = dm; S.psim = psim; S.tn = tn;
end

function y = wrap(x)
y = angle(exp(1j*x));
end
