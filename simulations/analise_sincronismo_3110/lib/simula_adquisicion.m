function [D, T, TX, M] = simula_adquisicion(P, seed)
% SIMULA_ADQUISICION  Genera y procesa (en el lado del Rx2) una adquisicion.
%
%   [D, T, TX, M] = simula_adquisicion(P, seed)
%
% Salidas
%   TX : log del Radar 1 (lo unico que se intercambia en post-proceso)
%        .t (N x 1) instantes de emision, .pos (N x 3) navegacion nominal,
%        .rx (N x 3) navegacion nominal del Rx2 en esos instantes
%   D  : medidas del Rx2, una fila por ventana observada (solo usa datos
%        disponibles en el receptor: muestras + marcas de su reloj local)
%        .m      ranura de escucha (0 = ventana de adquisicion)
%        .a_hat  llegada del directo en tiempo local (s)
%        .r_hat  pico complejo del filtro adaptado del directo
%        .seg    J x L datos comprimidos (sobremuestreados xU) tras cancelar
%                el directo; .seg_t0 tiempo local de su primera muestra
%   T  : verdad terreno (solo para evaluar): .n indice de pulso, .delta_d
%        error de reloj en la llegada del directo, .psi correccion de fase
%        exacta, parametros del reloj y posiciones verdaderas del Tx
%   M  : dato monostatico del Radar 1 (si P.mono), mismo formato que D
%
% Modelo de senal (banda base, "stop-and-hop", intra-pulso despreciable):
%   y2(t2) = sum_p A_p * s(t2 - a_p) * exp(-j(2*pi*fc*(tau_p + delta_p) + phi0 + phi_lo))
%   a_p = t_n + tau_p + delta_p     (llegada en tiempo local del Rx2)

rs = RandStream('mt19937ar', 'Seed', seed);
rn = RandStream('mt19937ar', 'Seed', seed + 7777);   % navegacion (independiente)

c = P.c; fc = P.fc; Ts = 1/P.fs; U = P.U;
Np = round(P.Tp*P.fs); K = P.B/P.Tp;
pref = chirp_bb((0:Np-1).'*Ts, P.Tp, K);

% ---------------- Log del Tx ------------------------------------------
[TX, pos_true, rx_true] = genera_log_tx(P, rs, rn);
N = P.N_tx; tt = TX.t; pos_nom = TX.pos;

% ---------------- Reloj del Rx2 ---------------------------------------
o = P.osc;
clk.eps   = o.eps_max*(2*rand(rs) - 1);
clk.gamma = o.gamma_max*(2*rand(rs) - 1);
clk.phi0  = 2*pi*rand(rs);
clk.d0    = P.delta0_rango(1) + diff(P.delta0_rango)*rand(rs);
ddet = @(t) clk.d0 + clk.eps*t + 0.5*clk.gamma*t.^2;

Mslots = P.rx.M_slots;
obs = P.rx.obs;
if isempty(obs), obs = true(Mslots,1); end
obs = logical(obs(:)); obs(1) = true;              % la adquisicion siempre
mlist = find(obs) - 1;  J = numel(mlist);

% Ventana de adquisicion: el Rx2 graba una PRT completa desde t_rx0
t_rx0 = P.rx.t_start + P.PRT0*rand(rs);
td_all = tt + vecnorm(pos_true - rx_true, 2, 2)/c;     % llegada verdadera del directo
n1 = find(td_all >= t_rx0 + 2e-6, 1);
assert(~isempty(n1) && n1 + Mslots - 1 <= N, 'N_tx insuficiente para la escucha');
L0 = t_rx0 + ddet(t_rx0);                           % hora local de la 1a muestra

% Semiancho de la puerta de seguimiento (calculable por el receptor): se
% centra en la ultima llegada medida del directo + m*PRT0 y cubre el stagger
% y la deriva maxima del reloj durante el hueco Tg desde esa medida.
gpre = @(Tg) 2*P.stagger_max + o.eps_max*Tg + 0.5*o.gamma_max*Tg^2 + ...
       5*o.sigma_y1*sqrt(Tg) + 3e-6;
gpost = (P.seg_post + Np)*Ts + 1e-6;

Ad = sqrt(10^(P.SNR_d_dB/10)/Np);
Ae = sqrt(10^(P.SNR_e_dB/10)/Np);
Ktg = size(P.targets,1);

Lseg = (P.seg_pre + P.seg_post)*U + 1;
D.m = mlist; D.a_hat = zeros(J,1); D.r_hat = zeros(J,1);
D.seg = complex(zeros(J, Lseg, 'single')); D.seg_t0 = zeros(J,1);
D.Ts_up = Ts/U; D.det_ok = true(J,1);
T.n = zeros(J,1); T.delta_d = zeros(J,1); T.psi = zeros(J,1);

w_prev = 0; t_prev = NaN;
for j = 1:J
    m = mlist(j); n = n1 + m;
    txp = pos_true(n,:); rxp = rx_true(n,:);
    tau_d = norm(txp - rxp)/c;
    tau_e = (vecnorm(P.targets - txp, 2, 2) + vecnorm(P.targets - rxp, 2, 2))/c;
    t_d = tt(n) + tau_d;
    % paseo aleatorio (FM blanca) en las llegadas, incl. dentro de la ventana
    if j == 1, w_d = 0; else, w_d = w_prev + o.sigma_y1*sqrt(t_d - t_prev)*randn(rs); end
    w_prev = w_d; t_prev = t_d;
    t_e = tt(n) + tau_e;
    [ts_e, ord] = sort(t_e);
    w_e = zeros(Ktg,1);
    w_e(ord) = w_d + cumsum(o.sigma_y1*sqrt(diff([t_d; ts_e])).*randn(rs, Ktg, 1));
    jw  = o.jit_rms*randn(rs);
    plo = o.phi_lo_rms*randn(rs);
    del_d = ddet(t_d) + w_d + jw;
    del_e = ddet(t_e) + w_e + jw;
    a_p   = [t_d + del_d; t_e + del_e];
    tau_p = [tau_d; tau_e]; del_p = [del_d; del_e];
    ph_p  = -2*pi*(mod(fc*tau_p,1) + mod(fc*del_p,1)) - clk.phi0 - plo;
    amp_p = [Ad; Ae*ones(Ktg,1)].*exp(1j*ph_p);

    % --- puerta de grabacion (solo informacion del receptor) -----------
    if j == 1
        g0 = L0; Nw = ceil((P.PRT0 + 3e-6)/Ts);
    else
        Tg = (m - mlist(j-1))*P.PRT0; Gpre = gpre(Tg);
        cm = D.a_hat(j-1) + (m - mlist(j-1))*P.PRT0;
        g0 = L0 + round((cm - Gpre - L0)/Ts)*Ts;  Nw = ceil((2*Gpre + gpost)/Ts);
    end
    y = sintetiza(g0, Nw, Ts, a_p, amp_p, P.Tp, K, rs);

    % --- camino directo: deteccion + ML fino de retardo y fase ----------
    [ah, rh, m0, ok] = estima_directo(y, g0, Ts, pref, P.det_umbral_dB);
    D.a_hat(j) = ah; D.r_hat(j) = rh; D.det_ok(j) = ok;
    % cancelacion del directo y compresion para BP
    y = y - (rh/Np)*sintetiza_sin_ruido(g0, Nw, Ts, ah, P.Tp, K);
    [D.seg(j,:), D.seg_t0(j)] = rc_segmento(y, g0, Ts, pref, m0, P);

    T.n(j) = n; T.delta_d(j) = del_d;
    T.psi(j) = angle(exp(1j*(2*pi*mod(fc*del_d,1) + clk.phi0 + plo)));
end
T.clk = clk; T.pos_true = pos_true; T.rx_true = rx_true; T.n1 = n1;

% ---------------- Dato monostatico (Radar 1, reloj de referencia) ------
M = [];
if P.mono
    Am = sqrt(10^(P.SNR_m_dB/10)/Np);
    nlist = n1 + (0:Mslots-1).';
    M.n = nlist; Jm = numel(nlist);
    M.seg = complex(zeros(Jm, Lseg, 'single')); M.seg_t0 = zeros(Jm,1); M.Ts_up = Ts/U;
    for j = 1:Jm
        n = nlist(j); txp = pos_true(n,:);
        Rmin = min(vecnorm(P.targets - pos_nom(n,:), 2, 2));   % de la navegacion propia
        tau_m = 2*vecnorm(P.targets - txp, 2, 2)/c;
        amp = Am*exp(-1j*2*pi*mod(fc*tau_m,1));
        g0 = tt(n) + round((2*Rmin/c - 2e-6)/Ts)*Ts;     % rejilla del propio Tx
        Nw = 2*(P.seg_pre + P.seg_post) + Np + 400;
        y = sintetiza(g0, Nw, Ts, tt(n) + tau_m, amp, P.Tp, K, rs);
        m0 = round(2e-6/Ts) - 20 + 1;                    % "directo" ficticio
        [M.seg(j,:), M.seg_t0(j)] = rc_segmento(y, g0, Ts, pref, m0, P);
    end
end
end

% =====================================================================
function y = sintetiza(g0, Nw, Ts, a, amp, Tp, K, rs)
y = (randn(rs, Nw, 1) + 1j*randn(rs, Nw, 1))/sqrt(2);
for p = 1:numel(a)
    i1 = max(1, ceil((a(p) - g0)/Ts) + 1);
    i2 = min(Nw, floor((a(p) + Tp - g0)/Ts) + 1);
    if i2 < i1, continue; end
    ii = (i1:i2).';
    y(ii) = y(ii) + amp(p)*chirp_bb((g0 - a(p)) + (ii-1)*Ts, Tp, K);
end
end

function y = sintetiza_sin_ruido(g0, Nw, Ts, a, Tp, K)
y = zeros(Nw,1);
i1 = max(1, ceil((a - g0)/Ts) + 1);
i2 = min(Nw, floor((a + Tp - g0)/Ts) + 1);
ii = (i1:i2).';
y(ii) = chirp_bb((g0 - a) + (ii-1)*Ts, Tp, K);
end

function [seg, t0] = rc_segmento(y, g0, Ts, pref, m0, P)
% Compresion en rango con ponderacion Hann, sobremuestreo xU por relleno
% espectral y extraccion de [m0-seg_pre, m0+seg_post] muestras.
U = P.U; Lsub = 512; Nf = 1024;
i0 = max(1, m0 - 64); i0 = min(i0, max(1, numel(y) - Lsub + 1));
ys = y(i0:min(end, i0 + Lsub - 1));
f = [0:Nf/2-1, -Nf/2:-1].'/(Nf*Ts);
W = 0.5 + 0.5*cos(2*pi*f/P.B); W(abs(f) > P.B/2) = 0;
X = fft(ys, Nf).*conj(fft(pref, Nf)).*W;
Xu = complex(zeros(U*Nf,1));
Xu(1:Nf/2) = X(1:Nf/2); Xu(end-Nf/2+1:end) = X(Nf/2+1:end);
r = ifft(Xu)*U;
l0 = (m0 - i0 - P.seg_pre)*U + 1;
idx = l0 + (0:(P.seg_pre + P.seg_post)*U);
ok = idx >= 1 & idx <= numel(r);
seg = complex(zeros(1, numel(idx), 'single'));
seg(ok) = single(r(idx(ok)));
t0 = g0 + (i0 - 1)*Ts + (l0 - 1)*Ts/U;
end
