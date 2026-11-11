function p = chirp_bb(t, Tp, K)
% LFM centrado en banda base con flancos coseno alzado (10 % de Tp): emula
% el filtrado anti-alias y evita sesgos de retardo fraccional por los flancos
r = 0.1*Tp;
env = ones(size(t));
i1 = t < r;        env(i1) = 0.5 - 0.5*cos(pi*t(i1)/r);
i2 = t > Tp - r;   env(i2) = 0.5 - 0.5*cos(pi*(Tp - t(i2))/r);
p = env.*exp(1j*pi*K*(t - Tp/2).^2).*(t >= 0 & t < Tp);
end
