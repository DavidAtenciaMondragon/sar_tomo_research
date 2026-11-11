function [sig_t, sig_phi, beta] = crb_directo(P, SNR_dB)
% CRB_DIRECTO  Cota de Cramer-Rao del retardo y la fase del camino directo.
%   y = A*s(t - a)*exp(j*phi) + ruido blanco complejo; con SNR = |A|^2*E_s/sigma^2
%   (SNR tras filtro adaptado) y espectro centrado (retardo y fase desacoplados):
%     var(a)   >= 1/(2*SNR*(2*pi*beta)^2),   var(phi) >= 1/(2*SNR)
%   beta = ancho de banda rms del pulso realmente emitido (con flancos).
%   SNR_dB se refiere a A^2*Np (convencion de P.SNR_d_dB); se corrige por la
%   energia real del pulso con flancos.
Ts = 1/P.fs; Np = round(P.Tp*P.fs);
p = chirp_bb((0:Np-1).'*Ts, P.Tp, P.B/P.Tp);
Nf = 8192; Pf = abs(fft(p, Nf)).^2;
f = [0:Nf/2-1, -Nf/2:-1].'/(Nf*Ts);
beta = sqrt(sum(f.^2.*Pf)/sum(Pf));
snr = 10.^(SNR_dB/10)*sum(abs(p).^2)/Np;
sig_t = 1./(2*pi*beta*sqrt(2*snr));
sig_phi = 1./sqrt(2*snr);
end
