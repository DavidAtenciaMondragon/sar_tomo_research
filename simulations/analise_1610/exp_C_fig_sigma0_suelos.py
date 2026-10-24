"""Exp. C - Mapas de sigma0 biestatico SPM (VV y HH) para los cuatro suelos S0..S3.

Rx a la misma elevacion que el Tx (theta_s = theta_i = theta) y separacion acimutal dphi
(phi_s = 180 - dphi: 0 = especular hacia adelante = Ref1, 180 = retrodispersion = mono).
Marca: banda de incidencias de la helice base, theta_pB de cada suelo (linea discontinua),
lugar del minimo de VV (nulo exacto si eps es real) y nulo de HH (dphi = 90).
Uso: python exp_C_fig_sigma0_suelos.py
"""
import numpy as np
import matplotlib.pyplot as plt

from bisar import config, em, geometry, plotting as P

NAME = "expC"
p = config.load(NAME)
out = config.outdir(NAME)
_, P0 = geometry.tx_trajectory(p)
th0 = geometry.incidence_deg(P0)
k0 = 2 * np.pi * p["fc"] / p["c"]
ru = p["rugosidad"]
TH, DP = np.meshgrid(np.linspace(20, 80, 241), np.linspace(0, 180, 361), indexing="ij")
phis = np.radians(180 - DP)

fig, axs = plt.subplots(2, 4, figsize=(15, 6.6), sharex=True, sharey=True)
for j, sk in enumerate(["S0", "S1", "S2", "S3"]):
    e = em.eps_complex(p["suelos"][sk]["eps_r"])
    thB = em.pseudo_brewster_deg(e)[0]
    for i, pol in enumerate(["TM", "TE"]):
        a = axs[i, j]
        s0 = em.spm_sigma0(np.radians(TH), np.radians(TH), phis, e, pol, k0, ru["sigma_h"], ru["corr_len"])
        im = a.pcolormesh(DP, TH, 10 * np.log10(np.maximum(s0, 1e-12)), cmap=P.SEQ, vmin=-60, vmax=-10,
                          shading="auto", rasterized=True)
        a.grid(False)
        a.axhspan(th0.min(), th0.max(), color="w", alpha=0.18, lw=0)
        a.axhline(thB, color=P.C["s8"], lw=0.9, ls="--")
        if pol == "TM":
            # minimo de VV en cada incidencia (nulo exacto si eps real): busqueda numerica en dphi
            imin = np.argmin(s0, axis=1)
            depth = 10 * np.log10(s0[np.arange(len(imin)), imin] / s0[:, -1])
            ok = (imin > 0) & (imin < DP.shape[1] - 1)
            a.plot(DP[0, imin[ok]], TH[ok, 0], color=P.C["s8"], lw=1.3, label="mínimo de VV")
            a.set_title(f"{p['suelos'][sk]['nombre']}\n$\\varepsilon_r$={e.real:g}$-j${-e.imag:g},  "
                        f"$\\theta_{{pB}}$={thB:.1f}°   |   VV", fontsize=8.5)
        else:
            a.axvline(90, color=P.C["s8"], lw=1.3, label="nulo de HH")
            a.set_title("HH", fontsize=8.5)
        if j == 0:
            a.set_ylabel(f"incidencia $\\theta$ [°]")
            a.legend(fontsize=7, loc="lower left")
        if i == 1:
            a.set_xlabel("$\\Delta\\phi$ Tx-Rx [°]")
        a.set_xticks(range(0, 181, 45))
        # anotaciones mono / Ref1
        if i == 0 and j == 0:
            a.text(3, 77, "mono", fontsize=7, color="w"); a.text(150, 77, "Ref1", fontsize=7, color="w")
fig.colorbar(im, ax=axs, shrink=0.85, label="$\\sigma^0$ SPM [dB]")
fig.suptitle(f"$\\sigma^0$ biestático SPM (1er orden), Rx a la misma elevación que el Tx, $s$={ru['sigma_h']*100:g} cm, "
             f"$\\ell$={ru['corr_len']*100:g} cm. Banda clara: hélice base; discontinua: $\\theta_{{pB}}$ del suelo", y=0.995)
P.save(fig, "C_sigma0_mapa_suelos", out)

# profundidad del minimo de VV en la incidencia base y en theta_pB (para el texto)
for sk in ["S0", "S1", "S2", "S3"]:
    e = em.eps_complex(p["suelos"][sk]["eps_r"])
    for th in [th0.mean(), em.pseudo_brewster_deg(e)[0]]:
        dp = np.linspace(0, 180, 3601)
        s = em.spm_sigma0(np.radians(th), np.radians(th), np.radians(180 - dp), e, "TM", k0, ru["sigma_h"], ru["corr_len"])
        k = np.argmin(s)
        print(f"{sk} theta={th:5.2f}: minimo VV en dphi={dp[k]:6.1f}  "
              f"profundidad respecto de la mono {10*np.log10(s[k]/s[0]):6.1f} dB, "
              f"ancho a -10 dB del minimo: {np.ptp(dp[10*np.log10(s/s[k]) < 10]):.1f} deg")
