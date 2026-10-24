"""Figura de enlace del Exp. C (Fig. 7): reflectancia TM/TE de los suelos y robustez del minimo
frente a error angular y de atitud para el suelo sin perdidas (S0) y el de referencia (S1).
Solo usa formulas analiticas (bisar.em); rapido.  Uso: python exp_C_figuras_enlace.py"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bisar import config, em, plotting as P

p = config.load("expC")
out = config.outdir("expC")
eps = {k: em.eps_complex(v["eps_r"]) for k, v in p["suelos"].items()}
thB = {k: em.pseudo_brewster_deg(e)[0] for k, e in eps.items()}
th = np.linspace(0, 89.5, 1791)
fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.8))
fig.subplots_adjust(wspace=0.3)
for (k, e), col in zip(eps.items(), P.SERIES[:4]):
    axs[0].plot(th, 10 * np.log10(em.reflectance(np.radians(th), e, "TM")), color=col,
                label=f"{k}: $\\varepsilon_r$={e.real:g}-j{-e.imag:g}, $\\theta_{{pB}}$={thB[k]:.2f}°, "
                      f"min {10*np.log10(em.pseudo_brewster_deg(e)[1] + 1e-30):.1f} dB")
    axs[0].plot(th, 10 * np.log10(em.reflectance(np.radians(th), e, "TE")), color=col, ls=":")
axs[0].axvspan(46.27, 60.52, color=P.C["grid"], alpha=0.5, lw=0)
axs[0].text(53.4, -2.5, "helice base", ha="center", fontsize=7, color=P.C["ink2"])
axs[0].set_ylim(-50, 0); axs[0].set_xlim(0, 90)
axs[0].set_xlabel("incidencia $\\theta_i$ [°]"); axs[0].set_ylabel("reflectancia [dB]")
axs[0].set_title("TM (continuo) y TE (punteado)")
axs[0].legend(fontsize=6, loc="lower left")
d = np.linspace(-6, 6, 481)
rows = []
for a, k in zip(axs[1:], ["S0", "S1"]):
    for r, col in zip(p["roll_list_deg"], P.SERIES[:4]):
        G = 10 * np.log10(em.reflectance_mixed(np.radians(thB[k] + d), eps[k], r) + 1e-30)
        a.plot(d, G, color=col, label=f"roll {r}°")
        for e_ in [0, 0.5, 1, 2, 5]:
            for sgn in [-1, 1]:
                rows.append({"suelo": k, "roll_deg": r, "err_deg": sgn * e_,
                             "Gamma_TM_ef_dB": float(10 * np.log10(em.reflectance_mixed(
                                 np.radians(thB[k] + sgn * e_), eps[k], r) + 1e-30))})
    for e_ in [0.5, 1, 2, 5]:
        a.axvline(e_, color=P.C["grid"], lw=0.7); a.axvline(-e_, color=P.C["grid"], lw=0.7)
    a.set_ylim(-70, -10)
    a.set_xlabel("error respecto a $\\theta_{pB}$ [°]"); a.set_ylabel("reflectancia TM efectiva [dB]")
    a.set_title(f"{k}: error de incidencia y de atitud (roll)")
axs[1].legend(fontsize=7)
P.save(fig, "C_reflectancia", out)
pd.DataFrame(rows).drop_duplicates().to_csv(f"{out}/expC_enlace_error_roll.csv", index=False)
print(pd.DataFrame(rows).drop_duplicates().pivot_table(index=["suelo", "err_deg"], columns="roll_deg",
                                                      values="Gamma_TM_ef_dB").round(1).to_string())
