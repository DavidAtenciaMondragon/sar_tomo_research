"""Exp. B - Figura de geometria: trayectorias reales del Tx y del Rx para cada dphi.

Aclara que en el Exp. B el Rx recorre la MISMA helice que el Tx (mismos radios y alturas),
rotada rigidamente dphi en acimut y sincronizada; lo unico que cambia entre columnas es
el desfase acimutal inicial.  Fila 1: planta (x, y); fila 2: rho(t) y z(t) de ambos.
Uso: python exp_B_fig_geometria.py
"""
import numpy as np
import matplotlib.pyplot as plt

from bisar import config, geometry, plotting as P, runner

NAME = "expB"
p = config.load(NAME)
out = config.outdir(NAME)
t, Ptx = geometry.tx_trajectory(p)
sel = p["kspace_dphi"]
k_inst = [0, 400, 800, 1200]          # instantes marcados (pulsos)

fig, axs = plt.subplots(2, len(sel), figsize=(13, 5.8), gridspec_kw={"height_ratios": [1.3, 1]})
for j, dphi in enumerate(sel):
    _, Prx = runner.channels(p, dphi, Ptx=Ptx)
    a = axs[0, j]
    a.plot(Ptx[:, 0], Ptx[:, 1], color=P.C["mono"], lw=1.8, label="Tx (transceiver)")
    a.plot(Prx[:, 0], Prx[:, 1], color=P.C["bi"], lw=1.1, ls="--", label="Rx remoto")
    for k in k_inst:
        a.plot([Ptx[k, 0], Prx[k, 0]], [Ptx[k, 1], Prx[k, 1]], color=P.C["ink2"], lw=0.7, ls=":")
        a.plot(*Ptx[k, :2], "o", color=P.C["mono"], ms=5)
        a.plot(*Prx[k, :2], "s", color=P.C["bi"], ms=5)
    a.plot(*Ptx[0, :2], "o", mfc="none", mec=P.C["ink"], ms=9)
    a.plot(*Prx[0, :2], "s", mfc="none", mec=P.C["ink"], ms=9)
    a.plot(0, 0, "+", color=P.C["ink"], ms=8)
    a.set_aspect("equal"); a.set_xlim(-185, 185); a.set_ylim(-185, 185)
    a.set_xlabel("x [m]")
    tag = "Ref2 (mono)  " if dphi == 0 else ("Ref1  " if dphi == 180 else "")
    a.set_title(tag + r"$\Delta\phi$=" + f"{dphi}°")
    rt, rr = np.hypot(*Ptx[:, :2].T), np.hypot(*Prx[:, :2].T)
    b = axs[1, j]
    b.plot(t, rt, color=P.C["mono"], lw=2.6, label=r"Tx")
    b.plot(t, rr, color=P.C["bi"], lw=1.1, ls="--", label=r"Rx")
    b.plot(t, Ptx[:, 2], color=P.C["mono"], lw=2.6)
    b.plot(t, Prx[:, 2], color=P.C["bi"], lw=1.1, ls="--")
    b.text(7.0, 166, r"radio $\rho(t)$", fontsize=7.5, color=P.C["ink2"])
    b.text(0.4, 100, r"altura $z(t)$", fontsize=7.5, color=P.C["ink2"])
    b.set_ylim(85, 180); b.set_xlabel("t [s]")
    b.set_title(r"max$|\rho_{Rx}-\rho_{Tx}|$=" + f"{np.abs(rr - rt).max():.0e} m\n"
                r"max$|z_{Rx}-z_{Tx}|$=" + f"{np.abs(Prx[:, 2] - Ptx[:, 2]).max():.0f} m", fontsize=7.5)
axs[0, 0].set_ylabel("y [m]"); axs[1, 0].set_ylabel("[m]")
axs[0, 0].legend(fontsize=7, loc="lower left"); axs[1, 0].legend(fontsize=7, loc="center right")
fig.subplots_adjust(hspace=0.5, wspace=0.3)
fig.suptitle(r"Exp. B: el Rx recorre la misma hélice que el Tx (radio 137.4$\to$165.5 m, altura 131.5$\to$93.6 m), "
             r"desfasada $\Delta\phi$ en acimut" + "\nmarcas = mismos instantes (pulsos 0, 400, 800, 1200; "
             "contorno = inicio); línea punteada = línea de base Tx–Rx", y=1.0)
P.save(fig, "B_geometria", out)
print("ok")
