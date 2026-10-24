"""Experimento D2 - Separabilidad horizontal: rejilla 4x4 a z = -2 m con separacion variable.

La rejilla de 64 puntos del Exp. D (1 m en xy, 1.5 m en z) es resuelta por todos los canales
salvo Ref1 y no discrimina; una version densa en 4 profundidades queda dominada por la
interferencia vertical (FWHM_z = 1.39 m, igual en todos los canales).  Aqui se aisla la
separabilidad horizontal: 16 puntos en un plano, separacion s = 0.15..0.30 m (FWHM mono =
0.157 m), fases aleatorias fijas por punto (3 realizaciones), suelo S0.
Salidas: expD/expD2_separabilidad.csv, figura D2_separabilidad.
Uso: python exp_D2_rejilla_densa.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bisar import config, geometry, imaging, metrics, plotting as P, runner, scene

p = config.load("expD")
out = config.outdir("expD")
t, Ptx = geometry.tx_trajectory(p)
sys0 = scene.System(p, soil=p["suelos"]["S0"]["eps_r"])
g = imaging.Grid(np.arange(-0.65, 0.6501, 0.01), np.arange(-0.65, 0.6501, 0.01), np.arange(-2.4, -1.59, 0.1))
rng = np.random.default_rng(p["semilla"] + 64)
rows, keep = [], {}
for sp in [0.15, 0.20, 0.25, 0.30]:
    c = (np.arange(4) - 1.5) * sp
    pts = np.array([[x, y, -2.0] for x in c for y in c])
    for real in range(3):
        ph = np.exp(2j * np.pi * rng.random(len(pts)))
        tg = [(q, None, w) for q, w in zip(pts, ph)]
        s_m = scene.multi_target_echo(sys0, Ptx, Ptx, tg)
        I_m = runner.bp(sys0, Ptx, Ptx, s_m, g.vox)
        res = {("mono", 0): I_m}
        for dphi in [45, 90, 120, 135]:
            Prx = geometry.rx_from_tx(Ptx, dphi)
            s_b = scene.multi_target_echo(sys0, Ptx, Prx, tg)
            I_b = runner.bp(sys0, Ptx, Prx, s_b, g.vox)
            res[("bi", dphi)] = I_b
            res[("fus", dphi)] = 0.5 * (I_m + I_b)
        for (ch, dphi), I in res.items():
            s = metrics.separability(g, I, pts, tol_h=sp / 3, tol_z=0.25)
            rows.append({"sep_m": sp, "real": real, "canal": ch, "dphi": dphi,
                         "frac_resueltos": s["frac_resolved"], "loc_err_medio_m": s["loc_err_mean"],
                         "disp_amplitud_dB": s["amp_spread_dB"]})
        if real == 0:
            keep[sp] = res
    print(f"sep {sp}", flush=True)
d = pd.DataFrame(rows)
d.to_csv(f"{out}/expD2_separabilidad.csv", index=False)
sm = d.groupby(["canal", "dphi", "sep_m"]).frac_resueltos.mean().unstack("sep_m")
sm.to_csv(f"{out}/expD2_separabilidad_resumen.csv")
print(sm.round(3).to_string())

fig, axs = plt.subplots(1, 5, figsize=(14, 3.3))
fig.subplots_adjust(wspace=0.4)
a = axs[0]
seps = sorted(d.sep_m.unique())
for (ch, dphi), col in zip([("mono", 0), ("fus", 45), ("fus", 90), ("fus", 120), ("fus", 135)],
                           [P.C["mono"], P.SERIES[6], P.C["fus"], P.SERIES[3], P.SERIES[4]]):
    a.plot(seps, [sm.loc[(ch, dphi), s] * 100 for s in seps], marker="o", color=col,
           label="Mono" if ch == "mono" else f"Mono+Bi {dphi}°")
a.set_xlabel("separacion [m]"); a.set_ylabel("blancos resueltos [%]")
a.legend(fontsize=6.5)
iz = np.argmin(np.abs(g.z + 2.0))
sp = 0.20
ref = np.abs(keep[sp][("mono", 0)]).max()
c = (np.arange(4) - 1.5) * sp
for a, key, ti in zip(axs[1:], [("mono", 0), ("bi", 90), ("fus", 90), ("fus", 120)],
                      ["Mono", "Bi 90°", "Mono+Bi 90°", "Mono+Bi 120°"]):
    im = P.show_slice(a, g.x, g.y, np.abs(g.cube(keep[sp][key]))[:, :, iz], ref, ti + f" (s = {sp} m)")
    X, Y = np.meshgrid(c, c)
    a.scatter(X, Y, marker="+", color=P.C["s8"], s=30, lw=1)
fig.colorbar(im, ax=axs[1:], shrink=0.8, label="dB re pico mono")
P.save(fig, "D2_separabilidad", out)
