"""Experimento A - Referencias Ref2 (mono), Ref1 (diametral) y fusion Ref1+Ref2.

Salidas (io/analise_1610/expA):
  datos_complejos_expA.npz   ecos complejos mono y bi (mismo pulso) + trayectorias
  expA_metricas.csv          metricas de PSF por canal (mono, bi=Ref1, fus, fusN) + Ref intermedia
  figuras: A_geometria, A_psf_cortes, A_imagenes, A_kspace_ref
Uso: python exp_A_referencias.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bisar import config, geometry, imaging, metrics, plotting as P, runner, scene

NAME = "expA"
p = config.load(NAME)
config.save_resolved(p, NAME)
out = config.outdir(NAME)
sysm = scene.System(p)
tgt = np.array(p["blanco"])
t, Ptx = geometry.tx_trajectory(p)
M = len(Ptx)
g = imaging.Grid.around(tgt, p["rejilla_psf"]["hx"], p["rejilla_psf"]["dx"],
                        p["rejilla_psf"]["z_lo"], p["rejilla_psf"]["z_hi"], p["rejilla_psf"]["dz"])

s_m, inf_m = scene.target_echo(sysm, Ptx, Ptx, tgt)
I_m = runner.bp(sysm, Ptx, Ptx, s_m, g.vox)
odd = (np.arange(M) % 2 == 1).astype(float)
I_m_odd = runner.bp(sysm, Ptx, Ptx, s_m, g.vox, odd)

rows, store, cases = [], {"Ptx": Ptx, "t": t, "s_mono": s_m, "freqs": sysm.freqs}, {}
for lab, dphi in [("Ref1", p["dphi_ref1"]), ("Intermedia", p["dphi_intermedia"])]:
    _, Prx = runner.channels(p, dphi, Ptx=Ptx)
    s_b, inf_b = scene.target_echo(sysm, Ptx, Prx, tgt)
    imgs = runner.four_images(sysm, Ptx, Prx, s_m, s_b, g.vox, I_m, I_m_odd)
    ref = {"mono": np.mean(np.abs(s_m)), "bi": np.mean(np.abs(s_b)),
           "fus": 0.5 * (np.mean(np.abs(s_m)) + np.mean(np.abs(s_b))),
           "fusN": 0.5 * (np.mean(np.abs(s_m)) + np.mean(np.abs(s_b)))}
    r = runner.psf_rows(sysm, g, imgs, tgt, M, {"caso": lab, "dphi": dphi}, ref)
    # anchos horizontales en cortes largos (Ref1 queda censurado en la rejilla de +-1 m)
    for row in r:
        ch = row["canal"]
        if ch in ("mono", "bi"):
            row.update(runner.cut_widths(sysm, Ptx, Ptx if ch == "mono" else Prx,
                                         s_m if ch == "mono" else s_b, tgt, h=p["corte_largo_h"], dx=0.01))
    rows += r
    store[f"Prx_{lab}"] = Prx
    store[f"s_bi_{lab}"] = s_b
    cases[lab] = (Prx, s_b, imgs, inf_b)
    # potencia media recibida por canal (solo informativa: NO es criterio)
    print(f"{lab}: Pr medio mono = {10*np.log10(np.mean(np.abs(s_m)**2))+30:.2f} dBm, "
          f"bi = {10*np.log10(np.mean(np.abs(s_b)**2))+30:.2f} dBm")
np.savez_compressed(f"{out}/datos_complejos_expA.npz", **store)
df = pd.DataFrame(rows)
df.to_csv(f"{out}/expA_metricas.csv", index=False)
cols = ["caso", "canal", "fwhm_x", "fwhm_y", "fwhm_z", "fwhm_x_largo", "vol3dB_m3", "vol3dB_cens", "aniso3d",
        "aniso_xy", "PSLR_dB", "ISLR_dB", "contrast_dB", "loc_err_m", "amp_err_dB", "SNR_dB"]
print(df[cols].to_string(index=False, float_format=lambda v: f"{v:.4g}"))

# ------------------------------------------------------------------ figura 1: geometria
Prx1, Prxi = cases["Ref1"][0], cases["Intermedia"][0]
fig = plt.figure(figsize=(11.5, 4.4))
fig.subplots_adjust(wspace=0.32)
ax = fig.add_subplot(1, 2, 1, projection="3d")
ax.plot(*Ptx.T, color=P.C["mono"], lw=1.4, label="Tx transceiver (canal mono = Ref2)")
ax.plot(*Prx1.T, color=P.C["bi"], lw=1.4, label="Rx remoto Ref1 ($\\Delta\\phi$=180°)")
ax.plot(*Prxi.T, color=P.C["s7"], lw=1.2, ls="--", label="Rx intermedio ($\\Delta\\phi$=90°)")
k = 500
for Q, col in [(Prx1[k], P.C["bi"]), (Prxi[k], P.C["s7"]), (Ptx[k], P.C["mono"])]:
    ax.plot([Q[0], tgt[0]], [Q[1], tgt[1]], [Q[2], tgt[2]], color=col, lw=0.8, alpha=0.7)
ax.scatter(*Ptx[k], color=P.C["mono"], s=25)
ax.scatter(*Prx1[k], color=P.C["bi"], s=25)
ax.scatter(*Prxi[k], color=P.C["s7"], s=25)
ax.scatter(*tgt, color="k", s=20)
xx, yy = np.meshgrid([-180, 180], [-180, 180])
ax.plot_surface(xx, yy, 0 * xx, color="#c9b38a", alpha=0.25)
ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_zlabel("z [m]")
ax.set_title("Geometria: helice del Tx y Rx a la misma altura")
ax.legend(loc="upper left", fontsize=7)
ax.view_init(24, -60)
ax2 = fig.add_subplot(1, 2, 2)
rho = np.hypot(Ptx[:, 0], Ptx[:, 1])
ax2.plot(rho, Ptx[:, 2], color=P.C["mono"], label="Tx y Rx: $\\rho(t)$, $z(t)$ (identicos)")
th = geometry.incidence_deg(Ptx)
ax2.set_xlabel("radio horizontal $\\rho$ [m]")
ax2.set_ylabel("altura z [m]")
ax2.set_title(f"Perfil radial: incidencia {th.min():.1f}°–{th.max():.1f}° (media {th.mean():.1f}°)")
for a in [th.min(), 63.45, th.max()]:
    r = np.linspace(0, 200, 2)
    ax2.plot(r, r / np.tan(np.radians(a)), color=P.C["ink2"] if a != 63.45 else P.C["s8"],
             lw=0.8, ls=":")
ax2.text(60, 18, "pseudo-Brewster S1 (63.45°)", color=P.C["s8"], fontsize=7, ha="left")
ax2.set_xlim(0, 200); ax2.set_ylim(0, 150)
ax2.legend(loc="upper left")
P.save(fig, "A_geometria", out)

# ------------------------------------------------------------- figura 3: cortes de PSF
I_b = cases["Ref1"][2]
fig, axs = plt.subplots(1, 3, figsize=(10, 3.0))
ix = np.argmin(np.abs(g.x - tgt[0])); iy = np.argmin(np.abs(g.y - tgt[1])); iz = np.argmin(np.abs(g.z - tgt[2]))
for ch, ls in [("mono", "-"), ("bi", "-"), ("fus", "-"), ("fusN", ":")]:
    A = np.abs(g.cube(I_b[ch]))
    ref = np.abs(g.cube(I_b["mono"])).max()
    for a, prof, axis in [(axs[0], A[:, iy, iz], g.x), (axs[1], A[ix, :, iz], g.y), (axs[2], A[ix, iy, :], g.z)]:
        a.plot(axis, 20 * np.log10(prof / ref + 1e-12), color=P.C[ch], ls=ls,
               lw=2 if ls == "-" else 1.6, label=P.LAB[ch] if ch != "bi" else "Bi Ref1 (180°)")
for a, l in zip(axs, ["x [m]", "y [m]", "z [m]"]):
    a.set_xlabel(l); a.set_ylim(-30, 2); a.axhline(-3, color=P.C["ink2"], lw=0.6, ls="--")
axs[0].set_ylabel("|I| / pico mono [dB]")
axs[0].legend(loc="lower left", fontsize=7)
fig.suptitle("Exp. A: cortes de la PSF por el blanco (escala comun: pico monoestatico = 0 dB)", y=1.03)
P.save(fig, "A_psf_cortes", out)

# --------------------------------------------------------------- figura: imagenes
fig, axs = plt.subplots(2, 4, figsize=(12.5, 6.0))
fig.subplots_adjust(wspace=0.55, hspace=0.5)
ref = np.abs(I_m).max()
for j, (ch, title) in enumerate([("mono", "Ref2 mono"), ("bi", "Ref1 bi (180°)"), ("fus", "Ref1+Ref2 (2N)"),
                                 ("fusN", "Ref1+Ref2 (N)")]):
    A = np.abs(g.cube(I_b[ch]))
    im = P.show_slice(axs[0, j], g.x, g.y, A[:, :, iz], ref, title + "\ncorte z = -2 m")
    P.show_slice(axs[1, j], g.x, g.z, A[:, iy, :], ref, "corte y = 0", "x [m]", "z [m]")
fig.colorbar(im, ax=axs, shrink=0.6, label="dB (pico mono = 0 dB, rango 30 dB)")
P.save(fig, "A_imagenes", out)

# -------------------------------------------- figura: soporte en k de Ref2 vs Ref1
fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
Km = metrics.kvectors(sysm, inf_m["gt"], inf_m["gr"])
Kb = metrics.kvectors(sysm, cases["Ref1"][3]["gt"], cases["Ref1"][3]["gr"])
sub = slice(None, None, 37)
for K, ch, lab in [(Km, "mono", "Ref2 mono"), (Kb, "bi", "Ref1 (180°)")]:
    axs[0].scatter(K[sub, 0], K[sub, 1], s=2, color=P.C[ch], label=lab)
    axs[1].scatter(np.hypot(K[sub, 0], K[sub, 1]), K[sub, 2], s=2, color=P.C[ch], label=lab)
axs[0].set_xlabel("$k_x$ [rad/m]"); axs[0].set_ylabel("$k_y$ [rad/m]"); axs[0].set_aspect("equal")
axs[1].set_xlabel("$|k_h|$ [rad/m]"); axs[1].set_ylabel("$k_z$ [rad/m]")
axs[0].legend(markerscale=4); axs[0].set_title("Plano horizontal: Ref1 colapsa en el origen")
axs[1].set_title("Plano radial-vertical")
P.save(fig, "A_kspace_ref", out)
print("Exp. A terminado ->", out)
