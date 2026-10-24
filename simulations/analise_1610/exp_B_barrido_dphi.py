"""Experimento B - Barrido de la separacion acimutal Tx-Rx (dphi = 0:15:180).

Para cada blanco (centro y descentrado) y cada dphi: imagenes mono, bi, fus (2N) y fusN (N),
metricas de PSF, soporte en k, condicionamiento.  CSV: io/analise_1610/expB/expB_metricas.csv
Uso: python exp_B_barrido_dphi.py
"""
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bisar import config, geometry, imaging, metrics, plotting as P, runner, scene

NAME = "expB"
p = config.load(NAME)
config.save_resolved(p, NAME)
out = config.outdir(NAME)
sysm = scene.System(p)
t, Ptx = geometry.tx_trajectory(p)
M = len(Ptx)
rg = p["rejilla_psf"]
odd = (np.arange(M) % 2 == 1).astype(float)

rows, krows, saved = [], [], {}
for tname, tpos in p["blancos"].items():
    tgt = np.array(tpos)
    g = imaging.Grid.around(tgt, rg["hx"], rg["dx"], rg["z_lo"], rg["z_hi"], rg["dz"])
    s_m, inf_m = scene.target_echo(sysm, Ptx, Ptx, tgt)
    I_m = runner.bp(sysm, Ptx, Ptx, s_m, g.vox)
    I_mo = runner.bp(sysm, Ptx, Ptx, s_m, g.vox, odd)
    Km = metrics.kvectors(sysm, inf_m["gt"], inf_m["gr"])
    for dphi in p["dphi_list"]:
        t0 = time.time()
        _, Prx = runner.channels(p, dphi, Ptx=Ptx)
        s_b, inf_b = scene.target_echo(sysm, Ptx, Prx, tgt)
        imgs = runner.four_images(sysm, Ptx, Prx, s_m, s_b, g.vox, I_m, I_mo)
        am, ab = np.mean(np.abs(s_m)), np.mean(np.abs(s_b))
        ref = {"mono": am, "bi": ab, "fus": 0.5 * (am + ab), "fusN": 0.5 * (am + ab)}
        r = runner.psf_rows(sysm, g, imgs, tgt, M, {"blanco": tname, "dphi": dphi}, ref)
        rho = abs(np.vdot(imgs["fus"], imgs["fusN"])) / (np.linalg.norm(imgs["fus"]) * np.linalg.norm(imgs["fusN"]))
        sep = np.linalg.norm(Ptx - Prx, axis=1)
        for row in r:
            row["corr_fus_fusN"] = rho
            row["sep_min_m"] = sep.min()
            row["Pr_medio_dBm"] = 10 * np.log10(np.mean(np.abs(s_b if row["canal"] == "bi" else s_m) ** 2)) + 30
            if row["canal"] == "bi" and (row["fwhm_x_cens"] or row["fwhm_y_cens"]):
                row.update(runner.cut_widths(sysm, Ptx, Prx, s_b, tgt, h=25.0, dx=0.01))
        rows += r
        # informacion: soporte en k y condicionamiento (Gram normalizado por canal)
        Kb = metrics.kvectors(sysm, inf_b["gt"], inf_b["gr"])
        for ch, Ks, ws in [("mono", [Km], [1.0]), ("bi", [Kb], [1.0]), ("fus", [Km, Kb], [0.5, 0.5])]:
            K = np.vstack(Ks)
            km = metrics.kspace_metrics(K)
            cm = metrics.conditioning(Ks, ws)
            krows.append({"blanco": tname, "dphi": dphi, "canal": ch, **km, **cm})
        if tname == "centro":
            saved[dphi] = {"imgs": {k: v for k, v in imgs.items() if not k.startswith("_")}, "Kb": Kb}
        print(f"{tname:11s} dphi={dphi:5.1f}  {time.time()-t0:5.1f}s  "
              f"fwhm_x bi={r[1]['fwhm_x']:.3f} fus={r[2]['fwhm_x']:.3f} | PSLR fus={r[2]['PSLR_dB']:.2f} "
              f"| ISLR fus={r[2]['ISLR_dB']:.2f} | corr fus/fusN={rho:.5f}", flush=True)
    if tname == "centro":
        saved["grid"] = g
        saved["Km"] = Km
        saved["I_m"] = I_m

df = pd.DataFrame(rows)
df.to_csv(f"{out}/expB_metricas.csv", index=False)
dk = pd.DataFrame(krows)
dk.to_csv(f"{out}/expB_kspace_condicionamiento.csv", index=False)

# ============================================================ figuras
D = df[df.blanco == "centro"]
dph = np.array(p["dphi_list"], float)
mono = D[D.canal == "mono"].iloc[0]


def series(col, ch, data=D):
    return data[data.canal == ch].sort_values("dphi")[col].values


# --- Fig. 4 resoluciones
fig, axs = plt.subplots(1, 3, figsize=(11, 3.3))
for a, lab in zip(axs, ["x", "y", "z"]):
    for ch, ls, mk in [("bi", "-", "o"), ("fus", "-", "s"), ("fusN", ":", None)]:
        y = series(f"fwhm_{lab}", ch)
        if lab != "z" and ch == "bi":
            y = np.where(series(f"fwhm_{lab}_cens", ch) == 1, np.nan, y)
        a.plot(dph, y, color=P.C[ch], ls=ls, marker=mk, label=P.LAB[ch] if ch != "bi" else "Bi")
    a.axhline(mono[f"fwhm_{lab}"], color=P.C["mono"], lw=2, label="Mono (Ref2)")
    if lab != "z":
        th = mono[f"fwhm_{lab}"] / np.abs(np.cos(np.radians(dph) / 2))
        a.plot(dph, np.where(dph < 170, th, np.nan), color=P.C["ink2"], lw=0.9, ls="--",
               label="teoria $\\delta_{xy}^{mono}/|\\cos(\\Delta\\phi/2)|$")
        a.set_ylim(0, 0.8)
    a.set_xlabel("$\\Delta\\phi$ [°]")
    a.set_ylabel(f"$\\delta_{lab}$ (-3 dB) [m]")
    a.set_xticks(range(0, 181, 30))
axs[0].legend(fontsize=7, loc="upper left")
axs[2].set_ylim(1.2, 1.6)
fig.suptitle("Exp. B: resolucion -3 dB vs separacion acimutal (blanco en (0,0,-2) m; Bi en 180° censurado: >17 m)", y=1.03)
P.save(fig, "B_resoluciones", out)

# --- Fig. 5 volumen y anisotropia
fig, axs = plt.subplots(1, 2, figsize=(9, 3.3))
for ch, ls, mk in [("bi", "-", "o"), ("fus", "-", "s"), ("fusN", ":", None)]:
    v = series("vol3dB_m3", ch)
    v = np.where(series("vol3dB_cens", ch) == 1, np.nan, v)
    axs[0].plot(dph, v * 1e3, color=P.C[ch], ls=ls, marker=mk, label=P.LAB[ch] if ch != "bi" else "Bi")
    axs[1].plot(dph, series("aniso3d", ch), color=P.C[ch], ls=ls, marker=mk)
axs[0].axhline(mono.vol3dB_m3 * 1e3, color=P.C["mono"], label="Mono (Ref2)")
axs[1].axhline(mono.aniso3d, color=P.C["mono"])
axs[0].set_ylabel("volumen del lobulo -3 dB [litros]")
axs[1].set_ylabel("anisotropia $\\sqrt{\\lambda_{max}/\\lambda_{min}}$")
for a in axs:
    a.set_xlabel("$\\Delta\\phi$ [°]"); a.set_xticks(range(0, 181, 30))
axs[0].legend(fontsize=7)
axs[0].set_ylim(0, 160)
fig.suptitle("Exp. B: volumen y anisotropia de la PSF 3D", y=1.02)
P.save(fig, "B_volumen_anisotropia", out)

# --- Fig. 6 PSLR / ISLR
fig, axs = plt.subplots(1, 3, figsize=(11, 3.3))
for ch, ls, mk in [("bi", "-", "o"), ("fus", "-", "s"), ("fusN", ":", None)]:
    axs[0].plot(dph, series("PSLR_dB", ch), color=P.C[ch], ls=ls, marker=mk, label=P.LAB[ch] if ch != "bi" else "Bi")
    axs[1].plot(dph, series("ISLR_dB", ch), color=P.C[ch], ls=ls, marker=mk)
    axs[2].plot(dph, series("contrast_dB", ch), color=P.C[ch], ls=ls, marker=mk)
for a, col in zip(axs, ["PSLR_dB", "ISLR_dB", "contrast_dB"]):
    a.axhline(mono[col], color=P.C["mono"], label="Mono (Ref2)")
    a.set_xlabel("$\\Delta\\phi$ [°]"); a.set_xticks(range(0, 181, 30))
axs[0].set_ylabel("PSLR 3D [dB]"); axs[1].set_ylabel("ISLR 3D [dB] (ventana fija)")
axs[2].set_ylabel("contraste pico/fondo [dB]")
axs[0].legend(fontsize=7)
fig.suptitle("Exp. B: lobulos secundarios y contraste (ventana 2 x 2 x 5 m alrededor del blanco)", y=1.02)
P.save(fig, "B_pslr_islr", out)

# --- Fig. 2 soporte en k
g = saved["grid"]
Km = saved["Km"]
fig, axs = plt.subplots(2, len(p["kspace_dphi"]), figsize=(13, 5.4))
sub = slice(None, None, 41)
for j, dphi in enumerate(p["kspace_dphi"]):
    Kb = saved[dphi]["Kb"]
    for K, ch, lab in [(Km, "mono", "mono"), (Kb, "bi", "bi")]:
        axs[0, j].scatter(K[sub, 0], K[sub, 1], s=1.5, color=P.C[ch], label=lab, rasterized=True)
        axs[1, j].scatter(np.hypot(K[sub, 0], K[sub, 1]), K[sub, 2], s=1.5, color=P.C[ch], rasterized=True)
    axs[0, j].set_aspect("equal"); axs[0, j].set_xlim(-17, 17); axs[0, j].set_ylim(-17, 17)
    axs[0, j].set_title(("Ref2 (mono)" if dphi == 0 else ("Ref1" if dphi == 180 else "")) + f"  $\\Delta\\phi$={dphi}°")
    axs[0, j].set_xlabel("$k_x$ [rad/m]")
    axs[1, j].set_xlabel("$|k_h|$ [rad/m]")
    axs[1, j].set_xlim(-0.5, 17)
axs[0, 0].set_ylabel("$k_y$ [rad/m]"); axs[1, 0].set_ylabel("$k_z$ [rad/m]")
axs[0, 0].legend(markerscale=5, loc="lower left", fontsize=7)
fig.subplots_adjust(hspace=0.35, wspace=0.35)
fig.suptitle("Exp. B: soporte en el dominio de k (azul: canal mono, naranja: canal bi); la fusion es la union", y=0.98)
P.save(fig, "B_kspace", out)

# --- informacion vs dphi
K = dk[dk.blanco == "centro"]
fig, axs = plt.subplots(1, 3, figsize=(11, 3.3))
for ch, mk in [("mono", None), ("bi", "o"), ("fus", "s")]:
    d = K[K.canal == ch].sort_values("dphi")
    axs[0].plot(d.dphi, d.k_occ_area_h, color=P.C[ch], marker=mk, label=P.LAB.get(ch, ch) if ch != "bi" else "Bi")
    axs[1].plot(d.dphi, d.erank, color=P.C[ch], marker=mk)
    axs[2].plot(d.dphi, d.coh_max, color=P.C[ch], marker=mk)
axs[0].set_ylabel("area ocupada en $(k_x,k_y)$ [rad$^2$/m$^2$]")
axs[1].set_ylabel("rango efectivo de A (de 245)")
axs[2].set_ylabel("coherencia mutua max.")
for a in axs:
    a.set_xlabel("$\\Delta\\phi$ [°]"); a.set_xticks(range(0, 181, 30))
axs[0].legend(fontsize=7)
fig.suptitle("Exp. B: cobertura del dominio de k y condicionamiento de la matriz de adquisicion", y=1.02)
P.save(fig, "B_informacion", out)

# --- Fig. 9 imagenes lado a lado (misma escala dinamica, misma referencia = pico mono)
sel = p["imagenes_dphi"]
iz = np.argmin(np.abs(g.z + 2.0)); iy = np.argmin(np.abs(g.y))
ref = np.abs(saved["I_m"]).max()
fig, axs = plt.subplots(3, len(sel) + 1, figsize=(13, 7.6))
fig.subplots_adjust(wspace=0.45, hspace=0.45)
A = np.abs(g.cube(saved["I_m"]))
P.show_slice(axs[0, 0], g.x, g.y, A[:, :, iz], ref, "Mono (Ref2)\nz = -2 m")
P.show_slice(axs[1, 0], g.x, g.z, A[:, iy, :], ref, "Mono, y = 0", "x [m]", "z [m]")
axs[2, 0].axis("off")
for j, dphi in enumerate(sel):
    for i, ch in enumerate(["bi", "fus"]):
        A = np.abs(g.cube(saved[dphi]["imgs"][ch]))
        im = P.show_slice(axs[i, j + 1], g.x, g.y, A[:, :, iz], ref,
                          f"{'Bi' if ch == 'bi' else 'Mono+Bi'} $\\Delta\\phi$={dphi}°\nz = -2 m")
    A = np.abs(g.cube(saved[dphi]["imgs"]["fus"]))
    P.show_slice(axs[2, j + 1], g.x, g.z, A[:, iy, :], ref, f"Mono+Bi {dphi}°, y = 0", "x [m]", "z [m]")
fig.colorbar(im, ax=axs, shrink=0.5, label="dB respecto al pico mono (rango 30 dB)")
P.save(fig, "B_imagenes", out)
print("Exp. B terminado ->", out)
