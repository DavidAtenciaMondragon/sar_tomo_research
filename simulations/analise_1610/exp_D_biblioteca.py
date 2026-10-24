"""Experimento D - Biblioteca de blancos P0/P1/P2: mono vs bi vs fusion por clase de blanco.

Salidas (io/analise_1610/expD): expD_objetos.csv, expD_rejilla64.csv, expD_escenaP2.csv
Figuras: D_mip_objetos, D_fusion_vs_mono (Fig. 10), D_rejilla64, D_escenaP2
Uso: python exp_D_biblioteca.py
"""
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bisar import config, geometry, imaging, metrics, plotting as P, rcs, runner, scene

NAME = "expD"
p = config.load(NAME)
config.save_resolved(p, NAME)
out = config.outdir(NAME)
t, Ptx = geometry.tx_trajectory(p)
M = len(Ptx)
lam_s = p["c"] / p["fc"] / scene.System(p).n_real
odd = (np.arange(M) % 2 == 1).astype(float)


def make_obj(spec):
    if spec["tipo"] == "punto":
        return None, None
    V, F = rcs.make_mesh(spec["tipo"], spec["dim"], lam_s)
    R = rcs.rot_xyz(*spec.get("rot", [0, 0, 0]))
    return rcs.PO(V, F, R), R


def region_mask(grid, obj, center, dil=0.3):
    X = grid.vox - center
    if obj is None:
        return np.linalg.norm(X, axis=1) <= dil
    Vc = obj.V
    lo, hi = Vc.min(0) - dil, Vc.max(0) + dil
    return np.all((X >= lo) & (X <= hi), axis=1)


def obj_rows(sysm, g, imgs, center, obj, extra):
    rows = runner.psf_rows(sysm, g, imgs, center, M, extra)
    mask = region_mask(g, obj, center)
    for r in rows:
        I = imgs[r["canal"]]
        E = np.abs(I) ** 2
        r["energia_en_region"] = E[mask].sum() / E.sum()
        r["loc_err_h_m"] = float(np.hypot(r["peak_x"] - center[0], r["peak_y"] - center[1]))
        r["peak_rel_mono_dB"] = 20 * np.log10(r["peak_amp"] / np.abs(imgs["mono"]).max())
    return rows


# =============================================================== objetos P0/P1
sysm = scene.System(p)
tgt = np.array(p["blanco"])
ro = p["rejilla_obj"]
g = imaging.Grid.around(tgt, ro["hx"], ro["dx"], ro["z_lo"], ro["z_hi"], ro["dz"])
rows, keep = [], {}
T0 = time.time()
for name, spec in p["objetos"].items():
    obj, R = make_obj(spec)
    s_m, _ = scene.target_echo(sysm, Ptx, Ptx, tgt, obj)
    I_m = runner.bp(sysm, Ptx, Ptx, s_m, g.vox)
    I_mo = runner.bp(sysm, Ptx, Ptx, s_m, g.vox, odd)
    keep[name] = {"mono": I_m}
    for dphi in p["dphi_list"]:
        Prx = geometry.rx_from_tx(Ptx, dphi)
        s_b, _ = scene.target_echo(sysm, Ptx, Prx, tgt, obj)
        imgs = runner.four_images(sysm, Ptx, Prx, s_m, s_b, g.vox, I_m, I_mo)
        r = obj_rows(sysm, g, imgs, tgt, obj, {"objeto": name, "prioridad": spec["prio"], "dphi": dphi,
                                               "sigma_medio_mono_dBsm": np.nan})
        rows += r
        keep[name][dphi] = {k: imgs[k] for k in ("bi", "fus")}
    print(f"{name:20s} {time.time()-T0:6.0f}s", flush=True)
do = pd.DataFrame(rows)
do.to_csv(f"{out}/expD_objetos.csv", index=False)

# --- figura MIP por objeto (misma escala: pico mono del objeto = 0 dB)
names = list(p["objetos"])
sel = [90, 135]
fig, axs = plt.subplots(len(names), 1 + 2 * len(sel), figsize=(12.5, 2.05 * len(names)))
fig.subplots_adjust(wspace=0.35, hspace=0.95)
for i, nm in enumerate(names):
    ref = np.abs(keep[nm]["mono"]).max()
    panels = [("Mono", keep[nm]["mono"])]
    for d in sel:
        panels += [(f"Bi {d}°", keep[nm][d]["bi"]), (f"Mono+Bi {d}°", keep[nm][d]["fus"])]
    for j, (ti, I) in enumerate(panels):
        A = np.abs(g.cube(I))
        im = P.show_slice(axs[i, j], g.x, g.z, A.max(axis=1), ref, f"{nm}\n{ti}" if j == 0 else ti,
                          "x [m]", "z [m]" if j == 0 else "")
        axs[i, j].tick_params(labelsize=6)
fig.colorbar(im, ax=axs, shrink=0.4, label="MIP en y, dB re pico mono (rango 30 dB)")
P.save(fig, "D_mip_objetos", out)

# --- Fig. 10: metricas de la fusion relativas a la imagen mono, por clase
fig, axs = plt.subplots(1, 4, figsize=(13, 3.8), sharey=True)
mets = [("vol3dB_m3", "volumen -3 dB: fus/mono [dB]", lambda f, m: 10 * np.log10(f / m)),
        ("PSLR_dB", "$\\Delta$PSLR fus-mono [dB]", lambda f, m: f - m),
        ("contrast_dB", "$\\Delta$contraste fus-mono [dB]", lambda f, m: f - m),
        ("energia_en_region", "$\\Delta$energia en region [p.p.]", lambda f, m: 100 * (f - m))]
yy = np.arange(len(names))
for a, (col, lab, fn) in zip(axs, mets):
    for k, (dphi, colr) in enumerate(zip(p["dphi_list"], [P.SERIES[i] for i in (6, 2, 3, 1)])):
        v = []
        for nm in names:
            d = do[(do.objeto == nm) & (do.dphi == dphi)]
            v.append(fn(d[d.canal == "fusN"][col].values[0], d[d.canal == "mono"][col].values[0]))
        a.scatter(v, yy + (k - 1.5) * 0.15, color=colr, s=22, label=f"$\\Delta\\phi$={dphi}°", zorder=3)
    a.axvline(0, color=P.C["ink2"], lw=0.8)
    a.set_xlabel(lab)
axs[0].set_yticks(yy); axs[0].set_yticklabels(names)
axs[0].invert_yaxis()
axs[0].legend(fontsize=7, loc="lower left")
fig.suptitle("Exp. D: fusion (mismo presupuesto N) relativa a la imagen monoestatica, por blanco "
             "(volumen < 0, PSLR < 0: mejor; contraste > 0, energia > 0: mejor)", y=1.02)
P.save(fig, "D_fusion_vs_mono", out)

# =============================================================== rejilla de 64 puntos
r64 = p["rejilla64"]
sys0 = scene.System(p, soil=p["suelos"][r64["suelo"]]["eps_r"])
pts = np.array([[x, y, z] for z in r64["z"] for x in r64["xy"] for y in r64["xy"]])
gm = r64["malla"]
g64 = imaging.Grid(np.arange(-gm["hx"], gm["hx"] + 1e-9, gm["dx"]), np.arange(-gm["hx"], gm["hx"] + 1e-9, gm["dx"]),
                   np.arange(gm["z_lo"], gm["z_hi"] + 1e-9, gm["dz"]))
tg = [(q, None, 1.0) for q in pts]
s_m = scene.multi_target_echo(sys0, Ptx, Ptx, tg)
I_m = runner.bp(sys0, Ptx, Ptx, s_m, g64.vox)
I_mo = runner.bp(sys0, Ptx, Ptx, s_m, g64.vox, odd)
rows64, keep64 = [], {}
for dphi in p["dphi_list"]:
    Prx = geometry.rx_from_tx(Ptx, dphi)
    s_b = scene.multi_target_echo(sys0, Ptx, Prx, tg)
    imgs = runner.four_images(sys0, Ptx, Prx, s_m, s_b, g64.vox, I_m, I_mo)
    keep64[dphi] = imgs
    for ch in ["mono", "bi", "fus", "fusN"]:
        s = metrics.separability(g64, imgs[ch], pts)
        per_depth = {f"resueltos_z{z:g}": float(s["resolved"][pts[:, 2] == z].mean()) for z in r64["z"]}
        rows64.append({"dphi": dphi, "canal": ch, "frac_resueltos": s["frac_resolved"],
                       "loc_err_medio_m": s["loc_err_mean"], "disp_amplitud_dB": s["amp_spread_dB"], **per_depth})
    print(f"rejilla64 dphi={dphi}", flush=True)
d64 = pd.DataFrame(rows64)
d64.to_csv(f"{out}/expD_rejilla64.csv", index=False)
fig, axs = plt.subplots(2, 4, figsize=(12.5, 6.2))
fig.subplots_adjust(wspace=0.4, hspace=0.45)
ref = np.abs(I_m).max()
iy = np.argmin(np.abs(g64.y - 0.5))
for j, (ti, I) in enumerate([("Mono", I_m), ("Bi 90°", keep64[90]["bi"]), ("Mono+Bi 90°", keep64[90]["fus"]),
                             ("Mono+Bi 135°", keep64[135]["fus"])]):
    A = np.abs(g64.cube(I))
    im = P.show_slice(axs[0, j], g64.x, g64.y, A.max(axis=2), ref, ti + "\nMIP en z")
    P.show_slice(axs[1, j], g64.x, g64.z, A[:, iy, :], ref, "corte y = 0.5 m", "x [m]", "z [m]")
    for z in r64["z"]:
        axs[1, j].axhline(z, color=P.C["s8"], lw=0.4, ls=":")
fig.colorbar(im, ax=axs, shrink=0.5, label="dB re pico mono (rango 30 dB)")
P.save(fig, "D_rejilla64", out)

# =============================================================== escena P2
objs = []
for o in p["escena_P2"]:
    ob, _ = make_obj(o)
    objs.append((np.array(o["pos"]), ob, 1.0, o["nombre"]))
rp = p["rejilla_P2"]
gP = imaging.Grid.around(tgt, rp["hx"], rp["dx"], rp["z_lo"], rp["z_hi"], rp["dz"])
tg = [(q, ob, w) for q, ob, w, _ in objs]
s_m = scene.multi_target_echo(sysm, Ptx, Ptx, tg)
I_m = runner.bp(sysm, Ptx, Ptx, s_m, gP.vox)
I_mo = runner.bp(sysm, Ptx, Ptx, s_m, gP.vox, odd)
masks = [region_mask(gP, ob, q) for q, ob, _, _ in objs]
anymask = np.any(masks, axis=0)
rowsP, keepP = [], {}
for dphi in p["dphi_list"]:
    Prx = geometry.rx_from_tx(Ptx, dphi)
    s_b = scene.multi_target_echo(sysm, Ptx, Prx, tg)
    imgs = runner.four_images(sysm, Ptx, Prx, s_m, s_b, gP.vox, I_m, I_mo)
    keepP[dphi] = imgs
    for ch in ["mono", "bi", "fus", "fusN"]:
        A = np.abs(imgs[ch])
        side = A[~anymask].max()
        row = {"dphi": dphi, "canal": ch}
        pk = [A[m].max() for m in masks]
        for (q, ob, _, nm), v in zip(objs, pk):
            row[f"{nm}_pico_rel_max_dB"] = 20 * np.log10(v / max(pk))
            row[f"{nm}_OSR_dB"] = 20 * np.log10(v / side)
        row["rango_dinamico_dB"] = 20 * np.log10(max(pk) / min(pk))
        row["min_OSR_dB"] = min(row[f"{nm}_OSR_dB"] for *_, nm in objs)
        rowsP.append(row)
dP = pd.DataFrame(rowsP)
dP.to_csv(f"{out}/expD_escenaP2.csv", index=False)
fig, axs = plt.subplots(1, 4, figsize=(13, 3.6))
fig.subplots_adjust(wspace=0.35)
ref = np.abs(I_m).max()
for j, (ti, I) in enumerate([("Mono", I_m), ("Bi 135°", keepP[135]["bi"]), ("Mono+Bi 90°", keepP[90]["fus"]),
                             ("Mono+Bi 135°", keepP[135]["fus"])]):
    A = np.abs(gP.cube(I))
    im = P.show_slice(axs[j], gP.x, gP.y, A.max(axis=2), ref, ti + " (MIP en z)")
    for q, _, _, nm in objs:
        axs[j].text(q[0], q[1] + 0.45, nm, ha="center", fontsize=7, color=P.C["s8"])
fig.colorbar(im, ax=axs, shrink=0.8, label="dB re pico mono")
P.save(fig, "D_escenaP2", out)
print(dP.to_string())
print("Exp. D terminado ->", out)
