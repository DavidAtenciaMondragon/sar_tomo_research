"""Experimento F - Diseno de la trayectoria del Rx para imagen (frontera de Pareto) y robustez.

Parte 1: evaluacion de candidatos (dphi0, A, kappa) -> metricas separadas:
   Q_img  = ISLR_mono - ISLR_fus [dB]          (>0: la fusion concentra mejor la energia)
   A_PSF  = 10 log10(V_fus / V_mono) [dB]       (>0: lobulo principal 3D mas grande)
   SCR    = SCR_fus - SCR_mono [dB]             (supresion del retorno de superficie, S1, s = 2 cm, TM)
   C_voo  = (v_max,Rx / v_Tx - 1)_+ + |L_Rx/L_Tx - 1|   (+ inviable si separacion < sep_min)
   PSLR, FWHM, anisotropia, aceleracion: se reportan por separado.
Parte 2: frontera de Pareto, J ponderado (pesos explicitos y alternativos), trayectoria final.
Parte 3: robustez de la trayectoria final (posicion, fase, sincronizacion, suelo).
Uso: python exp_F_trayectoria.py
"""
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bisar import config, geometry, imaging, metrics, plotting as P, runner, scene

NAME = "expF"
p = config.load(NAME)
config.save_resolved(p, NAME)
out = config.outdir(NAME)
sysm = scene.System(p)
tgt = np.array(p["blanco"])
t, Ptx = geometry.tx_trajectory(p)
M = len(Ptx)
T = t[-1]
rg = p["rejilla_psf"]
g = imaging.Grid.around(tgt, rg["hx"], rg["dx"], rg["z_lo"], rg["z_hi"], rg["dz"])
vn = p["vecindad_scr"]
gN = imaging.Grid.around(tgt, vn["h"], vn["dx"], -vn["hz"], vn["hz"], vn["dz"])
i0 = np.argmin(np.linalg.norm(gN.vox - tgt, axis=1))
rough = dict(p["rugosidad"])
fac = scene.surface_facets(rough["patch_half"], rough["facet"], rough["taper"])
seeds = [p["semilla"] + 5000 + r for r in range(p["realizaciones"])]
vtx, _, _ = geometry.kinematics(t, Ptx)
Ltx = np.sum(np.linalg.norm(np.diff(Ptx, axis=0), axis=1))
scale = vtx.mean() / p["v_ref_tx"]           # factor de tiempo para un UAV real (v_Tx = v_ref)


def clutter_images(Ptx_, Prx_):
    s_t, _ = scene.target_echo(sysm, Ptx_, Prx_, tgt)
    It = runner.bp(sysm, Ptx_, Prx_, s_t, gN.vox, gate=True)[i0]
    s_s, _ = scene.specular_echo(sysm, Ptx_, Prx_, rough["sigma_h"])
    Ic = [runner.bp(sysm, Ptx_, Prx_, s_s + scene.diffuse_echo(sysm, Ptx_, Prx_, rough, np.random.default_rng(sd), fac),
                    gN.vox, gate=True) for sd in seeds]
    return It, np.array(Ic)


def scr(It, Ic):
    return 10 * np.log10(abs(It) ** 2 / np.mean(np.abs(Ic) ** 2))


# ---------------------------------------------------------------- referencia mono
s_m, _ = scene.target_echo(sysm, Ptx, Ptx, tgt)
I_m = runner.bp(sysm, Ptx, Ptx, s_m, g.vox, gate=True)
mm = metrics.psf_metrics(g, I_m, tgt)
It_m, Ic_m = clutter_images(Ptx, Ptx)
SCR_m = scr(It_m, Ic_m)
print(f"mono: ISLR {mm['ISLR_dB']:.2f}  PSLR {mm['PSLR_dB']:.2f}  V {mm['vol3dB_m3']*1e3:.2f} l  SCR {SCR_m:.2f} dB")

# ---------------------------------------------------------------- Parte 1
rows = []
T0 = time.time()
cand = [(d0, A, k) for d0 in p["dphi0_list"] for A in p["barrido_A_list"] for k in p["kappa_list"]]
for ic, (d0, A, kap) in enumerate(cand):
    dphi_t = d0 + A * (2 * t / T - 1)
    Prx = geometry.rx_from_tx(Ptx, dphi_t, kap)
    s_b, _ = scene.target_echo(sysm, Ptx, Prx, tgt)
    I_b = runner.bp(sysm, Ptx, Prx, s_b, g.vox, gate=True)
    mf = metrics.psf_metrics(g, 0.5 * (I_m + I_b), tgt)
    mb = metrics.psf_metrics(g, I_b, tgt)
    It_b, Ic_b = clutter_images(Ptx, Prx)
    SCR_f = scr(0.5 * (It_m + It_b), 0.5 * (Ic_m + Ic_b))
    v, a, cv = geometry.kinematics(t, Prx)
    Lrx = np.sum(np.linalg.norm(np.diff(Prx, axis=0), axis=1))
    sep = np.linalg.norm(Prx - Ptx, axis=1)
    vr = v.max() / vtx.mean()
    C = max(vr - 1, 0) + abs(Lrx / Ltx - 1)
    rows.append({"dphi0": d0, "A": A, "kappa": kap, "dphi_min": dphi_t.min(), "dphi_max": dphi_t.max(),
                 "Q_img_dB": mm["ISLR_dB"] - mf["ISLR_dB"], "A_PSF_dB": 10 * np.log10(mf["vol3dB_m3"] / mm["vol3dB_m3"]),
                 "SCR_gain_dB": SCR_f - SCR_m, "C_voo": C, "SCR_fus_dB": SCR_f, "SCR_bi_dB": scr(It_b, Ic_b),
                 "PSLR_fus_dB": mf["PSLR_dB"], "ISLR_fus_dB": mf["ISLR_dB"], "vol_fus_l": mf["vol3dB_m3"] * 1e3,
                 "fwhm_x_fus": mf["fwhm_x"], "fwhm_y_fus": mf["fwhm_y"], "fwhm_z_fus": mf["fwhm_z"],
                 "aniso3d_fus": mf["aniso3d"], "aniso_xy_fus": mf["aniso_xy"], "contrast_fus_dB": mf["contrast_dB"],
                 "PSLR_bi_dB": mb["PSLR_dB"], "fwhm_x_bi": mb["fwhm_x"],
                 "v_ratio_max": vr, "L_ratio": Lrx / Ltx, "sep_min_m": sep.min(), "sep_max_m": sep.max(),
                 "a_max_ms2_vref": a.max() / scale ** 2, "curv_max": cv.max(),
                 "viable": bool(sep.min() >= p["sep_min_m"] and vr <= 1.5)})
    print(f"[{ic+1}/{len(cand)}] d0={d0} A={A} k={kap}: Q={rows[-1]['Q_img_dB']:.2f} A={rows[-1]['A_PSF_dB']:.2f} "
          f"SCR={rows[-1]['SCR_gain_dB']:.2f} C={C:.3f}  ({time.time()-T0:.0f}s)", flush=True)
df = pd.DataFrame(rows)

# ---------------------------------------------------------------- Parte 2: Pareto y J
obj = np.c_[-df.Q_img_dB, -df.SCR_gain_dB, df.A_PSF_dB, df.C_voo]       # todos a minimizar
viable = df.viable.values
par = np.zeros(len(df), bool)
for i in range(len(df)):
    if not viable[i]:
        continue
    dom = np.all(obj[viable] <= obj[i], axis=1) & np.any(obj[viable] < obj[i], axis=1)
    par[i] = not dom.any()
df["pareto"] = par
nz = lambda x: (x - x[viable].min()) / (x[viable].max() - x[viable].min() + 1e-12)
Qn, Sn, An, Cn = nz(df.Q_img_dB.values), nz(df.SCR_gain_dB.values), nz(df.A_PSF_dB.values), nz(df.C_voo.values)


def J(w):
    j = w["w1_Q"] * Qn + w["w2_SCR"] * Sn - w["w3_A"] * An - w["w4_C"] * Cn
    return np.where(viable, j, -np.inf)


df["J"] = J(p["pesos"])
for k, w in enumerate(p["pesos_alternativos"]):
    df[f"J_alt{k+1}"] = J(w)
df.to_csv(f"{out}/expF_candidatos.csv", index=False)
best = df.loc[df.J.idxmax()]
choice = [{"pesos": "principal", **p["pesos"], **best[["dphi0", "A", "kappa", "Q_img_dB", "SCR_gain_dB", "A_PSF_dB", "C_voo"]].to_dict()}]
for k, w in enumerate(p["pesos_alternativos"]):
    b = df.loc[df[f"J_alt{k+1}"].idxmax()]
    choice.append({"pesos": f"alt{k+1}", **w, **b[["dphi0", "A", "kappa", "Q_img_dB", "SCR_gain_dB", "A_PSF_dB", "C_voo"]].to_dict()})
pd.DataFrame(choice).to_csv(f"{out}/expF_eleccion_pesos.csv", index=False)
print(pd.DataFrame(choice).to_string())

# --- Fig. 11 Pareto
fig, axs = plt.subplots(1, 3, figsize=(13, 3.9))
cmap = plt.get_cmap("viridis")
for a, (xc, yc, xl, yl) in zip(axs, [("C_voo", "Q_img_dB", "costo de vuelo $C_{voo}$", "$Q_{img}$ = ISLR$_{mono}$-ISLR$_{fus}$ [dB]"),
                                     ("SCR_gain_dB", "Q_img_dB", "$\\Delta$SCR fus-mono [dB]", "$Q_{img}$ [dB]"),
                                     ("A_PSF_dB", "Q_img_dB", "$A_{PSF}$ = vol fus/mono [dB]", "$Q_{img}$ [dB]")]):
    sc = a.scatter(df[xc][viable], df[yc][viable], c=df.dphi0[viable], cmap=cmap, s=18, vmin=60, vmax=165)
    a.scatter(df[xc][~viable], df[yc][~viable], marker="x", color=P.C["ink2"], s=14, label="inviable")
    a.scatter(df[xc][par], df[yc][par], facecolor="none", edgecolor=P.C["s8"], s=55, lw=1.2, label="Pareto")
    a.scatter([best[xc]], [best[yc]], marker="*", color=P.C["s8"], s=160, label="J max (pesos principales)")
    a.set_xlabel(xl); a.set_ylabel(yl)
axs[0].legend(fontsize=7)
fig.colorbar(sc, ax=axs, label="$\\Delta\\phi_0$ [°]")
fig.suptitle("Exp. F: frontera de Pareto (sin ningun termino de potencia del eco)", y=1.02)
P.save(fig, "F_pareto", out)

# --- Fig. 12 trayectoria final vs Ref1 y Ref2
d0, A, kap = best.dphi0, best.A, best.kappa
Pbest = geometry.rx_from_tx(Ptx, d0 + A * (2 * t / T - 1), kap)
Pref1 = geometry.rx_from_tx(Ptx, 180.0)
fig = plt.figure(figsize=(11.5, 4.4))
fig.subplots_adjust(wspace=0.3)
ax = fig.add_subplot(1, 2, 1, projection="3d")
ax.plot(*Ptx.T, color=P.C["mono"], lw=1.4, label="Tx transceiver = Ref2")
ax.plot(*Pref1.T, color=P.C["bi"], lw=1.0, ls=":", label="Rx Ref1 (180°)")
ax.plot(*Pbest.T, color=P.C["fus"], lw=1.6, label=f"Rx propuesto ($\\Delta\\phi_0$={d0:g}°, A={A:g}°, $\\kappa$={kap:g})")
ax.scatter(*tgt, color="k", s=12)
ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_zlabel("z [m]")
ax.legend(fontsize=7, loc="upper left"); ax.view_init(28, -55)
ax2 = fig.add_subplot(1, 2, 2)
ax2.plot(t * scale, (d0 + A * (2 * t / T - 1)) * np.ones_like(t), color=P.C["fus"], label="$\\Delta\\phi(t)$ propuesto")
ax2.axhline(180, color=P.C["bi"], ls=":", label="Ref1")
ax2.axhline(0, color=P.C["mono"], label="Ref2")
ax2.set_xlabel(f"tiempo [s] (Tx a {p['v_ref_tx']:g} m/s)"); ax2.set_ylabel("$\\Delta\\phi$ [°]")
ax2.set_ylim(-10, 190); ax2.legend(fontsize=7)
P.save(fig, "F_trayectoria_final", out)

# ---------------------------------------------------------------- Parte 3: robustez
rb = p["robustez"]
dphi_t = d0 + A * (2 * t / T - 1)
s_b0, _ = scene.target_echo(sysm, Ptx, Pbest, tgt)
rng = np.random.default_rng(p["semilla"] + 9)
nom_f = metrics.psf_metrics(g, 0.5 * (I_m + runner.bp(sysm, Ptx, Pbest, s_b0, g.vox, gate=True)), tgt)
rrows = [{"perturbacion": "nominal", "valor": 0, "canal": "fus", **nom_f}, {"perturbacion": "nominal", "valor": 0, "canal": "mono", **mm}]


def add(label, val, s_b, Prx_bp=Pbest, n_bp=None, s_mono=None, Imono=None):
    Ib = runner.bp(sysm, Ptx, Prx_bp, s_b, g.vox, n_bp=n_bp, gate=True)
    Im = I_m if Imono is None else Imono
    r = metrics.psf_metrics(g, 0.5 * (Im + Ib), tgt)
    rrows.append({"perturbacion": label, "valor": val, "canal": "fus", **r})
    if Imono is not None:
        rrows.append({"perturbacion": label, "valor": val, "canal": "mono", **metrics.psf_metrics(g, Imono, tgt)})


for sgm in rb["pos_sigma_m"]:                       # ruido de posicion del Rx (independiente por pulso)
    Ptrue = Pbest + rng.normal(0, sgm, Pbest.shape)
    add("pos_ruido_m", sgm, scene.target_echo(sysm, Ptx, Ptrue, tgt)[0])
for sgm in rb["pos_paseo_m"]:                       # deriva suave (paseo aleatorio, rms = sgm)
    w = np.cumsum(rng.normal(0, 1, Pbest.shape), axis=0)
    w = w - w.mean(0)
    w *= sgm / np.sqrt(np.mean(w ** 2))
    add("pos_deriva_m", sgm, scene.target_echo(sysm, Ptx, Pbest + w, tgt)[0])
for sg in rb["fase_sigma_deg"]:                     # ruido de fase del oscilador del Rx
    ph = np.exp(1j * np.radians(rng.normal(0, sg, M)))[:, None]
    add("fase_deg", sg, s_b0 * ph)
for dtn in rb["sinc_offset_ns"]:                    # error de sincronizacion constante
    add("sinc_offset_ns", dtn, s_b0 * np.exp(-2j * np.pi * np.outer(np.full(M, dtn * 1e-9), sysm.freqs)))
for dtn in rb["sinc_deriva_ns"]:                    # deriva lineal del reloj 0 -> dtn
    add("sinc_deriva_ns", dtn, s_b0 * np.exp(-2j * np.pi * np.outer(np.linspace(0, dtn * 1e-9, M), sysm.freqs)))
for er in rb["n_error_rel"]:                        # indice del suelo mal conocido en el BP (afecta a ambos canales)
    nb = sysm.n_real * (1 + er)
    Imono_e = runner.bp(sysm, Ptx, Ptx, s_m, g.vox, n_bp=nb, gate=True)
    add("n_error_rel", er, s_b0, n_bp=nb, Imono=Imono_e)
dr = pd.DataFrame(rrows)
dr.to_csv(f"{out}/expF_robustez.csv", index=False)

fig, axs = plt.subplots(1, 3, figsize=(13, 3.4))
groups = [("pos_ruido_m", "ruido posicion Rx [cm]", 100), ("pos_deriva_m", "deriva posicion Rx [cm]", 100),
          ("fase_deg", "ruido de fase Rx [°]", 1), ("sinc_offset_ns", "offset de sincronismo [ns]", 1),
          ("sinc_deriva_ns", "deriva de sincronismo [ns]", 1)]
for (gname, lab, fct), col in zip(groups, P.SERIES):
    d = dr[(dr.perturbacion == gname) & (dr.canal == "fus")].sort_values("valor")
    x = np.arange(len(d))
    for a, met in zip(axs, ["PSLR_dB", "ISLR_dB", "loc_err_m"]):
        a.plot(d.valor * fct, d[met], marker="o", color=col, label=lab)
for a, met, yl in zip(axs, ["PSLR_dB", "ISLR_dB", "loc_err_m"], ["PSLR fus [dB]", "ISLR fus [dB]", "error de localizacion [m]"]):
    a.axhline(nom_f[met], color=P.C["fus"], ls="--", lw=0.8, label="fusion nominal" if met == "PSLR_dB" else None)
    a.axhline(mm[met], color=P.C["mono"], lw=1.2, label="mono" if met == "PSLR_dB" else None)
    a.set_xscale("symlog", linthresh=1); a.set_xlabel("magnitud de la perturbacion"); a.set_ylabel(yl)
axs[0].legend(fontsize=6.5)
fig.suptitle("Exp. F: robustez de la fusion con la trayectoria propuesta (errores solo en el canal bi)", y=1.03)
P.save(fig, "F_robustez", out)
print(dr[["perturbacion", "valor", "canal", "PSLR_dB", "ISLR_dB", "vol3dB_m3", "loc_err_m", "peak_amp"]].to_string())
print("Exp. F terminado ->", out)
