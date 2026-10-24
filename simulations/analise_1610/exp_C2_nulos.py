"""Exp. C2 - ¿Se pueden usar los nulos SPM de la superficie para ubicar el Rx?

Para cada polarizacion (TM, TE co-polares) y separacion dphi, con el Tx en la helice base:
  - imagen de la superficie (especular + difusa SPM, superficie congelada y compartida por
    mono y bi, R realizaciones) en una rejilla alrededor del blanco que incluye la superficie;
  - imagen de cuatro blancos: punto escalar (S = 1, como el Exp. C), esfera PEC, cubo PEC
    nominal y cubo con orientacion aleatoria SO(3) (n_orient_so3 orientaciones).
Metricas: SCR en el voxel del blanco (como el Exp. C) y PSCR = pico del blanco / pico de la
superficie en la rejilla, para el canal bi y para la fusion mono+bi (N), y su diferencia
respecto de la mono (Delta SCR).
CSV: io/analise_1610/expC2/expC2_scr_vs_dphi.csv ; figura C2_nulos_scr.
Uso: python exp_C2_nulos.py
"""
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bisar import config, geometry, imaging, plotting as P, rcs, runner, scene

NAME = "expC2"
p = config.load(NAME)
config.save_resolved(p, NAME)
out = config.outdir(NAME)
tgt = np.array(p["blanco"])
t, Ptx = geometry.tx_trajectory(p)
rg = p["rejilla_scr"]
g = imaging.Grid.around(tgt, rg["hx"], rg["dx"], rg["z_lo"], rg["z_hi"], rg["dz"])
i0 = int(np.argmin(np.linalg.norm(g.vox - tgt, axis=1)))
below = g.vox[:, 2] < -0.5                     # voxeles bajo la superficie (para el pico del blanco)
ru = dict(p["rugosidad"])
fac = scene.surface_facets(ru["patch_half"], ru["facet"], ru["taper"])
R = p["realizaciones"]
seeds = [p["semilla"] + 500 + r for r in range(R)]
lam_s = p["c"] / p["fc"] / np.sqrt(p["suelo"]["eps_r"][0])
objs = {k: rcs.PO(*rcs.make_mesh(v["tipo"], v["dim"], lam_s)) for k, v in p["objetos"].items()}
Rso3 = rcs.random_so3(p["n_orient_so3"], np.random.default_rng(p["semilla"] + 7))
Vc, Fc = rcs.make_mesh("cubo", p["objetos"]["cubo"]["dim"], lam_s)
cubes = [rcs.PO(Vc, Fc, Rk) for Rk in Rso3]


def bp(sysm, Prx, s):
    return runner.bp(sysm, Ptx, Prx, s, g.vox, gate=True)


def channel(sysm, Prx):
    """Imagenes complejas de la superficie (R) y de los blancos para un par Tx-Rx."""
    s_sp, _ = scene.specular_echo(sysm, Ptx, Prx, ru["sigma_h"])
    I_sp = bp(sysm, Prx, s_sp)
    Ic = [I_sp + bp(sysm, Prx, scene.diffuse_echo(sysm, Ptx, Prx, ru, np.random.default_rng(sd), fac))
          for sd in seeds]
    It = {"punto": bp(sysm, Prx, scene.target_echo(sysm, Ptx, Prx, tgt)[0])}
    for k, o in objs.items():
        It[k] = bp(sysm, Prx, scene.target_echo(sysm, Ptx, Prx, tgt, scatter=o)[0])
    It["cubo_so3"] = [bp(sysm, Prx, scene.target_echo(sysm, Ptx, Prx, tgt, scatter=o)[0]) for o in cubes]
    return {"Ic": np.array(Ic), "Isp": I_sp, "It": It}


def metrics(ch):
    Ic = ch["Ic"]
    Pc_vox = np.mean(np.abs(Ic[:, i0]) ** 2)
    Pc_pk = np.mean(np.max(np.abs(Ic) ** 2, axis=1))
    out_ = {"Pc_vox_dB": 10 * np.log10(Pc_vox), "Pc_pico_dB": 10 * np.log10(Pc_pk),
            "Pesp_vox_dB": 10 * np.log10(np.abs(ch["Isp"][i0]) ** 2 + 1e-300)}
    for k, I in ch["It"].items():
        if k == "cubo_so3":
            pv = np.mean([np.abs(Ii[i0]) ** 2 for Ii in I])
            pp = np.mean([np.max(np.abs(Ii[below]) ** 2) for Ii in I])
        else:
            pv, pp = np.abs(I[i0]) ** 2, np.max(np.abs(I[below]) ** 2)
        out_[f"Pt_{k}_dB"] = 10 * np.log10(pv)
        out_[f"SCR_{k}_dB"] = 10 * np.log10(pv / Pc_vox)
        out_[f"PSCR_{k}_dB"] = 10 * np.log10(pp / Pc_pk)
    return out_


def fuse(a, b):
    """Fusion coherente con el mismo presupuesto: media de las dos imagenes (superficie compartida)."""
    It = {k: ([0.5 * (x + y) for x, y in zip(a["It"][k], b["It"][k])] if k == "cubo_so3"
              else 0.5 * (a["It"][k] + b["It"][k])) for k in a["It"]}
    return {"Ic": 0.5 * (a["Ic"] + b["Ic"]), "Isp": 0.5 * (a["Isp"] + b["Isp"]), "It": It}


rows = []
T0 = time.time()
for pol in p["pols"]:
    sysm = scene.System(p, pol=pol)
    mono = channel(sysm, Ptx)
    for dphi in p["dphi_list"]:
        bi = mono if dphi == 0 else channel(sysm, geometry.rx_from_tx(Ptx, dphi))
        for canal, ch in [("bi", bi), ("fus", fuse(mono, bi))]:
            if dphi == 0 and canal == "fus":
                continue
            rows.append({"pol": pol, "dphi": dphi, "canal": "mono" if dphi == 0 else canal, **metrics(ch)})
        print(f"{pol} dphi={dphi:3d}  ({time.time() - T0:.0f}s)", flush=True)

df = pd.DataFrame(rows)
TGT = ["punto", "esfera", "cubo", "cubo_so3"]
for pol in p["pols"]:
    m = df[(df.pol == pol) & (df.canal == "mono")].iloc[0]
    sel = df.pol == pol
    for k in TGT:
        for met in ["SCR", "PSCR", "Pt"]:
            df.loc[sel, f"d{met}_{k}_dB"] = df.loc[sel, f"{met}_{k}_dB"] - m[f"{met}_{k}_dB"]
    df.loc[sel, "dPc_vox_dB"] = df.loc[sel, "Pc_vox_dB"] - m["Pc_vox_dB"]
df.to_csv(f"{out}/expC2_scr_vs_dphi.csv", index=False)

# ------------------------------------------------------------------ figura
LAB = {"punto": "punto escalar (Exp. C)", "esfera": "esfera PEC r=0.5 m", "cubo": "cubo PEC 1 m (nominal)",
       "cubo_so3": f"cubo SO(3) ({p['n_orient_so3']} orient.)"}
COL = {"punto": P.C["ink2"], "esfera": P.C["mono"], "cubo": P.C["bi"], "cubo_so3": P.C["fus"]}
fig, axs = plt.subplots(2, 3, figsize=(14, 7.6), sharex=True)
fig.subplots_adjust(hspace=0.3, wspace=0.28)
for i, pol in enumerate(p["pols"]):
    d = df[(df.pol == pol) & (df.canal.isin(["mono", "bi"]))].sort_values("dphi")
    f_ = df[(df.pol == pol) & (df.canal == "fus")].sort_values("dphi")
    a = axs[i, 0]
    a.plot(d.dphi, d.dPc_vox_dB, color=P.C["s8"], lw=2, label="superficie en el vóxel del blanco")
    for k in TGT:
        a.plot(d.dphi, d[f"dPt_{k}_dB"], color=COL[k], lw=1.4, marker="o", ms=3, label=LAB[k])
    a.set_ylabel(f"{pol}: potencia en imagen respecto de la mono [dB]")
    a.set_title(f"{pol}: lo que cae la superficie frente a lo que cae el blanco", fontsize=9)
    a.set_ylim(-45, 10)
    b = axs[i, 1]
    for k in TGT:
        b.plot(d.dphi, d[f"dSCR_{k}_dB"], color=COL[k], lw=1.6, marker="o", ms=3, label=LAB[k])
        b.plot(f_.dphi, f_[f"dSCR_{k}_dB"], color=COL[k], lw=1.0, ls="--")
    b.axhline(0, color=P.C["ink"], lw=0.8)
    b.set_ylabel("$\\Delta$SCR (vóxel) respecto de la mono [dB]")
    b.set_title(f"{pol}: ganancia de SCR por ubicar el Rx (continua: bi; discontinua: mono+bi)", fontsize=9)
    c = axs[i, 2]
    for k in TGT:
        c.plot(d.dphi, d[f"dPSCR_{k}_dB"], color=COL[k], lw=1.6, marker="o", ms=3)
        c.plot(f_.dphi, f_[f"dPSCR_{k}_dB"], color=COL[k], lw=1.0, ls="--")
    c.axhline(0, color=P.C["ink"], lw=0.8)
    c.set_ylabel("$\\Delta$PSCR (pico/pico) respecto de la mono [dB]")
    c.set_title(f"{pol}: pico del blanco / pico de la superficie", fontsize=9)
    for a_ in axs[i]:
        a_.axvline(90, color=P.C["grid"], lw=0.8); a_.axvline(140, color=P.C["grid"], lw=0.8)
for a_ in axs[1]:
    a_.set_xlabel("$\\Delta\\phi$ Tx-Rx [°]"); a_.set_xticks(range(0, 181, 30))
axs[0, 0].legend(fontsize=7, loc="lower left")
fig.suptitle("Exp. C2: uso de los nulos SPM (HH a 90°, VV ≈ 140°) para ubicar el Rx. Hélice base, S1, "
             "$s$=2 cm, $\\ell$=20 cm; líneas grises en 90° y 140°", y=0.995)
P.save(fig, "C2_nulos_scr", out)
print("Exp. C2 terminado ->", out)
