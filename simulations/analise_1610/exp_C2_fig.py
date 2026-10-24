"""Exp. C2 - figura a partir del CSV (sin recalcular). Uso: python exp_C2_fig.py"""
import pandas as pd
import matplotlib.pyplot as plt

from bisar import config, plotting as P

out = config.outdir("expC2")
df = pd.read_csv(f"{out}/expC2_scr_vs_dphi.csv")
n_so3 = config.load("expC2")["n_orient_so3"]
TGT = ["punto", "esfera", "cubo", "cubo_so3"]
LAB = {"punto": "punto escalar (Exp. C)", "esfera": "esfera PEC r=0.5 m", "cubo": "cubo PEC 1 m (nominal)",
       "cubo_so3": f"cubo SO(3) ({n_so3} orient.)"}
COL = {"punto": P.C["ink2"], "esfera": P.C["mono"], "cubo": P.C["bi"], "cubo_so3": P.C["fus"]}

fig, axs = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
fig.subplots_adjust(hspace=0.32, wspace=0.3)
for i, pol in enumerate(["TM", "TE"]):
    d = df[(df.pol == pol) & (df.canal.isin(["mono", "bi"]))].sort_values("dphi")
    f_ = df[(df.pol == pol) & (df.canal == "fus")].sort_values("dphi")
    a = axs[i, 0]
    a.plot(d.dphi, d.dPc_vox_dB, color=P.C["s8"], lw=2, label="superficie (vóxel del blanco)")
    for k in TGT:
        a.plot(d.dphi, d[f"dPt_{k}_dB"], color=COL[k], lw=1.4, marker="o", ms=3, label=LAB[k])
    a.set_ylabel("potencia bi / mono [dB]")
    a.set_title(f"({'ad'[i]}) {pol}: caída de superficie y blanco", fontsize=10)
    a.set_ylim(-45, 10)
    b = axs[i, 1]
    for k in TGT:
        b.plot(d.dphi, d[f"dSCR_{k}_dB"], color=COL[k], lw=1.6, marker="o", ms=3)
        b.plot(f_.dphi, f_[f"dSCR_{k}_dB"], color=COL[k], lw=1.0, ls="--")
    b.axhline(0, color=P.C["ink"], lw=0.8)
    b.set_ylabel("$\\Delta$SCR vóxel [dB]")
    b.set_title(f"({'be'[i]}) {pol}: SCR respecto de la mono", fontsize=10)
    c = axs[i, 2]
    for k in TGT:
        c.plot(d.dphi, d[f"dPSCR_{k}_dB"], color=COL[k], lw=1.6, marker="o", ms=3)
        c.plot(f_.dphi, f_[f"dPSCR_{k}_dB"], color=COL[k], lw=1.0, ls="--")
    c.axhline(0, color=P.C["ink"], lw=0.8)
    c.set_ylabel("$\\Delta$PSCR pico/pico [dB]")
    c.set_title(f"({'cf'[i]}) {pol}: pico blanco / pico superficie", fontsize=10)
    for a_ in axs[i]:
        a_.axvline(90, color=P.C["grid"], lw=0.8)
        a_.axvline(140, color=P.C["grid"], lw=0.8)
for a_ in axs[1]:
    a_.set_xlabel("$\\Delta\\phi$ Tx-Rx [°]")
    a_.set_xticks(range(0, 181, 30))
axs[0, 0].legend(fontsize=7, loc="lower left")
P.save(fig, "C2_nulos_scr", out)
print("figura C2 regenerada ->", out)
