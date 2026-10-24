"""Experimento D3 - Biblioteca P0/P1 en orientacion nominal con Rx de polarizacion ALINEADA.

Repite la parte de objetos del Exp. D con pol = 'TMAL' (Tx TM; Rx de doble polarizacion
combinado para recibir la componente del campo del Tx transversal a la direccion de recepcion).
Compara con el Rx TM co-polar del Exp. D (expD_objetos.csv).
Salidas: expD/expD3_objetos_rx_alineado.csv, figura D3_rx_alineado.  Uso: python exp_D3_rx_alineado.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bisar import config, geometry, imaging, plotting as P, rcs, runner, scene

p = config.load("expD")
out = config.outdir("expD")
t, Ptx = geometry.tx_trajectory(p)
M = len(Ptx)
sysm = scene.System(p, pol="TMAL")
tgt = np.array(p["blanco"])
ro = p["rejilla_obj"]
g = imaging.Grid.around(tgt, ro["hx"], ro["dx"], ro["z_lo"], ro["z_hi"], ro["dz"])
lam_s = p["c"] / p["fc"] / sysm.n_real
rows = []
for name, spec in p["objetos"].items():
    if spec["tipo"] == "punto":
        continue
    V, F = rcs.make_mesh(spec["tipo"], spec["dim"], lam_s)
    obj = rcs.PO(V, F, rcs.rot_xyz(*spec.get("rot", [0, 0, 0])))
    s_m, _ = scene.target_echo(sysm, Ptx, Ptx, tgt, obj)
    I_m = runner.bp(sysm, Ptx, Ptx, s_m, g.vox)
    for dphi in [90, 135]:
        Prx = geometry.rx_from_tx(Ptx, dphi)
        s_b, _ = scene.target_echo(sysm, Ptx, Prx, tgt, obj)
        imgs = {"mono": I_m, "bi": runner.bp(sysm, Ptx, Prx, s_b, g.vox)}
        imgs["fus"] = 0.5 * (imgs["mono"] + imgs["bi"])
        for r in runner.psf_rows(sysm, g, imgs, tgt, M, {"objeto": name, "dphi": dphi}):
            r["peak_rel_mono_dB"] = 20 * np.log10(r["peak_amp"] / np.abs(I_m).max())
            if r["canal"] == "fus":           # mismo presupuesto N: misma imagen, mitad de muestras
                r["SNR_dB"] -= 10 * np.log10(2)
            rows.append(r)
    print(name, flush=True)
d = pd.DataFrame(rows)
d.to_csv(f"{out}/expD3_objetos_rx_alineado.csv", index=False)
dTM = pd.read_csv(f"{out}/expD_objetos.csv")
comp = []
for name in d.objeto.unique():
    for dphi in [90, 135]:
        for tag, dd, fch in [("TM", dTM, "fusN"), ("alineado", d, "fus")]:
            x = dd[(dd.objeto == name) & (dd.dphi == dphi)].set_index("canal")
            comp.append({"objeto": name, "dphi": dphi, "rx": tag,
                         "bi_pico_rel_mono_dB": x.loc["bi", "peak_rel_mono_dB"],
                         "dSNR_fusN_dB": x.loc[fch, "SNR_dB"] - x.loc["mono", "SNR_dB"],
                         "dPSLR_dB": x.loc[fch, "PSLR_dB"] - x.loc["mono", "PSLR_dB"],
                         "dContraste_dB": x.loc[fch, "contrast_dB"] - x.loc["mono", "contrast_dB"],
                         "dVol_dB": 10 * np.log10(x.loc[fch, "vol3dB_m3"] / x.loc["mono", "vol3dB_m3"])})
c = pd.DataFrame(comp)
c.to_csv(f"{out}/expD3_comparacion_TM_vs_alineado.csv", index=False)
print(c.round(2).to_string())
fig, axs = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
names = list(d.objeto.unique())
yy = np.arange(len(names))
for a, dphi in zip(axs, [90, 135]):
    for k, (tag, col) in enumerate([("TM", P.C["bi"]), ("alineado", P.C["fus"])]):
        v = [c[(c.objeto == n) & (c.dphi == dphi) & (c.rx == tag)].bi_pico_rel_mono_dB.values[0] for n in names]
        a.scatter(v, yy + (k - 0.5) * 0.2, color=col, s=26, label=f"Rx {tag}", zorder=3)
    a.axvline(0, color=P.C["ink2"], lw=0.8)
    a.set_xlabel("pico del canal bi respecto del mono [dB]")
    a.set_title(f"$\\Delta\\phi$ = {dphi}°")
axs[0].set_yticks(yy); axs[0].set_yticklabels(names); axs[0].invert_yaxis()
axs[0].legend(fontsize=7)
fig.suptitle("Exp. D3: Rx TM co-polar frente a Rx de polarizacion alineada (orientacion nominal)", y=1.03)
P.save(fig, "D3_rx_alineado", out)
