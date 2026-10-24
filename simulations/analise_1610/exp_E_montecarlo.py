"""Experimento E - Monte Carlo de orientacion en SO(3), una orientacion por ADQUISICION.

Para cada repeticion n: R_n ~ U(SO(3)); el mismo objeto rotado se usa en toda la apertura y en
los canales mono y bi (diseno pareado).  Metricas por adquisicion: SNR (pico), contraste, PSLR,
volumen -3 dB, FWHM, error de localizacion horizontal y energia en la region del objeto.
Salidas: expE_mc_<objeto>.csv (todas las repeticiones), expE_resumen.csv (medianas, IC 95 %
bootstrap y probabilidad de mejora), figura E_montecarlo.
Uso: python exp_E_montecarlo.py [objeto ...]
"""
import os
import sys
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bisar import config, geometry, imaging, metrics, plotting as P, rcs, runner, scene

NAME = os.environ.get("BISAR_EXP", "expE")   # expE (Rx TM) o expE2 (Rx alineado)
p = config.load(NAME)
config.save_resolved(p, NAME)
out = config.outdir(NAME)
t, Ptx = geometry.tx_trajectory(p)
M = len(Ptx)
sysm = scene.System(p)
tgt = np.array(p["blanco"])
rg = p["rejilla_mc"]
g = imaging.Grid.around(tgt, rg["hx"], rg["dx"], rg["z_lo"], rg["z_hi"], rg["dz"])
lam_s = p["c"] / p["fc"] / sysm.n_real
Prx = {d: geometry.rx_from_tx(Ptx, d) for d in p["dphi_list"]}
METS = ["SNR_dB", "contrast_dB", "PSLR_dB", "vol3dB_m3", "fwhm_x", "fwhm_z", "loc_err_h_m", "energia_en_region"]
BETTER = {"SNR_dB": +1, "contrast_dB": +1, "PSLR_dB": -1, "vol3dB_m3": -1, "fwhm_x": -1, "fwhm_z": -1,
          "loc_err_h_m": -1, "energia_en_region": +1}


def row_for(I, ch, mask, Imono_pk, M_eff):
    m = metrics.psf_metrics(g, I, tgt)
    E = np.abs(I) ** 2
    return {"canal": ch, "SNR_dB": 10 * np.log10(m["peak_amp"] ** 2 / sysm.noise_var_image(M * M_eff)),
            "contrast_dB": m["contrast_dB"], "PSLR_dB": m["PSLR_dB"], "vol3dB_m3": m["vol3dB_m3"],
            "fwhm_x": m["fwhm_x"], "fwhm_z": m["fwhm_z"],
            "loc_err_h_m": float(np.hypot(m["peak_x"] - tgt[0], m["peak_y"] - tgt[1])),
            "energia_en_region": E[mask].sum() / E.sum(), "peak_rel_mono_dB": 20 * np.log10(m["peak_amp"] / Imono_pk)}


objs = sys.argv[1:] or list(p["objetos"])
rng = np.random.default_rng(p["semilla"])
for name in objs:
    spec = p["objetos"][name]
    V, F = rcs.make_mesh(spec["tipo"], spec["dim"], lam_s)
    Rs = rcs.random_so3(spec["N"], rng)
    rows = []
    T0 = time.time()
    for n, R in enumerate(Rs):
        obj = rcs.PO(V, F, R)
        X = g.vox - tgt
        mask = np.all((X >= obj.V.min(0) - 0.3) & (X <= obj.V.max(0) + 0.3), axis=1)
        s_m, _ = scene.target_echo(sysm, Ptx, Ptx, tgt, obj)
        I_m = runner.bp(sysm, Ptx, Ptx, s_m, g.vox)
        pk = np.abs(I_m).max()
        base = {"objeto": name, "rep": n, "R": R.ravel().round(6).tolist(),
                "Pr_mono_dBm": 10 * np.log10(np.mean(np.abs(s_m) ** 2)) + 30}
        rows.append({**base, "dphi": 0, **row_for(I_m, "mono", mask, pk, 1)})
        for d in p["dphi_list"]:
            s_b, _ = scene.target_echo(sysm, Ptx, Prx[d], tgt, obj)
            I_b = runner.bp(sysm, Ptx, Prx[d], s_b, g.vox)
            b = {**base, "dphi": d, "Pr_bi_dBm": 10 * np.log10(np.mean(np.abs(s_b) ** 2)) + 30}
            rows.append({**b, **row_for(I_b, "bi", mask, pk, 1)})
            rf = row_for(0.5 * (I_m + I_b), "fus", mask, pk, 2)
            rows.append({**b, **rf})
            # mismo presupuesto (N): misma imagen (verificado en Exp. B, corr > 0.999), mitad de muestras
            rows.append({**b, **rf, "canal": "fusN", "SNR_dB": rf["SNR_dB"] - 10 * np.log10(2)})
        if (n + 1) % 25 == 0:
            print(f"{name} {n+1}/{len(Rs)}  {time.time()-T0:.0f}s", flush=True)
    pd.DataFrame(rows).to_csv(f"{out}/expE_mc_{name}.csv", index=False)

# ===================================================================== resumen
bs = np.random.default_rng(1)
summ = []
for name in p["objetos"]:
    try:
        d = pd.read_csv(f"{out}/expE_mc_{name}.csv")
    except FileNotFoundError:
        continue
    mono = d[d.canal == "mono"].set_index("rep")
    for dphi in p["dphi_list"]:
        for ch in ["bi", "fusN"]:
            x = d[(d.canal == ch) & (d.dphi == dphi)].set_index("rep")
            for met in METS:
                diff = (x[met] - mono[met]).values
                if met == "vol3dB_m3":
                    diff = 10 * np.log10(x[met].values / mono[met].values)
                diff = diff[np.isfinite(diff)]
                b = np.median(diff[bs.integers(0, len(diff), (p["n_boot"], len(diff)))], axis=1)
                summ.append({"objeto": name, "dphi": dphi, "canal": ch, "metrica": met, "N": len(diff),
                             "mediana_dif": np.median(diff), "IC95_lo": np.percentile(b, 2.5),
                             "IC95_hi": np.percentile(b, 97.5), "media_dif": diff.mean(),
                             "P_mejora": float(np.mean(BETTER[met] * diff > 0))})
        x = d[(d.canal == "bi") & (d.dphi == dphi)].set_index("rep")
        r = 10 * np.log10(np.mean(10 ** (x.Pr_bi_dBm / 10)) / np.mean(10 ** (mono.Pr_mono_dBm / 10)))
        summ.append({"objeto": name, "dphi": dphi, "canal": "bi", "metrica": "potencia_media_bi_vs_mono_dB",
                     "N": len(x), "mediana_dif": np.median(x.Pr_bi_dBm.values - mono.Pr_mono_dBm.values),
                     "IC95_lo": np.nan, "IC95_hi": np.nan, "media_dif": r,
                     "P_mejora": float(np.mean(x.Pr_bi_dBm.values > mono.Pr_mono_dBm.values))})
ds = pd.DataFrame(summ)
ds.to_csv(f"{out}/expE_resumen.csv", index=False)

# ===================================================================== figura
show = [("SNR_dB", "$\\Delta$SNR fus(N)-mono [dB]"), ("contrast_dB", "$\\Delta$contraste [dB]"),
        ("PSLR_dB", "$\\Delta$PSLR [dB]"), ("vol3dB_m3", "volumen fus/mono [dB]")]
names = [n for n in p["objetos"] if (ds.objeto == n).any()]
fig, axs = plt.subplots(len(names), len(show), figsize=(13, 2.3 * len(names)), squeeze=False)
fig.subplots_adjust(hspace=0.6, wspace=0.3)
for i, nm in enumerate(names):
    d = pd.read_csv(f"{out}/expE_mc_{nm}.csv")
    mono = d[d.canal == "mono"].set_index("rep")
    for j, (met, lab) in enumerate(show):
        data = []
        for dphi in p["dphi_list"]:
            x = d[(d.canal == "fusN") & (d.dphi == dphi)].set_index("rep")
            v = (10 * np.log10(x[met] / mono[met])) if met == "vol3dB_m3" else (x[met] - mono[met])
            data.append(v.values[np.isfinite(v.values)])
        a = axs[i, j]
        bp_ = a.boxplot(data, positions=range(len(data)), widths=0.55, showfliers=False, patch_artist=True)
        for patch in bp_["boxes"]:
            patch.set_facecolor("#a8e3cf"); patch.set_edgecolor(P.C["ink2"])
        for med in bp_["medians"]:
            med.set_color(P.C["ink"])
        a.axhline(0, color=P.C["s8"], lw=0.8)
        a.set_xticks(range(len(data))); a.set_xticklabels([f"{x}°" for x in p["dphi_list"]])
        a.set_title(f"{nm}: {lab}", fontsize=8)
        rr = ds[(ds.objeto == nm) & (ds.canal == "fusN") & (ds.metrica == met)].sort_values("dphi")
        for k, pm in enumerate(rr.P_mejora.values):
            a.text(k, a.get_ylim()[1], f"{100*pm:.0f}%", ha="center", va="bottom", fontsize=6.5, color=P.C["ink2"])
fig.suptitle("Exp. E: Monte Carlo SO(3) por adquisicion; caja = cuartiles de (fusion N - mono); "
             "% = probabilidad de mejora", y=1.0)
P.save(fig, "E_montecarlo" if NAME == "expE" else "E2_montecarlo_rx_alineado", out)
print("Exp. E terminado ->", out)
