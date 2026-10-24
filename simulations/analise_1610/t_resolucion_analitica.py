"""Verificacion de las formas cerradas de resolucion (bisar/resolucion.py, Apendice A).

Compara delta_x, delta_y, delta_z cerradas con las medidas en la imagen (back-projection
con refraccion exacta, cortes finos por el pico) para:
  (1) el barrido del Exp. B (dphi = 0:15:180, misma altura y radio), blanco en el eje y descentrado;
  (2) geometrias del Rx fuera del Exp. B (otro radio u otra altura), donde el bi puede superar a la mono.
CSV: io/analise_1610/expB/expB_resolucion_analitica.csv
Uso: python t_resolucion_analitica.py
"""
import numpy as np
import pandas as pd

from bisar import config, geometry, metrics, resolucion, runner, scene

p = config.load("expB")
out = config.outdir("expB")
sysm = scene.System(p)
t, Ptx = geometry.tx_trajectory(p)
fc, B, n = p["fc"], p["B"], sysm.n_real
ven = p.get("ventana_espectral", "rect")


def medidas(Prx, tgt):
    """delta_x, delta_y, delta_z medidas en cortes finos de la imagen (BP exacta)."""
    s, _ = scene.target_echo(sysm, Ptx, Prx, tgt)
    ax, az = np.arange(-3.0, 3.0, 0.002), np.arange(-3.0, 3.0, 0.005)
    out = []
    for v, a in [(np.c_[tgt[0] + ax, 0 * ax + tgt[1], 0 * ax + tgt[2]], ax),
                 (np.c_[0 * ax + tgt[0], tgt[1] + ax, 0 * ax + tgt[2]], ax),
                 (np.c_[0 * az + tgt[0], 0 * az + tgt[1], tgt[2] + az], az)]:
        w, cens = metrics.fwhm_1d(a, np.abs(runner.bp(sysm, Ptx, Prx, s, v)), return_censored=True)
        out.append(np.nan if cens else w)
    return out


def rx(dphi, kappa=1.0, dz=0.0):
    P = geometry.rx_from_tx(Ptx, dphi, kappa)
    P[:, 2] += dz
    return P


casos = [(f"B dphi={d}", d, rx(d)) for d in p["dphi_list"]]
casos += [("radio 1.3x, dphi=0", 0, rx(0, 1.3)), ("radio 1.3x, dphi=30", 30, rx(30, 1.3)),
          ("Rx 60 m mas bajo, dphi=0", 0, rx(0, 1.0, -60.0)), ("Rx casi nadir (0.2x), dphi=0", 0, rx(0, 0.2)),
          ("radio 0.85x, dphi=0", 0, rx(0, 0.85))]

# incidencia extrema de la helice vista desde el blanco en el eje (para la forma 100 % explicita)
_D, _, _h = resolucion.geometria_aparente(Ptx, np.array(p["blancos"]["centro"]))
_psi = np.degrees(resolucion.psi_explicito(_D, _h, -p["blancos"]["centro"][2], n))
PSI_MIN, PSI_MAX = _psi.min(), _psi.max()
print(f"helice vista desde el eje: psi = {PSI_MIN:.2f} .. {PSI_MAX:.2f} deg")

rows = []
for tname, tpos in p["blancos"].items():
    tgt = np.array(tpos)
    for name, dphi, Prx in casos:
        mx, my, mz = medidas(Prx, tgt)
        c1 = resolucion.cerrada(Ptx, Prx, tgt, fc, B, n, ven, iters=1)
        cx = resolucion.cerrada(Ptx, Prx, tgt, fc, B, n, ven, iters=30)
        r = {"blanco": tname, "caso": name, "dphi": dphi,
             "delta_x_medido": mx, "delta_y_medido": my, "delta_z_medido": mz,
             "delta_xy_cerrada": c1["delta_xy"], "delta_z_cerrada": c1["delta_z"],
             "delta_xy_cerrada_psi_exacto": cx["delta_xy"], "delta_z_cerrada_psi_exacto": cx["delta_z"],
             **{k: v for k, v in c1.items() if k not in ("delta_xy", "delta_z")}}
        if tname == "centro" and name.startswith("B "):
            e = resolucion.eje_misma_orbita(PSI_MIN, PSI_MAX, dphi, fc, B, n, ven)
            r["delta_xy_eje_explicita"], r["delta_z_eje_explicita"] = e["delta_xy"], e["delta_z"]
        rows.append(r)
        print(f"{tname:11s} {name:30s} dx med={mx:.3f} cerr={c1['delta_xy']:.3f} ({100 * (c1['delta_xy'] / mx - 1):+5.1f}%) | "
              f"dz med={mz:.3f} cerr={c1['delta_z']:.3f} ({100 * (c1['delta_z'] / mz - 1):+5.1f}%)", flush=True)

df = pd.DataFrame(rows)
df.to_csv(f"{out}/expB_resolucion_analitica.csv", index=False)
ok = df.dphi != 180
for med, cer in [("delta_x_medido", "delta_xy_cerrada"), ("delta_y_medido", "delta_xy_cerrada"),
                 ("delta_z_medido", "delta_z_cerrada")]:
    e = 100 * (df[cer] / df[med] - 1)
    print(f"{cer} vs {med}: max |e| = {np.nanmax(np.abs(e[ok])):.1f} %  (dphi < 180)")
print("->", f"{out}/expB_resolucion_analitica.csv")
