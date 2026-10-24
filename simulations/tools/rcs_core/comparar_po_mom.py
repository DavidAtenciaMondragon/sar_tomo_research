"""Compara PO (rápido, aproximado) con MoM (lento, exacto) en cubo y cilindro PEC."""
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rcs import geometry as g
from rcs.po import POSolver, C0, dbsm
from rcs.mom import MoMSolver

lam = 1.0
f = C0 / lam
ang = np.arange(0, 180.1, 1.0)

casos = [
    ("Cubo 1λ, corte az (el=0°)", g.cube(1.0, lam / 10), "az", 0),
    ("Cubo 2λ, corte az (el=0°)", g.cube(2.0, lam / 10), "az", 0),
    ("Cubo 2λ, corte az (el=30°)", None, "az", 30),
    ("Cilindro r=0.5λ L=2λ, corte el (az=0°)", g.cylinder(0.5, 2.0, lam / 10), "el", 0),
]

fig, axs = plt.subplots(2, 2, figsize=(13, 8.5))
prev = None
resumen = []
for ax, (titulo, mesh, corte, fijo) in zip(axs.ravel(), casos):
    if mesh is None:
        mesh, mom = prev
    else:
        t = time.time()
        mom = MoMSolver(*mesh, f)
        print(f"{titulo}: montaje+LU {time.time() - t:.1f}s")
    prev = (mesh, mom)
    if corte == "az":
        az, el = ang, fijo
        x = ang
    else:
        az, el = fijo, ang - 90
        x = ang - 90
    t = time.time(); r = mom.rcs(az, el); tm = time.time() - t
    t = time.time(); po = POSolver(*mesh).rcs(az, el, f).ravel(); tp = time.time() - t
    vv, hh, vh = r["VV"].ravel(), r["HH"].ravel(), r["VH"].ravel()
    ax.plot(x, dbsm(vv), label="MoM VV", lw=1.6)
    ax.plot(x, dbsm(hh), label="MoM HH", lw=1.6)
    ax.plot(x, dbsm(po), "k--", label="PO", lw=1.3)
    ax.plot(x, dbsm(vh), color="0.6", lw=0.8, label="MoM VH (cross)")
    ax.set_title(titulo, fontsize=10)
    ax.set_xlabel("azimut [°]" if corte == "az" else "elevación [°]")
    ax.set_ylabel("RCS [dBsm]  (λ = 1 m)")
    top = dbsm(max(vv.max(), po.max()))
    ax.set_ylim(top - 50, top + 5)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower center", ncol=4)
    lit = dbsm(po) > top - 25  # zona "fuerte"
    e_all = np.median(np.abs(dbsm(po) - dbsm(vv)))
    e_lit = np.median(np.abs(dbsm(po[lit]) - dbsm(vv[lit])))
    resumen.append((titulo, e_lit, e_all, tm, tp, len(x)))
    print(f"   MoM {len(x)} ángulos x2 pol: {tm:.1f}s | PO: {tp*1e3:.0f} ms | "
          f"|PO-MoM| mediana: {e_lit:.2f} dB (lóbulos fuertes), {e_all:.2f} dB (todo)")

fig.suptitle("RCS monoestática PEC: Óptica Física vs Método de Momentos (bempp-cl)")
fig.tight_layout()
fig.savefig("resultados/po_vs_mom.png", dpi=130)
