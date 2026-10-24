"""RCS vs frecuencia: muestra las regiones Rayleigh / resonancia / óptica y dónde vale PO."""
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rcs import geometry as g
from rcs.po import POSolver, C0, dbsm
from rcs.mom import MoMSolver
from rcs.analytic import sphere_pec_monostatic

fs = np.geomspace(0.1e9, 20e9, 300)
f_mom = np.array([0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.25]) * 1e9
fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))

# Esfera r = 0.1 m
r = 0.1
V, F = g.sphere(r, r, curv=int(np.ceil(2 * np.pi / np.sqrt(0.08 * C0 / fs[-1] / r))))
po = POSolver(V, F).rcs_freq(0, 0, fs).ravel()
mie = [sphere_pec_monostatic(r, C0 / f) for f in fs]
ax = axs[0]
ax.plot(fs / 1e9, dbsm(np.array(mie)), label="Mie (exacto)")
ax.plot(fs / 1e9, dbsm(po), "k--", label="PO")
ax.set_title("Esfera PEC r = 0.1 m")

# Cubo 0.3 m: cara frontal (az=0, el=0) y fuera de especular (az=20, el=30)
a = 0.3
V, F = g.cube(a, a / 2)
dirs = [((0,), (0,), "cara (az 0°, el 0°)"), ((20,), (30,), "az 20°, el 30°")]
for ax_i, (az, el, name) in zip([1, 2], dirs):
    po = POSolver(V, F).rcs_freq(az, el, fs).ravel()
    axs[ax_i].plot(fs / 1e9, dbsm(po), "k--", label="PO")
    axs[ax_i].set_title(f"Cubo PEC 0.3 m — {name}")
mom = {1: {"VV": [], "HH": []}, 2: {"VV": [], "HH": []}}
for f in f_mom:
    t = time.time()
    Vm, Fm = g.cube(a, min(C0 / f / 10, a / 6))
    m = MoMSolver(Vm, Fm, f, verbose=False)
    for ax_i, (az, el, _) in zip([1, 2], dirs):
        rr = m.rcs(az, el)
        for p in ("VV", "HH"):
            mom[ax_i][p].append(rr[p][0, 0])
    print(f"MoM cubo f={f/1e9:.2f} GHz ({m.ndof} incógnitas): {time.time()-t:.0f}s")
for ax_i in (1, 2):
    axs[ax_i].plot(f_mom / 1e9, dbsm(np.array(mom[ax_i]["VV"])), "o", label="MoM VV")
    axs[ax_i].plot(f_mom / 1e9, dbsm(np.array(mom[ax_i]["HH"])), "s", mfc="none", label="MoM HH")

for ax in axs:
    ax.set_xscale("log")
    ax.set_xlabel("frecuencia [GHz]")
    ax.set_ylabel("RCS [dBsm]")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
axs[0].set_ylim(-45, -5)
fig.tight_layout()
fig.savefig("resultados/barrido_frecuencia.png", dpi=120)
