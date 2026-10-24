import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from bisar import rcs

lam = 3e8 / 425e6 / 2
k = 2 * np.pi / lam
db = lambda s: 10 * np.log10(s)
nad = np.array([[1e-9, 0, 1.0]])
nad /= np.linalg.norm(nad)
print("Tabla 2 relatorio 2809: cubo 20.04, caja 14.02, cilindro 8.94, esfera -1.54, cono -15.86 dBsm")
for kind, dim in [("cubo", 1.0), ("caja", [1, 0.5, 0.3]), ("cilindro", [0.3, 1.0]), ("esfera", 0.5), ("cono", [0.3, 1.0])]:
    V, F = rcs.make_mesh(kind, dim, lam)
    o = rcs.PO(V, F)
    print(f"  {kind:9s} {len(F):4d} tri  {db(o.sigma(nad, nad, k)[0]):7.2f} dBsm")
V, F = rcs.make_mesh("cubo", 1.0, lam)
o = rcs.PO(V, F)
print("Placa especular (cara superior): PO vs 4 pi A^2 cos^2 / lambda^2")
for th in [0, 15, 30, 45]:
    t = np.radians(th)
    i = np.array([[np.sin(t), 0, np.cos(t)]])
    s = np.array([[-np.sin(t), 0, np.cos(t)]])
    print(f"  {th:2d} deg  {db(o.sigma(i, s, k)[0]):.3f}  {db(4*np.pi*np.cos(t)**2/lam**2):.3f}")
print("Continuidad mono->bi (amplitud compleja TM / TE)")
t = np.radians(20)
i = np.array([[np.sin(t), 0, np.cos(t)]])
for dp in [0, 1e-3, 30]:
    a = np.radians(dp)
    s = np.array([[np.sin(t) * np.cos(a), np.sin(t) * np.sin(a), np.cos(t)]])
    print("  ", dp, np.round(o.amplitude(i, s, k, "TM"), 4), np.round(o.amplitude(i, s, k, "TE"), 4))
