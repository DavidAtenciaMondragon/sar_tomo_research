"""Validacion cruzada del motor BP contra el pipeline MATLAB GS/PROC existente.

Caso: parametros *_espiral_mono0 (fc = 10 GHz, B = 50 MHz, n = 2, blanco (0,0,-5),
helice Tx = Rx).  Resultado MATLAB (io/snell/processed/bistatic_sar_snell_resolution_3dB_espiral_mono0.csv):
    dX = 0.01555 m, dY = 0.01571 m, dZ = 0.1997 m
y offset180 (Rx con otra helice, NorthOffset 180): dX = 0.1912 m, dZ = 0.2056 m.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from bisar import geometry, imaging, metrics

c = 3e8
fc = 10e9
n = 2.0
tgt = np.array([0, 0, -5.0])
freqs = fc + np.linspace(-25e6, 25e6, 64)


def run(rx_par, label, matlab):
    t, Ptx = geometry.helix(2, 137.40, 165.53, 131.45, 93.58, 113.50, 400, 0)
    _, Prx = geometry.helix(2, *rx_par, 400, rx_par[-1] if False else 0)
    return t, Ptx, Prx


def echo_point(Ptx, Prx):
    # misma amplitud constante que GS_snell_test_script.m (rcs = 1, sin ecuacion radar)
    from bisar.scene import refract
    a = refract(Ptx, tgt, n)
    b = refract(Prx, tgt, n)
    tau = (a["R1"] + n * a["R2"] + b["R1"] + n * b["R2"]) / c
    return np.exp(-2j * np.pi * np.outer(tau, freqs))


def widths(Ptx, Prx):
    s = echo_point(Ptx, Prx)
    ax = np.arange(-0.06, 0.0601, 0.0005)
    az = np.arange(-0.4, 0.4001, 0.002)
    out = []
    for lab, vox, axis in [("X", np.c_[ax, 0 * ax, tgt[2] + 0 * ax], ax),
                           ("Y", np.c_[0 * ax, ax, tgt[2] + 0 * ax], ax),
                           ("Z", np.c_[0 * az, 0 * az, tgt[2] + az], az)]:
        I = imaging.backproject(Ptx, Prx, s, freqs, vox, n, c)
        out.append(metrics.fwhm_1d(axis, np.abs(I)))
    return out


t0 = time.time()
_, Ptx = geometry.helix(2, 137.40, 165.53, 131.45, 93.58, 113.50, 400, 0)
w = widths(Ptx, Ptx)
print(f"mono0 (Rx = Tx): dX = {w[0]:.5f}  dY = {w[1]:.5f}  dZ = {w[2]:.4f} m   | MATLAB 0.01555 0.01571 0.1997")
_, Prx = geometry.helix(2, 147.5, 172.5, 120, 80, 120, 400, 180)
m = min(len(Ptx), len(Prx))
w = widths(Ptx[:m], Prx[:m])
print(f"offset180      : dX = {w[0]:.5f}  dY = {w[1]:.5f}  dZ = {w[2]:.4f} m   | MATLAB 0.1912 0.1948 0.2056 (rejilla 0.04 m)")
_, Prx = geometry.helix(2, 147.5, 172.5, 120, 80, 120, 400, 0)
m = min(len(Ptx), len(Prx))
w = widths(Ptx[:m], Prx[:m])
print(f"offset0        : dX = {w[0]:.5f}  dY = {w[1]:.5f}  dZ = {w[2]:.4f} m   | MATLAB 0.01518 0.01519 0.2056")
print(f"tiempo {time.time()-t0:.1f} s, pulsos {len(Ptx)}")
