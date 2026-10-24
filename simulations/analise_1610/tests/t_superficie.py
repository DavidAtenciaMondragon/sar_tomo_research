"""Validaciones del modelo de superficie.
1) Suma de facetas de Kirchhoff (plano liso, coherente) == teoria de imagenes (especular).
2) Potencia difusa media (realizaciones) == integral de sigma0 (consistencia del generador).
3) Tiempo del nucleo difuso.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from bisar import config, geometry, scene, em

p = config.load("common")
sysm = scene.System(p)
_, Ptx = geometry.tx_trajectory(p)
m = 400
Tx = Ptx[m:m + 1]
Rx = geometry.rx_from_tx(Tx, 180.0)

# --- 1) plano liso: suma directa de facetas de Kirchhoff en fc
k0 = 2 * np.pi * sysm.fc / sysm.c
lam = sysm.c / sysm.fc
h = 0.04
x = np.arange(-60, 60, h)
y = np.arange(-45, 45, h)
acc = 0j
for xi in x:                      # por filas para limitar memoria
    P = np.c_[np.full(len(y), xi), y, np.zeros(len(y))]
    dT = Tx - P
    RT = np.linalg.norm(dT, axis=1)
    dR = Rx - P
    RR = np.linalg.norm(dR, axis=1)
    i = dT / RT[:, None]
    s = dR / RR[:, None]
    thl = np.arccos(np.clip(np.linalg.norm(i + s, axis=1) / 2, -1, 1))
    R = em.fresnel_r(thl, sysm.eps, "TM")
    cosf = 0.5 * (i[:, 2] + s[:, 2])
    gT = scene.gain(np.repeat(Tx, len(y), 0), -i, sysm.bw_az, sysm.bw_el, sysm.G0)
    gR = scene.gain(np.repeat(Rx, len(y), 0), -s, sysm.bw_az, sysm.bw_el, sysm.G0)
    S = 1j * (2 * np.sqrt(np.pi) / lam) * R * cosf * h * h
    acc += np.sum(np.sqrt(gT * gR) * S * np.exp(-1j * k0 * (RT + RR)) / (RT * RR))
acc *= np.sqrt(sysm.Pt) * lam / (4 * np.pi) ** 1.5
sp, inf = scene.specular_echo(sysm, Tx, Rx)
k_c = np.argmin(np.abs(sysm.freqs - sysm.fc))
print("1) Especular Ref1, theta = %.2f deg" % np.degrees(inf["theta"][0]))
print("   suma de facetas |s|^2 = %.3f dBm ; teoria de imagenes = %.3f dBm"
      % (10 * np.log10(abs(acc) ** 2) + 30, 10 * np.log10(abs(sp[0, k_c]) ** 2) + 30))
print("   (frecuencia de la teoria %.3f MHz vs fc %.3f MHz)" % (sysm.freqs[k_c] / 1e6, sysm.fc / 1e6))

# --- 2) difuso: potencia media vs integral de sigma0
rough = dict(p["rugosidad"])
rng = np.random.default_rng(1)
Tx4 = Ptx[::400]
Rx4 = geometry.rx_from_tx(Tx4, 135.0)
fac = scene.surface_facets(rough["patch_half"], rough["facet"], 0.0)
t0 = time.time()
pw = np.mean([np.mean(np.abs(scene.diffuse_echo(sysm, Tx4, Rx4, rough, rng, fac)) ** 2, axis=1)
              for _ in range(40)], axis=0)
dt = time.time() - t0
F, w, dA = fac
ref = []
for a in range(len(Tx4)):
    dT = Tx4[a] - np.c_[F, 0 * F[:, 0]]
    RT = np.linalg.norm(dT, axis=1)
    dR = Rx4[a] - np.c_[F, 0 * F[:, 0]]
    RR = np.linalg.norm(dR, axis=1)
    i = dT / RT[:, None]
    s = dR / RR[:, None]
    thi, ths = np.arccos(i[:, 2]), np.arccos(s[:, 2])
    cphi = -(i[:, 0] * s[:, 0] + i[:, 1] * s[:, 1]) / (np.hypot(i[:, 0], i[:, 1]) * np.hypot(s[:, 0], s[:, 1]))
    sig0 = em.spm_sigma0(thi, ths, np.arccos(np.clip(cphi, -1, 1)), sysm.eps, "TM", k0, rough["sigma_h"], rough["corr_len"])
    gT = scene.gain(np.repeat(Tx4[a:a + 1], len(F), 0), -i, sysm.bw_az, sysm.bw_el, sysm.G0)
    gR = scene.gain(np.repeat(Rx4[a:a + 1], len(F), 0), -s, sysm.bw_az, sysm.bw_el, sysm.G0)
    lamk = sysm.c / sysm.freqs
    ref.append(np.mean(sysm.Pt * lamk ** 2 / (4 * np.pi) ** 3) * np.sum(sig0 * dA * gT * gR / (RT * RR) ** 2))
print("2) Difuso SPM (dphi = 135): media de 40 realizaciones vs integral de sigma0 [dB]:")
print("   ", np.round(10 * np.log10(pw / np.array(ref)), 2))
print("3) nucleo difuso: %d facetas x %d pulsos x %d frec: %.2f s por realizacion"
      % (len(F), len(Tx4), sysm.Nf, dt / 40))
t0 = time.time()
Rx1600 = geometry.rx_from_tx(Ptx, 135.0)
scene.diffuse_echo(sysm, Ptx, Rx1600, rough, rng, fac)
print("   1600 pulsos: %.1f s" % (time.time() - t0))
