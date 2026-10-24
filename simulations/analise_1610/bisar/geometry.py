"""Trayectorias del Tx (transceiver) y del Rx remoto.

La hélice del Tx es la de funcao_espiral.m (gs/src) con los parámetros óptimos de
run_otimiza_tx.m (radarTx_espiral_mono0.json).  El Rx remoto vuela a la MISMA altura
que el Tx en cada instante, z_Rx(t) = z_Tx(t), con radio horizontal kappa(t)*rho_Tx(t)
y acimut phi_Tx(t) + dphi(t)  (Orientações §2.2).
"""
import numpy as np


def helix(N, R_top, R_bottom, H_high, H_low, Vt, PRF, north_offset_deg=0.0):
    """Puerto exacto de gs/src/funcao_espiral.m. Devuelve t, P (M,3)."""
    Vr = Vt * (np.log(R_top) - np.log(R_bottom)) / (N * 2 * np.pi)
    T = (R_top - R_bottom) / Vr
    Vd = (H_high - H_low) / T
    t = np.arange(0.0, T + 0.5 / PRF, 1.0 / PRF)
    t = t[t <= T]
    theta = (Vt / Vr) * (np.log(R_top) - np.log(R_top - Vr * t)) + np.radians(north_offset_deg)
    R = R_top - Vr * t
    P = np.stack([R * np.cos(theta), R * np.sin(theta), H_high - Vd * t], axis=1)
    return t, P


def tx_trajectory(p):
    """Trayectoria del Tx diezmada a p['n_pulsos'] pulsos uniformes en el tiempo.
    IMPORTANTE: igual que GS_snell_test_script.m, se llama funcao_espiral(N, RaioMenor,
    RaioMaior, AltMaior, AltMenor, ...): el radio menor esta ARRIBA.  Asi la linea de base
    de la helice queda perpendicular a la linea de vista (condicion beta = psi0 de otimiza_tx)."""
    h = p["tx_helice"]
    t, P = helix(h["NumVoltas"], h["RaioMenor"], h["RaioMaior"], h["AltMaior"], h["AltMenor"],
                 h["Vt"], h["PRF"], h.get("NorthOffset", 0.0))
    M = p["n_pulsos"]
    idx = np.round(np.linspace(0, len(t) - 1, M)).astype(int)
    return t[idx], P[idx]


def rotate_incidence(P, dtheta_deg, center=(0.0, 0.0, 0.0)):
    """Gira cada punto de la trayectoria en su plano vertical (alrededor de 'center'),
    conservando la distancia al centro y el acimut: theta_i -> theta_i + dtheta.
    Sirve para construir hélices equivalentes centradas en otra incidencia (Exp. C)."""
    c = np.asarray(center, float)
    d = P - c
    rho = np.hypot(d[:, 0], d[:, 1])
    az = np.arctan2(d[:, 1], d[:, 0])
    r = np.linalg.norm(d, axis=1)
    th = np.arctan2(rho, d[:, 2]) + np.radians(dtheta_deg)
    return c + np.stack([r * np.sin(th) * np.cos(az), r * np.sin(th) * np.sin(az), r * np.cos(th)], 1)


def incidence_deg(P, center=(0.0, 0.0, 0.0)):
    d = P - np.asarray(center, float)
    return np.degrees(np.arctan2(np.hypot(d[:, 0], d[:, 1]), d[:, 2]))


def rx_from_tx(Ptx, dphi_deg, kappa=1.0, center=(0.0, 0.0)):
    """Rx a la misma altura que el Tx; acimut desplazado dphi(t) y radio kappa(t)*rho_Tx(t).
    dphi_deg y kappa pueden ser escalares o vectores (M,)."""
    cx, cy = center
    x, y = Ptx[:, 0] - cx, Ptx[:, 1] - cy
    rho = np.hypot(x, y)
    az = np.arctan2(y, x) + np.radians(np.broadcast_to(dphi_deg, rho.shape))
    rr = rho * np.broadcast_to(kappa, rho.shape)
    return np.stack([cx + rr * np.cos(az), cy + rr * np.sin(az), Ptx[:, 2]], 1)


def kinematics(t, P):
    """Velocidad, aceleración y curvatura de una trayectoria muestreada en t."""
    v = np.gradient(P, t, axis=0)
    a = np.gradient(v, t, axis=0)
    sp = np.linalg.norm(v, axis=1)
    acc = np.linalg.norm(a, axis=1)
    curv = np.linalg.norm(np.cross(v, a), axis=1) / np.maximum(sp, 1e-12) ** 3
    return sp, acc, curv


def familia_C(P0, nombre, theta_c):
    """Familias de trayectoria del Tx del Exp. C, construidas a partir de la helice base P0.
    F1: cada punto de P0 se gira en su plano vertical (mismo acimut y misma distancia al centro)
        theta_c - <theta_0>; conserva la dispersion propia de incidencia de la helice (+-7 deg).
    F2: cada punto se coloca a la misma distancia al centro y mismo acimut, pero con incidencia
        theta_c constante (cono de incidencia constante = restriccion pseudo-Brewster)."""
    if nombre == "F1":
        return rotate_incidence(P0, theta_c - incidence_deg(P0).mean())
    th = np.radians(theta_c)
    az = np.arctan2(P0[:, 1], P0[:, 0])
    r = np.linalg.norm(P0, axis=1)
    return np.stack([r * np.sin(th) * np.cos(az), r * np.sin(th) * np.sin(az), r * np.cos(th)], 1)
