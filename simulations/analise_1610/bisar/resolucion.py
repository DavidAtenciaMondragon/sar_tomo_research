"""Resolucion analitica (formas cerradas) de la PSF mono / biestatica de un blanco enterrado.

Notacion (Apendice A del relatorio_1610):
  P = (x_P, y_P, z_P), d = -z_P > 0 profundidad;  plataforma i in {Tx, Rx} en (x_i, y_i, h_i).
  D_i      distancia horizontal plataforma-blanco
  phi'_i   acimut aparente de la plataforma visto desde el blanco
  psi_i    incidencia en el aire en el punto de refraccion (Snell: sin psi = n sin theta)
  S        = |sin psi_T e(phi'_T) + sin psi_R e(phi'_R)| / 2        (mono: S = sin psi)
  q        = sqrt(n^2 - sin^2 psi_T) + sqrt(n^2 - sin^2 psi_R)      (= n (cos th_T + cos th_R))

  delta_x = delta_y = 2 x_J0 / (k0 <2S>) = x_J0 lambda0 / (2 pi <S>),   x_J0 = 1.1264 (J0(x) = 1/sqrt 2)
  delta_z = c_w c / sqrt( (B <q>)^2 + 12 f_c^2 var(q) ),               c_w = 0.886 (rect), 1.443 (Hann)

Las resoluciones son anchos a -3 dB de |PSF| (amplitud -3 dB = mitad de potencia).
"""
import numpy as np

X_J0 = 1.1264                      # J0(X_J0) = 1/sqrt(2)
C_W = {"rect": 0.886, "hann": 1.443}   # ancho -3 dB de sinc / Hann en unidades de 1/B


def psi_explicito(D, h, d, n, iters=1):
    """Incidencia en el aire por Snell en forma explicita: D = h tan psi + d tan theta.
    iters=1 es la forma cerrada del apendice; iters grandes -> solucion exacta."""
    psi = np.arctan2(D, h)
    for _ in range(iters):
        s = np.sin(psi)
        tan_th = s / np.sqrt(n ** 2 - s ** 2)
        psi = np.arctan2(D - d * tan_th, h)
    return psi


def geometria_aparente(Pplat, tgt):
    """D_i, phi'_i y h_i de cada pulso, vistos desde el blanco."""
    dx, dy = Pplat[:, 0] - tgt[0], Pplat[:, 1] - tgt[1]
    return np.hypot(dx, dy), np.arctan2(dy, dx), Pplat[:, 2]


def cerrada(Ptx, Prx, tgt, fc, B, n, ventana="rect", c=3e8, iters=1):
    """Resoluciones -3 dB cerradas para un blanco en tgt.  Devuelve dict con los intermedios."""
    d = -float(tgt[2])
    DT, aT, hT = geometria_aparente(Ptx, tgt)
    DR, aR, hR = geometria_aparente(Prx, tgt)
    pT = psi_explicito(DT, hT, d, n, iters)
    pR = psi_explicito(DR, hR, d, n, iters)
    sT, sR = np.sin(pT), np.sin(pR)
    S = 0.5 * np.sqrt(np.maximum(sT ** 2 + sR ** 2 + 2 * sT * sR * np.cos(aR - aT), 0.0))
    q = np.sqrt(n ** 2 - sT ** 2) + np.sqrt(n ** 2 - sR ** 2)
    lam0 = c / fc
    dxy = X_J0 * lam0 / (2 * np.pi * S.mean()) if S.mean() > 1e-9 else np.inf
    dz = C_W[ventana] * c / np.sqrt((B * q.mean()) ** 2 + 12 * fc ** 2 * q.var())
    dap = np.degrees(np.abs(np.angle(np.exp(1j * (aR - aT)))))
    return {"delta_xy": dxy, "delta_z": dz, "S_media": S.mean(), "q_media": q.mean(), "q_std": q.std(),
            "psi_T_media_deg": np.degrees(pT).mean(), "psi_R_media_deg": np.degrees(pR).mean(),
            "dphi_ap_min": dap.min(), "dphi_ap_max": dap.max()}


def eje_misma_orbita(psi_min, psi_max, dphi_deg, fc, B, n, ventana="rect", c=3e8):
    """Caso particular 100 % explicito (ecs. A.12-A.13): blanco en el eje, Rx a la misma altura y radio
    que el Tx, incidencia barrida linealmente de psi_1 = psi_min a psi_2 = psi_max (helice):
      <S>     = |cos(dphi/2)| (cos psi_1 - cos psi_2) / (psi_2 - psi_1)
      q_media = Q(psi_1) + Q(psi_2),  Delta q = 2 [Q(psi_1) - Q(psi_2)],  Q(psi) = sqrt(n^2 - sin^2 psi)
      delta_z = c_w c / sqrt( (B q_media)^2 + (f_c Delta q)^2 )   (q lineal en psi: 12 var(q) = Delta q^2)"""
    p1, p2 = np.radians(psi_min), np.radians(psi_max)
    S = abs(np.cos(np.radians(dphi_deg) / 2)) * (np.cos(p1) - np.cos(p2)) / (p2 - p1)
    Q = lambda s: np.sqrt(n ** 2 - np.sin(s) ** 2)
    qm, dq = Q(p1) + Q(p2), 2 * (Q(p1) - Q(p2))
    lam0 = c / fc
    return {"delta_xy": X_J0 * lam0 / (2 * np.pi * S) if S > 1e-9 else np.inf,
            "delta_z": C_W[ventana] * c / np.sqrt((B * qm) ** 2 + (fc * dq) ** 2),
            "S_media": S, "q_media": qm, "Delta_q": dq}
