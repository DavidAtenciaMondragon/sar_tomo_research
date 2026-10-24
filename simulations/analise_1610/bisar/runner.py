"""Utilidades compartidas por los experimentos: canales, imagenes mono/bi/fusion y metricas.

Reglas de normalizacion (declaradas en el relatorio, Orientacoes §5.1 y §11):
  * cada imagen es la MEDIA coherente sobre sus muestras (no la suma);
  * 'fus'  = (I_mono + I_bi)/2 : fusion coherente con 2N muestras (N de cada canal);
  * 'fusN' = mono en pulsos impares + bi en pulsos pares (N muestras en total): mismo
    presupuesto de muestras que la imagen monoestatica -> aisla el efecto geometrico;
  * la SNR se calcula con la varianza del ruido termico propagada por el BP (M_eff pulsos),
    de modo que el aumento de SNR por duplicar muestras queda separado y explicito.
"""
import numpy as np

from . import geometry, imaging, metrics, scene


def channels(p, dphi_deg, kappa=1.0, Ptx=None):
    if Ptx is None:
        _, Ptx = geometry.tx_trajectory(p)
    Prx = Ptx.copy() if np.all(np.asarray(dphi_deg) == 0) and np.all(np.asarray(kappa) == 1) else \
        geometry.rx_from_tx(Ptx, dphi_deg, kappa)
    return Ptx, Prx


def gate_for(sysm, Ptx, Prx, lo_m=-30.0, hi_m=60.0):
    """Ventana de recepcion: camino Tx -> centro de la escena (0,0,0) -> Rx  + [lo_m, hi_m] m.
    Excluye el reflejo de nadir (~150 m antes) y conserva blancos y superficie cercana."""
    L = np.linalg.norm(Ptx, axis=1) + np.linalg.norm(Prx, axis=1)
    return (L + lo_m) / sysm.c, (L + hi_m) / sysm.c


def bp(sysm, Ptx, Prx, s, vox, w=None, n_bp=None, upsample=32, gate=None):
    """gate: None, True (ventana por defecto de gate_for) o tupla (tau_lo, tau_hi)."""
    n = sysm.n_real if n_bp is None else n_bp
    if gate is True:
        gate = gate_for(sysm, Ptx, Prx, *sysm.p.get("ventana_rango_m", [-30.0, 60.0]))
    return imaging.backproject(Ptx, Prx, s, sysm.freqs, vox, n, sysm.c, w, upsample, gate,
                               sysm.p.get("ventana_espectral", "rect"))


def four_images(sysm, Ptx, Prx, s_m, s_b, vox, I_mono=None, I_mono_odd=None, equal_budget=True):
    """Devuelve dict con imagenes complejas mono, bi, fus (2N) y fusN (N)."""
    M = Ptx.shape[0]
    odd = (np.arange(M) % 2 == 1).astype(float)
    even = 1.0 - odd
    out = {}
    out["mono"] = bp(sysm, Ptx, Ptx, s_m, vox) if I_mono is None else I_mono
    out["bi"] = bp(sysm, Ptx, Prx, s_b, vox)
    out["fus"] = 0.5 * (out["mono"] + out["bi"])
    if equal_budget:
        Im_odd = bp(sysm, Ptx, Ptx, s_m, vox, odd) if I_mono_odd is None else I_mono_odd
        Ib_even = bp(sysm, Ptx, Prx, s_b, vox, even)
        out["fusN"] = 0.5 * (Im_odd + Ib_even)
        out["_mono_odd"] = Im_odd
    return out


M_EFF = {"mono": 1.0, "bi": 1.0, "fus": 2.0, "fusN": 1.0}


def snr_db(sysm, I_peak, M, ch):
    return 10 * np.log10(np.abs(I_peak) ** 2 / sysm.noise_var_image(M * M_EFF[ch]))


def psf_rows(sysm, grid, imgs, tgt, M, extra=None, ref_amp=None):
    """Una fila de metricas por canal.  ref_amp: amplitud ideal (media de |a|) para el error de amplitud."""
    rows = []
    for ch in ["mono", "bi", "fus", "fusN"]:
        if ch not in imgs:
            continue
        m = metrics.psf_metrics(grid, imgs[ch], tgt)
        m["canal"] = ch
        m["SNR_dB"] = snr_db(sysm, m["peak_amp"], M, ch)
        if ref_amp is not None:
            m["amp_err_dB"] = 20 * np.log10(m["peak_amp"] / ref_amp[ch])
        if extra:
            m.update(extra)
        rows.append(m)
    return rows


def cut_widths(sysm, Ptx, Prx, s, tgt, h=25.0, dx=0.01, hz=4.0, dz=0.01):
    """FWHM horizontal en cortes x, y largos (para casos censurados, p.ej. Ref1) y z."""
    out = {}
    ax = np.arange(-h, h + dx / 2, dx)
    for lab, vox in [("x", np.c_[tgt[0] + ax, tgt[1] + 0 * ax, tgt[2] + 0 * ax]),
                     ("y", np.c_[tgt[0] + 0 * ax, tgt[1] + ax, tgt[2] + 0 * ax])]:
        I = np.abs(bp(sysm, Ptx, Prx, s, vox))
        w, c = metrics.fwhm_1d(ax, I, return_censored=True)
        out[f"fwhm_{lab}_largo"] = w if not c else np.nan
    return out
