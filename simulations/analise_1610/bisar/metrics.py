"""Metricas obligatorias (Orientacoes §6): PSF, lobulos, imagen, informacion.

Definiciones (todas sobre |I| con la imagen normalizada por el numero de muestras):
  FWHM_x,y,z   ancho a -3 dB en cortes por el pico (interpolacion lineal en dB).
  V_-3dB       volumen de la componente conexa con |I|^2 >= pico/2 que contiene el pico.
  anisotropia  sqrt(lambda_max/lambda_min) del tensor de segundo momento del lobulo -3 dB
               (3D) y en el plano horizontal (xy).
  lobulo ppal  elipsoide con semiejes = distancia al primer nulo en x, y, z.
  PSLR         max |I|^2 fuera del lobulo principal / pico  [dB]  (ventana de la rejilla)
  ISLR         energia fuera / energia dentro del lobulo principal [dB] (misma ventana)
  contraste    pico / potencia media fuera del lobulo principal [dB]
  loc. error   |r_pico - r_verdadero| con refinamiento parabolico en cada eje
"""
import numpy as np
from scipy import ndimage

DB = lambda x: 10 * np.log10(np.maximum(x, 1e-300))


# ----------------------------------------------------------------------- 1D
def _cross(x0, x1, y0, y1, thr):
    return x0 + (thr - y0) * (x1 - x0) / (y1 - y0) if y1 != y0 else x0


def fwhm_1d(axis, amp, ipk=None, thr_db=-3.0, return_censored=False):
    a = np.asarray(amp, float)
    ydb = 20 * np.log10(np.maximum(a / a.max(), 1e-300))
    ipk = int(np.argmax(a)) if ipk is None else ipk
    cens = False
    L = np.where(ydb[:ipk] < thr_db)[0]
    if len(L) == 0:
        xl, cens = axis[0], True
    else:
        i = L[-1]
        xl = _cross(axis[i], axis[i + 1], ydb[i], ydb[i + 1], thr_db)
    R = np.where(ydb[ipk:] < thr_db)[0]
    if len(R) == 0:
        xr, cens = axis[-1], True
    else:
        i = ipk + R[0]
        xr = _cross(axis[i - 1], axis[i], ydb[i - 1], ydb[i], thr_db)
    w = abs(xr - xl)
    return (w, cens) if return_censored else (w if not cens else np.nan)


def first_nulls_1d(axis, amp, ipk=None):
    """Distancias (izq, der) desde el pico al primer minimo local."""
    a = np.asarray(amp, float)
    ipk = int(np.argmax(a)) if ipk is None else ipk
    i = ipk
    while i > 0 and a[i - 1] < a[i]:
        i -= 1
    j = ipk
    while j < len(a) - 1 and a[j + 1] < a[j]:
        j += 1
    return axis[ipk] - axis[i], axis[j] - axis[ipk], i, j


def pslr_1d(axis, amp):
    a = np.asarray(amp, float) ** 2
    ipk = int(np.argmax(a))
    _, _, i, j = first_nulls_1d(axis, np.sqrt(a), ipk)
    side = np.r_[a[:i], a[j + 1:]]
    return DB(side.max() / a[ipk]) if len(side) else np.nan


def _parabolic(y, i):
    if i <= 0 or i >= len(y) - 1:
        return 0.0
    d = y[i - 1] - 2 * y[i] + y[i + 1]
    return 0.5 * (y[i - 1] - y[i + 1]) / d if d != 0 else 0.0


# ----------------------------------------------------------------------- 3D
def psf_metrics(grid, I, true_pos=None):
    A = np.abs(grid.cube(I))
    P = A ** 2
    ipk = np.unravel_index(np.argmax(A), A.shape)
    ix, iy, iz = ipk
    x, y, z = grid.x, grid.y, grid.z
    pk = P[ipk]
    out = {}
    # pico refinado y error de localizacion
    rx = x[ix] + _parabolic(A[:, iy, iz], ix) * (x[1] - x[0])
    ry = y[iy] + _parabolic(A[ix, :, iz], iy) * (y[1] - y[0])
    rz = z[iz] + _parabolic(A[ix, iy, :], iz) * (z[1] - z[0])
    out["peak_amp"] = A[ipk]
    out["peak_x"], out["peak_y"], out["peak_z"] = rx, ry, rz
    if true_pos is not None:
        out["loc_err_m"] = float(np.linalg.norm(np.array([rx, ry, rz]) - np.asarray(true_pos)))
    # FWHM por ejes
    for lab, ax, prof, ii in [("x", x, A[:, iy, iz], ix), ("y", y, A[ix, :, iz], iy), ("z", z, A[ix, iy, :], iz)]:
        w, c = fwhm_1d(ax, prof, ii, return_censored=True)
        out[f"fwhm_{lab}"] = w
        out[f"fwhm_{lab}_cens"] = int(c)
        out[f"pslr1d_{lab}_dB"] = pslr_1d(ax, prof)
    # lobulo -3 dB conexo: volumen y anisotropia
    lab, _ = ndimage.label(P >= pk / 2)
    ml = lab == lab[ipk]
    out["vol3dB_m3"] = ml.sum() * grid.dv
    out["vol3dB_cens"] = int(ml[0].any() or ml[-1].any() or ml[:, 0].any() or ml[:, -1].any()
                             or ml[:, :, 0].any() or ml[:, :, -1].any())
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    w = P * ml
    pts = np.stack([X[ml], Y[ml], Z[ml]], 1)
    ww = w[ml] / w[ml].sum()
    mu = (pts * ww[:, None]).sum(0)
    C = ((pts - mu).T * ww) @ (pts - mu)
    ev = np.sort(np.linalg.eigvalsh(C))
    out["aniso3d"] = float(np.sqrt(ev[-1] / max(ev[0], 1e-30)))
    evh = np.sort(np.linalg.eigvalsh(C[:2, :2]))
    out["aniso_xy"] = float(np.sqrt(evh[-1] / max(evh[0], 1e-30)))
    out["aniso_z_xy"] = float(np.sqrt(C[2, 2] / max(0.5 * (C[0, 0] + C[1, 1]), 1e-30)))
    # lobulo principal elipsoidal por primeros nulos
    sem = []
    for ax, prof, ii in [(x, A[:, iy, iz], ix), (y, A[ix, :, iz], iy), (z, A[ix, iy, :], iz)]:
        l, r, _, _ = first_nulls_1d(ax, prof, ii)
        sem.append(max(min(l, r), ax[1] - ax[0]))
    out["null_x"], out["null_y"], out["null_z"] = sem
    E = ((X - x[ix]) / sem[0]) ** 2 + ((Y - y[iy]) / sem[1]) ** 2 + ((Z - z[iz]) / sem[2]) ** 2 <= 1.0
    side = P[~E]
    out["PSLR_dB"] = DB(side.max() / pk) if side.size else np.nan
    out["ISLR_dB"] = DB(side.sum() / P[E].sum()) if side.size else np.nan
    out["contrast_dB"] = DB(pk / side.mean()) if side.size else np.nan
    return out


# --------------------------------------------------------------- k-espacio
def kvectors(sysm, info_t, info_r):
    """K = (2 pi f n / c)(u_t + u_r) para cada pulso y frecuencia (gradiente de fase)."""
    u = info_t["u_soil"] + info_r["u_soil"]
    ks = 2 * np.pi * sysm.freqs * sysm.n_real / sysm.c
    return (u[:, None, :] * ks[None, :, None]).reshape(-1, 3)


def kspace_metrics(K, w=None, dk=0.25):
    """Soporte en k: extension, ocupacion (celdas de dk^3) y dispersion (Fisher)."""
    w = np.ones(len(K)) if w is None else w
    kh = np.hypot(K[:, 0], K[:, 1])
    occ = np.unique(np.floor(K / dk).astype(np.int64), axis=0).shape[0]
    occ_h = np.unique(np.floor(K[:, :2] / dk).astype(np.int64), axis=0).shape[0]
    ww = w / w.sum()
    mu = (K * ww[:, None]).sum(0)
    C = ((K - mu).T * ww) @ (K - mu)
    return {"kh_min": kh.min(), "kh_max": kh.max(), "kz_min": K[:, 2].min(), "kz_max": K[:, 2].max(),
            "k_occ_vol": occ * dk ** 3, "k_occ_area_h": occ_h * dk ** 2,
            "sig_kx": np.sqrt(C[0, 0]), "sig_ky": np.sqrt(C[1, 1]), "sig_kz": np.sqrt(C[2, 2])}


def conditioning(Klist, wlist, hx=0.3, nx=7, hz=1.0, nz=5, sub=8):
    """Matriz de adquisicion linealizada A[i, v] = exp(-j K_i . (r_v - r0)) sobre una rejilla
    local nx x nx x nz (semiejes hx, hz).  Filas de cada canal escaladas por sqrt(w/N_canal):
    Gram normalizado (invariante al numero de muestras).  Devuelve numero de condicion,
    rango efectivo (Roy-Vetterli) y coherencia mutua."""
    g = np.linspace(-hx, hx, nx)
    gz = np.linspace(-hz, hz, nz)
    X, Y, Z = np.meshgrid(g, g, gz, indexing="ij")
    R = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
    rows = []
    for K, w in zip(Klist, wlist):
        Ks = K[::sub]
        rows.append(np.sqrt(w / len(Ks)) * np.exp(-1j * Ks @ R.T))
    A = np.vstack(rows)
    s = np.linalg.svd(A, compute_uv=False)
    p = s ** 2 / (s ** 2).sum()
    erank = float(np.exp(-(p * np.log(np.maximum(p, 1e-300))).sum()))
    An = A / np.linalg.norm(A, axis=0, keepdims=True)
    Gm = np.abs(An.conj().T @ An)
    np.fill_diagonal(Gm, 0)
    return {"cond": float(s[0] / s[-1]), "cond_dB": float(20 * np.log10(s[0] / s[-1])),
            "erank": erank, "erank_frac": erank / R.shape[0], "coh_max": float(Gm.max()),
            "coh_mean": float(Gm.sum() / (Gm.size - Gm.shape[0]))}


# ------------------------------------------------------------ multi-blanco
def separability(grid, I, targets, tol_h=0.15, tol_z=0.5, dip_db=3.0):
    """Para cada blanco verdadero: maximo local dentro de la tolerancia; 'resuelto' si el
    perfil recto hacia cada vecino mas cercano cae >= dip_db bajo el menor de los dos picos."""
    A = np.abs(grid.cube(I))
    x, y, z = grid.x, grid.y, grid.z
    from scipy.interpolate import RegularGridInterpolator
    f = RegularGridInterpolator((x, y, z), A, bounds_error=False, fill_value=0.0)
    T = np.asarray(targets)
    found, err, peaks = [], [], []
    for t in T:
        mx = (np.abs(x - t[0]) <= tol_h)
        my = (np.abs(y - t[1]) <= tol_h)
        mz = (np.abs(z - t[2]) <= tol_z)
        sub = A[np.ix_(mx, my, mz)]
        j = np.unravel_index(np.argmax(sub), sub.shape)
        p = np.array([x[mx][j[0]], y[my][j[1]], z[mz][j[2]]])
        peaks.append(p)
        err.append(np.linalg.norm(p - t))
    peaks = np.array(peaks)
    amp = f(peaks)
    D = np.linalg.norm(T[:, None] - T[None], axis=2)
    np.fill_diagonal(D, np.inf)
    res = []
    for i in range(len(T)):
        nb = np.where(D[i] <= D[i].min() * 1.01)[0]
        ok = True
        for j in nb:
            seg = peaks[i][None] + np.linspace(0, 1, 41)[:, None] * (peaks[j] - peaks[i])[None]
            prof = f(seg)
            ok &= 20 * np.log10(prof.min() / min(amp[i], amp[j]) + 1e-300) <= -dip_db
        res.append(ok)
    res = np.array(res)
    return {"frac_resolved": float(res.mean()), "loc_err_mean": float(np.mean(err)),
            "amp_spread_dB": float(20 * np.log10(amp.max() / max(amp.min(), 1e-300))),
            "resolved": res, "peaks": peaks, "amps": amp}
