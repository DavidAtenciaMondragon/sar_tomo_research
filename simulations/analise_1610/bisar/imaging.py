"""Back-projection 3D con refraccion en interfaz plana (Fermat/Snell), en numba.

Datos de entrada por canal: muestras complejas s[m, k] en frecuencias f_k = f0 + k*df
(senal compleja tras compresion en rango, equivalente al rootData de GS_snell_test_script.m).
Imagen (BP exacto en fase, envolvente interpolada sobre perfil sobremuestreado):

    I(v) = sum_m w_m sum_k s[m,k] exp(+j 2 pi f_k tau_m(v))  /  (sum_m w_m * Nf)

tau_m(v) = (R1_Tx + n R2_Tx + R1_Rx + n R2_Rx)/c con R1 (aire) y R2 (suelo) del rayo de
Fermat; para voxeles con z >= 0 el camino es recto (aire).  La normalizacion por el numero
de muestras (media, no suma) hace que fusionar canales NO aumente la amplitud de la PSF por
el simple hecho de tener mas muestras (Orientacoes §5.1, §11).
"""
import numpy as np
import numba as nb


@nb.njit(cache=True, inline="always")
def _path(xr, yr, zr, xv, yv, zv, n):
    """Longitud optica (R1 + n R2) del rayo radar -> voxel.  Newton sobre u1 (Snell)."""
    dx = xr - xv
    dy = yr - yv
    u = np.sqrt(dx * dx + dy * dy)
    if zv >= 0.0:
        dz = zr - zv
        return np.sqrt(u * u + dz * dz)
    h = zr
    d = -zv
    if u < 1e-9:
        return h + n * d
    # semilla: angulo recto hasta el voxel, refractado
    s_i = u / np.sqrt(u * u + (h + d) * (h + d))
    s_t = s_i / n
    u2 = d * s_t / np.sqrt(1.0 - s_t * s_t)
    u1 = u - u2
    if u1 < 0.0:
        u1 = 0.5 * u
    for _ in range(8):
        R1 = np.sqrt(u1 * u1 + h * h)
        u2 = u - u1
        R2 = np.sqrt(u2 * u2 + d * d)
        f = u1 / R1 - n * u2 / R2
        fp = h * h / (R1 * R1 * R1) + n * d * d / (R2 * R2 * R2)
        step = f / fp
        u1 -= step
        if u1 < 0.0:
            u1 = 0.0
        if u1 > u:
            u1 = u
        if abs(step) < 1e-10:
            break
    R1 = np.sqrt(u1 * u1 + h * h)
    u2 = u - u1
    R2 = np.sqrt(u2 * u2 + d * d)
    return R1 + n * R2


@nb.njit(parallel=True, cache=True)
def _bp_kernel(tx, rx, G, w, vox, n, f0, df, nfft, c):
    V = vox.shape[0]
    M = tx.shape[0]
    out = np.zeros(V, dtype=np.complex128)
    T = 1.0 / df
    wsum = 0.0
    for m in range(M):
        wsum += w[m]
    for iv in nb.prange(V):
        xv, yv, zv = vox[iv, 0], vox[iv, 1], vox[iv, 2]
        acc = 0j
        for m in range(M):
            if w[m] == 0.0:
                continue
            L = _path(tx[m, 0], tx[m, 1], tx[m, 2], xv, yv, zv, n) + \
                _path(rx[m, 0], rx[m, 1], rx[m, 2], xv, yv, zv, n)
            tau = L / c
            q = (tau % T) * nfft * df
            i0 = int(q)
            fr = q - i0
            i0 = i0 % nfft
            i1 = (i0 + 1) % nfft
            g = (1.0 - fr) * G[m, i0] + fr * G[m, i1]
            ph = 2.0 * np.pi * f0 * tau
            acc += w[m] * g * (np.cos(ph) + 1j * np.sin(ph))
        out[iv] = acc / wsum
    return out


def spectral_window(Nf, kind="rect"):
    """Ponderacion espectral en rango normalizada a media 1 (conserva la amplitud del pico)."""
    if kind in (None, "rect"):
        return np.ones(Nf)
    if kind == "hann":
        w = np.hanning(Nf + 2)[1:-1]
    elif kind == "hamming":
        w = np.hamming(Nf)
    else:
        raise ValueError(kind)
    return w / w.mean()


def range_profiles(s, upsample=32, window="rect"):
    """G_m(tau_q) = (1/Nf) sum_k w_k s[m,k] exp(j 2 pi k df tau_q), tau_q = q/(Nfft df)."""
    M, Nf = s.shape
    s = s * spectral_window(Nf, window)[None, :]
    nfft = Nf * upsample
    G = np.fft.ifft(s, n=nfft, axis=1) * (nfft / Nf)
    return np.ascontiguousarray(G), nfft


def apply_gate(G, nfft, df, lo, hi):
    """Ventana de recepcion en rango: anula el perfil fuera de [lo_m, hi_m] (s) en cada pulso
    (perfil periodico de periodo 1/df)."""
    T = 1.0 / df
    tq = np.arange(nfft) / (nfft * df)
    inwin = ((tq[None, :] - (np.asarray(lo)[:, None] % T)) % T) <= (np.asarray(hi) - np.asarray(lo))[:, None]
    return np.ascontiguousarray(np.where(inwin, G, 0.0))


def backproject(tx, rx, s, freqs, vox, n_real, c, weights=None, upsample=32, gate=None, window="rect"):
    """Imagen compleja en los voxeles vox (V,3). tx, rx: (M,3); s: (M,Nf); freqs uniformes.
    gate: None o (tau_lo (M,), tau_hi (M,)) ventana de recepcion en segundos."""
    freqs = np.asarray(freqs, float)
    df = freqs[1] - freqs[0]
    G, nfft = range_profiles(s, upsample, window)
    if gate is not None:
        G = apply_gate(G, nfft, df, gate[0], gate[1])
    if weights is None:
        weights = np.ones(tx.shape[0])
    return _bp_kernel(np.ascontiguousarray(tx, float), np.ascontiguousarray(rx, float), np.ascontiguousarray(G),
                      np.ascontiguousarray(weights, float), np.ascontiguousarray(vox, float),
                      float(n_real), float(freqs[0]), float(df), int(nfft), float(c))


class Grid:
    """Rejilla cartesiana regular (ejes x, y, z) con conversion a lista de voxeles."""

    def __init__(self, x, y, z):
        self.x, self.y, self.z = np.asarray(x, float), np.asarray(y, float), np.asarray(z, float)
        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing="ij")
        self.shape = X.shape
        self.vox = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)

    @classmethod
    def around(cls, c, hx, dx, hz_lo, hz_hi, dz, hy=None):
        hy = hx if hy is None else hy
        x = c[0] + np.arange(-hx, hx + dx / 2, dx)
        y = c[1] + np.arange(-hy, hy + dx / 2, dx)
        z = c[2] + np.arange(hz_lo, hz_hi + dz / 2, dz)
        return cls(x, y, z)

    def cube(self, I):
        return np.asarray(I).reshape(self.shape)

    @property
    def dv(self):
        d = lambda a: (a[1] - a[0]) if len(a) > 1 else 1.0
        return d(self.x) * d(self.y) * d(self.z)
