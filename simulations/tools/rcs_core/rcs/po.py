"""RCS monoestática de un conductor perfecto (PEC) por Óptica Física (PO).

Formulación
-----------
Corriente PO sobre caras iluminadas: J = 2 n x H_inc.  Para retrodispersión la
componente co-polar del campo lejano queda

    sqrt(sigma) ~ (2 sqrt(pi) / lambda) * sum_{caras iluminadas} (n . r) * I_T(2k r)

    I_T(w) = integral sobre el triángulo T de exp(j w . r') dS'

I_T se evalúa de forma exacta (fase lineal sobre un triángulo), así que las
caras NO necesitan ser menores que lambda: basta con representar bien la
geometría (un cubo son 12 triángulos).

Limitaciones de PO (importante):
  * Válida para objetos grandes frente a lambda (ka >~ 5-10).
  * No modela difracción de bordes/aristas, ondas de superficie ni reflexiones
    múltiples -> falla cerca de incidencia rasante y en la región resonante.
  * En monoestática PO no distingue polarizaciones (HH == VV, sin cross-pol).
  * Ocultamiento: solo "cara de espaldas" (n . r <= 0). Exacto para cuerpos
    convexos (cubo, cilindro, esfera, cono); para cóncavos usar SBR o MoM.
"""

import numpy as np
import numba as nb

C0 = 299_792_458.0


def directions(az_deg, el_deg):
    """Vectores unitarios de observación. az en plano xy desde +x, el desde el plano xy."""
    az, el = np.meshgrid(np.radians(az_deg), np.radians(el_deg), indexing="ij")
    return np.stack([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)], axis=-1)


def _phi1(t):
    """(exp(jt) - 1) / (jt), estable cerca de 0."""
    out = np.empty_like(t, dtype=complex)
    small = np.abs(t) < 1e-4
    ts = t[small]
    out[small] = 1 + 1j * ts / 2 - ts**2 / 6
    tl = t[~small]
    out[~small] = np.expm1(1j * tl) / (1j * tl)
    return out


def _tri_integral_normalized(x, y):
    """G(x,y) = int_0^1 int_0^{1-u} exp(j(xu + yv)) dv du  (triángulo de referencia).

    Simétrica en (x,y). Se usa la forma cerrada dividiendo por la mayor de |x|,|y|
    y una serie de Taylor cuando ambas son pequeñas.
    """
    swap = np.abs(x) > np.abs(y)
    x, y = np.where(swap, y, x), np.where(swap, x, y)  # ahora |y| >= |x|
    out = np.empty(x.shape, dtype=complex)
    small = np.abs(y) < 1e-3
    xs, ys = x[small], y[small]
    h1 = xs + ys
    h2 = xs**2 + xs * ys + ys**2
    h3 = xs**3 + xs**2 * ys + xs * ys**2 + ys**3
    out[small] = 0.5 + 1j * h1 / 6 - h2 / 24 - 1j * h3 / 120
    xl, yl = x[~small], y[~small]
    out[~small] = (np.exp(1j * yl) * _phi1(xl - yl) - _phi1(xl)) / (1j * yl)
    return out


@nb.njit(cache=True, inline="always")
def _phi1_s(t):
    if abs(t) < 1e-4:
        return 1 + 0.5j * t - t * t / 6
    return (np.exp(1j * t) - 1) / (1j * t)


@nb.njit(cache=True, inline="always")
def _g_s(x, y):
    """Versión escalar de _tri_integral_normalized."""
    if abs(x) > abs(y):
        x, y = y, x
    if abs(y) < 1e-3:
        return (0.5 + 1j * (x + y) / 6 - (x * x + x * y + y * y) / 24
                - 1j * (x**3 + x * x * y + x * y * y + y**3) / 120)
    return (np.exp(1j * y) * _phi1_s(x - y) - _phi1_s(x)) / (1j * y)


@nb.njit(parallel=True, cache=True, fastmath=False)
def _po_kernel(D, k, p0, p1, p2, n, area2, two_sided):
    nd, m = D.shape[0], p0.shape[0]
    out = np.empty(nd, dtype=np.complex128)
    for i in nb.prange(nd):
        dx, dy, dz = D[i, 0], D[i, 1], D[i, 2]
        acc = 0j
        for j in range(m):
            c = dx * n[j, 0] + dy * n[j, 1] + dz * n[j, 2]
            if two_sided:
                c = abs(c)
            if c <= 0.0:
                continue
            w = 2 * k
            a = w * (dx * p0[j, 0] + dy * p0[j, 1] + dz * p0[j, 2])
            b = w * (dx * p1[j, 0] + dy * p1[j, 1] + dz * p1[j, 2])
            cc = w * (dx * p2[j, 0] + dy * p2[j, 1] + dz * p2[j, 2])
            acc += c * area2[j] * np.exp(1j * a) * _g_s(b - a, cc - a)
        out[i] = acc
    return out


class POSolver:
    """Precalcula la geometría de la malla y evalúa sigma(az, el, f)."""

    def __init__(self, V, F, two_sided=False):
        V = np.asarray(V, float)
        F = np.asarray(F, int)
        self.p0, self.p1, self.p2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
        cr = np.cross(self.p1 - self.p0, self.p2 - self.p0)
        self.area2 = np.linalg.norm(cr, axis=1)  # 2 * área
        self.n = cr / self.area2[:, None]
        self.two_sided = two_sided  # para superficies abiertas (placas)

    def sqrt_sigma(self, rhat, freq):
        """Amplitud compleja sqrt(sigma) [m] para un bloque de direcciones rhat (Nd,3)."""
        lam = C0 / freq
        k = 2 * np.pi / lam
        w = 2 * k * rhat  # (Nd,3)
        cos_t = rhat @ self.n.T  # (Nd,M)
        if self.two_sided:
            lit = np.abs(cos_t)
        else:
            lit = np.where(cos_t > 0, cos_t, 0.0)
        a = w @ self.p0.T
        b = w @ self.p1.T
        c = w @ self.p2.T
        I = self.area2 * np.exp(1j * a) * _tri_integral_normalized(b - a, c - a)
        return (2 * np.sqrt(np.pi) / lam) * np.sum(lit * I, axis=1)

    def rcs(self, az_deg, el_deg, freq):
        """sigma [m^2] con forma (len(az), len(el)). Núcleo compilado (numba, multihilo)."""
        az_deg, el_deg = np.atleast_1d(az_deg), np.atleast_1d(el_deg)
        D = np.ascontiguousarray(directions(az_deg, el_deg).reshape(-1, 3))
        lam = C0 / freq
        s = _po_kernel(D, 2 * np.pi / lam, self.p0, self.p1, self.p2, self.n,
                       self.area2, self.two_sided)
        s *= 2 * np.sqrt(np.pi) / lam
        return (np.abs(s) ** 2).reshape(len(az_deg), len(el_deg))

    def rcs_freq(self, az_deg, el_deg, freqs):
        """sigma [m^2] con forma (len(freqs), len(az), len(el))."""
        return np.stack([self.rcs(az_deg, el_deg, f) for f in np.atleast_1d(freqs)])


def dbsm(sigma):
    return 10 * np.log10(np.maximum(sigma, 1e-30))
