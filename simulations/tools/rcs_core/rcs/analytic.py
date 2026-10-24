"""Soluciones de referencia exactas."""
import numpy as np
from scipy.special import spherical_jn, spherical_yn


def sphere_pec_monostatic(radius, lam, nmax=None):
    """RCS monoestática exacta (serie de Mie) de una esfera PEC [m^2]."""
    x = 2 * np.pi * radius / lam
    nmax = nmax or int(x + 4 * x ** (1 / 3) + 10)
    n = np.arange(1, nmax + 1)
    jn, yn = spherical_jn(n, x), spherical_yn(n, x)
    jnp, ynp = spherical_jn(n, x, True), spherical_yn(n, x, True)
    hn, hnp = jn - 1j * yn, jnp - 1j * ynp
    a = jn / hn
    b = (jn + x * jnp) / (hn + x * hnp)
    s = np.sum((-1) ** n * (n + 0.5) * (b - a))
    return lam**2 / np.pi * np.abs(s) ** 2
