"""Interfaz aire-suelo: Fresnel TM/TE con permitividad compleja, pseudo-Brewster,
atenuación y modelos de superficie rugosa.

Convención: eps_r = eps' - j eps''  (dependencia temporal e^{+j w t}); n = sqrt(eps_r)
con Re(n) > 0 e Im(n) <= 0.  La geometría de rayos dentro del suelo usa Re(n)
(aproximación de bajas pérdidas, tan_delta <= 0.2 en todos los suelos simulados).
"""
import numpy as np
from scipy.optimize import minimize_scalar


def eps_complex(e):
    """[eps', eps''] o escalar -> complejo eps' - j eps''."""
    if np.isscalar(e):
        return complex(e)
    return complex(e[0], -abs(e[1]))


def n_complex(eps):
    n = np.sqrt(complex(eps))
    return complex(abs(n.real), -abs(n.imag))


def _kz(eps, th_i):
    """n cos(theta_t) = sqrt(eps - sin^2 theta_i), rama con Re > 0 (onda que entra al suelo)."""
    s = np.sqrt(eps - np.sin(th_i) ** 2 + 0j)
    return np.where(s.real < 0, -s, s)


def fresnel_r(th_i, eps, pol):
    """Coeficiente de reflexión de amplitud aire->suelo (theta_i desde la normal, rad)."""
    th_i = np.asarray(th_i, float)
    c = np.cos(th_i)
    kz = _kz(eps, th_i)
    if pol.upper() == "TM":
        return (eps * c - kz) / (eps * c + kz)
    if pol.upper() == "TE":
        return (c - kz) / (c + kz)
    raise ValueError(pol)


def reflectance(th_i, eps, pol):
    return np.abs(fresnel_r(th_i, eps, pol)) ** 2


def transmittance(th_i, eps, pol):
    """Transmitancia de potencia de la interfaz (conservación en la interfaz: T = 1 - |r|^2)."""
    return 1.0 - reflectance(th_i, eps, pol)


def reflectance_mixed(th_i, eps, roll_deg):
    """Reflectancia vista por una antena nominalmente TM con error de atitud (roll) gamma:
    E = cos(g) TM + sin(g) TE; recepción co-polar con la misma antena girada.
    |cos^2 g r_TM + sin^2 g r_TE|^2 (despolarización por atitud)."""
    g = np.radians(roll_deg)
    return np.abs(np.cos(g) ** 2 * fresnel_r(th_i, eps, "TM") + np.sin(g) ** 2 * fresnel_r(th_i, eps, "TE")) ** 2


def pseudo_brewster_deg(eps):
    """Ángulo (aire) que minimiza |r_TM|^2 para eps complejo."""
    f = lambda t: reflectance(np.radians(t), eps, "TM")
    r = minimize_scalar(f, bounds=(1.0, 89.9), method="bounded", options={"xatol": 1e-6})
    return float(r.x), float(r.fun)


def attenuation_np_per_m(eps, k0):
    """Atenuación de campo del suelo [Np/m] = k0 |Im(n)|."""
    return k0 * abs(n_complex(eps).imag)


# ---------------------------------------------------------------------------
# Superficie rugosa
# ---------------------------------------------------------------------------
def coherent_factor(k0, sigma_h, th_i, th_s=None):
    """Factor de Rayleigh/Ament para la componente coherente (amplitud):
    exp(-0.5 k0^2 sigma_h^2 (cos th_i + cos th_s)^2)  (= exp(-2 k^2 s^2 cos^2) en especular)."""
    if th_s is None:
        th_s = th_i
    return np.exp(-0.5 * (k0 * sigma_h) ** 2 * (np.cos(th_i) + np.cos(th_s)) ** 2)


def spm_alpha(th_i, th_s, phi_s, eps, pol):
    """Coeficientes de polarización de primer orden del SPM (Rice; Ulaby-Moore-Fung vol. II;
    Tsang-Kong).  phi_s se mide desde la dirección especular directa (phi_s = 0 -> especular,
    phi_s = pi -> retrodispersión).  Reducen a alpha_hh = R_TE y a la expresión clásica de
    alpha_vv en retrodispersión.  En dirección especular alpha_vv = 0 exactamente en Brewster
    (eps real): es el 'Brewster biestático' del clutter difuso."""
    si, ss = np.sin(th_i), np.sin(th_s)
    ci, cs = np.cos(th_i), np.cos(th_s)
    qi = _kz(eps, th_i)
    qs = _kz(eps, th_s)
    if pol.upper() == "TE":
        return (eps - 1) * np.cos(phi_s) / ((ci + qi) * (cs + qs))
    return (eps - 1) * (eps * si * ss - np.cos(phi_s) * qi * qs) / ((eps * ci + qi) * (eps * cs + qs))


def spm_sigma0(th_i, th_s, phi_s, eps, pol, k0, sigma_h, corr_len):
    """Coeficiente de dispersión biestático sigma0_pp (adimensional, por unidad de área)
    de primer orden del SPM con espectro gaussiano:
        sigma0 = 8 k^4 s^2 cos^2(th_i) cos^2(th_s) |alpha|^2 W(Kx, Ky),
        W(K) = (l^2 / 4 pi) exp(-K^2 l^2 / 4),
        K = k (sin th_s cos phi_s - sin th_i, sin th_s sin phi_s).
    Retrodispersión: 8 k^4 s^2 cos^4 |alpha|^2 W(2k sin th)  (Ulaby et al., 1982)."""
    a = spm_alpha(th_i, th_s, phi_s, eps, pol)
    Kx = k0 * (np.sin(th_s) * np.cos(phi_s) - np.sin(th_i))
    Ky = k0 * np.sin(th_s) * np.sin(phi_s)
    W = corr_len ** 2 / (4 * np.pi) * np.exp(-(Kx ** 2 + Ky ** 2) * corr_len ** 2 / 4)
    return 8 * k0 ** 4 * sigma_h ** 2 * np.cos(th_i) ** 2 * np.cos(th_s) ** 2 * np.abs(a) ** 2 * W
