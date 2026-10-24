"""Sintesis de ecos complejos por canal (mismo pulso para mono y bi).

Canal = (posiciones Tx (M,3), posiciones Rx (M,3)).  El canal monoestatico usa Rx = Tx
(el transceiver), el biestatico el Rx remoto.  Para cada pulso m y frecuencia f_k:

  blanco:  s = A_k sqrt(Gt Gr) tau_t tau_r e^{-alpha_k (R2t+R2r)} S(m,k) e^{-j 2 pi f_k tau} / (Rt Rr)
           A_k = sqrt(Pt) lambda_k / (4 pi)^{3/2};  S = sqrt(sigma) complejo (PO) o 1 (punto 1 m^2)
  superficie especular (teoria de imagenes, fase estacionaria exacta para plano infinito):
           s = sqrt(Pt Gt Gr) lambda_k/(4 pi) * R_eff(theta) chi_coh / L * e^{-j 2 pi f_k L/c}
  superficie difusa (facetas, SPM 1er orden, superficie 'congelada' en toda la apertura):
           s = A_k sqrt(Gt Gr) a_p sqrt(sigma0 dA) e^{-j 2 pi f_k tau_p} / (Rt Rr),  a_p ~ CN(0,1)

Ganancias: haz gaussiano (anchos a -3 dB AperturaAzimut x AperturaElev de radarTx*.json)
apuntado al centro de la escena.  Convencion del radar: Pr = |s|^2 [W].
"""
import numpy as np
import numba as nb

from . import em

C_DEFAULT = 3e8


# ------------------------------------------------------------------ rayos al blanco
def refract(P, tgt, n):
    """Rayo de Fermat radar (M,3) -> blanco enterrado tgt (3,), interfaz z=0, indice n real.
    Devuelve dict con R1 (aire), R2 (suelo), Q (punto de refraccion), th_i, th_t (rad),
    u_soil (unit. blanco -> Q, direccion de observacion en el suelo), u_air (unit. Q -> radar)."""
    P = np.atleast_2d(P).astype(float)
    tgt = np.asarray(tgt, float)
    d = -tgt[2]
    h = P[:, 2]
    dxy = P[:, :2] - tgt[:2]
    u = np.hypot(dxy[:, 0], dxy[:, 1])
    s_i = u / np.sqrt(u ** 2 + (h + d) ** 2)
    s_t = s_i / n
    u1 = u - d * s_t / np.sqrt(1 - s_t ** 2)
    u1 = np.clip(u1, 0, u)
    for _ in range(30):
        R1 = np.sqrt(u1 ** 2 + h ** 2)
        u2 = u - u1
        R2 = np.sqrt(u2 ** 2 + d ** 2)
        f = u1 / R1 - n * u2 / R2
        fp = h ** 2 / R1 ** 3 + n * d ** 2 / R2 ** 3
        st = f / fp
        u1 = np.clip(u1 - st, 0, u)
        if np.max(np.abs(st)) < 1e-11:
            break
    R1 = np.sqrt(u1 ** 2 + h ** 2)
    u2 = u - u1
    R2 = np.sqrt(u2 ** 2 + d ** 2)
    e = np.where(u[:, None] > 0, dxy / np.maximum(u, 1e-15)[:, None], 0.0)
    Q = np.c_[tgt[0] + e[:, 0] * u2, tgt[1] + e[:, 1] * u2, np.zeros_like(u)]
    u_soil = (Q - tgt) / R2[:, None]
    u_air = (P - Q) / R1[:, None]
    return dict(R1=R1, R2=R2, Q=Q, th_i=np.arccos(np.clip(h / R1, -1, 1)),
                th_t=np.arccos(np.clip(d / R2, -1, 1)), u_soil=u_soil, u_air=u_air)


# ------------------------------------------------------------------------- antena
def antenna_frames(P, aim=(0.0, 0.0, 0.0)):
    """Base del haz por pulso: b (apunte), h (horizontal), e (elevacion)."""
    b = np.asarray(aim, float) - P
    b /= np.linalg.norm(b, axis=1, keepdims=True)
    hz = np.cross(np.array([0.0, 0.0, 1.0]), b)
    hz /= np.linalg.norm(hz, axis=1, keepdims=True)
    e = np.cross(b, hz)
    return b, hz, e


def gain(P, dirs, bw_az_deg, bw_el_deg, G0_dBi, aim=(0.0, 0.0, 0.0)):
    """Ganancia de potencia (lineal) en las direcciones dirs (M,3) (desde la antena)."""
    b, hz, e = antenna_frames(P, aim)
    cb = np.einsum("ij,ij->i", dirs, b)
    daz = np.arctan2(np.einsum("ij,ij->i", dirs, hz), cb)
    delv = np.arctan2(np.einsum("ij,ij->i", dirs, e), cb)
    a, el = np.radians(bw_az_deg), np.radians(bw_el_deg)
    return 10 ** (G0_dBi / 10) * np.exp(-4 * np.log(2) * ((daz / a) ** 2 + (delv / el) ** 2))


# ----------------------------------------------------------------------- sistema
class System:
    """Parametros fisicos comunes (del JSON) + suelo + polarizacion."""

    def __init__(self, p, soil=None, pol=None, roll_deg=None):
        self.p = p
        self.c = p["c"]
        self.fc = p["fc"]
        self.B = p["B"]
        self.Nf = p["Nf"]
        self.freqs = self.fc + np.linspace(-self.B / 2, self.B / 2, self.Nf)
        self.Pt = p["Pt"]
        self.G0 = p["antena"]["G_dBi"]
        self.bw_az = p["antena"]["bw_az_deg"]
        self.bw_el = p["antena"]["bw_el_deg"]
        self.aim = tuple(p["antena"].get("apunte", [0, 0, 0]))
        self.eps = em.eps_complex(soil if soil is not None else p["suelo"]["eps_r"])
        self.n = em.n_complex(self.eps)
        self.n_real = self.n.real
        self.pol = (pol or p.get("pol", "TM")).upper()
        self.roll = p.get("roll_deg", 0.0) if roll_deg is None else roll_deg
        self.Tp = p["Tp"]
        self.kT = 1.380649e-23 * p["T_sys"] * 10 ** (p.get("NF_dB", 0.0) / 10)

    # pesos de mezcla de polarizacion por error de atitud (roll) igual en Tx y Rx
    def pol_weights(self):
        g = np.radians(self.roll)
        main, cross = np.cos(g) ** 2, np.sin(g) ** 2
        return (main, cross) if self.pol in ("TM", "TMAL") else (cross, main)   # (peso TM, peso TE)

    def amp_const(self):
        lam = self.c / self.freqs
        return np.sqrt(self.Pt) * lam / (4 * np.pi) ** 1.5           # (Nf,)

    def noise_var_image(self, M_eff):
        """Varianza del ruido termico en un voxel de la imagen normalizada (media):
        filtro adaptado por pulso (energia Pr Tp / kT) e integracion coherente de M_eff pulsos."""
        return self.kT / (self.Tp * M_eff)


# ---------------------------------------------------------------- eco de blancos
def target_echo(sysm, Ptx, Prx, tgt, scatter=None, R_override=None):
    """Eco (M,Nf) de un blanco enterrado centrado en tgt.
    scatter: None -> punto isotropo con sigma = 1 m^2 (S = 1);
             objeto rcs.PO -> S(m,k) por PO con las direcciones refractadas.
    Devuelve (s, info)."""
    gt = refract(Ptx, tgt, sysm.n_real)
    gr = refract(Prx, tgt, sysm.n_real)
    f = sysm.freqs
    k0 = 2 * np.pi * f / sysm.c
    alpha = k0 * abs(sysm.n.imag)                                    # (Nf,)
    wTM, wTE = sysm.pol_weights()
    # transmitancias de amplitud (sqrt de potencia) en bajada y subida (reciprocidad)
    tTM = np.sqrt(em.transmittance(gt["th_i"], sysm.eps, "TM") * em.transmittance(gr["th_i"], sysm.eps, "TM"))
    tTE = np.sqrt(em.transmittance(gt["th_i"], sysm.eps, "TE") * em.transmittance(gr["th_i"], sysm.eps, "TE"))
    M = Ptx.shape[0]
    if scatter is None:
        S_TM = np.ones((M, len(f)))
        S_TE = np.ones((M, len(f)))
    else:
        ks = 2 * np.pi * f * sysm.n_real / sysm.c
        I = np.repeat(gt["u_soil"], len(f), axis=0)
        Sd = np.repeat(gr["u_soil"], len(f), axis=0)
        K = np.tile(ks, M)
        S_TM = scatter.amplitude(I, Sd, K, "TM").reshape(M, len(f)) if wTM > 0 else 0.0
        S_TE = scatter.amplitude(I, Sd, K, "TE").reshape(M, len(f)) if wTE > 0 else 0.0
    if sysm.pol == "TMAL":
        # Tx TM; Rx de doble polarizacion combinado para alinearse con el campo del Tx
        from .rcs import aligned_rx
        _, _, cpsi = aligned_rx(gt["u_soil"], gr["u_soil"])
        TdTM = em.transmittance(gt["th_i"], sysm.eps, "TM")
        Tup = cpsi ** 2 * em.transmittance(gr["th_i"], sysm.eps, "TM") + \
            (1 - cpsi ** 2) * em.transmittance(gr["th_i"], sysm.eps, "TE")
        tAL = np.sqrt(TdTM * Tup)
        S_AL = np.ones((M, len(f))) if scatter is None else \
            scatter.amplitude(I, Sd, K, "TMAL").reshape(M, len(f))
        Seff = tAL[:, None] * S_AL
    else:
        Seff = wTM * tTM[:, None] * S_TM + wTE * tTE[:, None] * S_TE
    G = np.sqrt(gain(Ptx, -gt["u_air"], sysm.bw_az, sysm.bw_el, sysm.G0, sysm.aim) *
                gain(Prx, -gr["u_air"], sysm.bw_az, sysm.bw_el, sysm.G0, sysm.aim))
    Rt = gt["R1"] + gt["R2"]
    Rr = gr["R1"] + gr["R2"]
    tau = (gt["R1"] + sysm.n_real * gt["R2"] + gr["R1"] + sysm.n_real * gr["R2"]) / sysm.c
    loss = np.exp(-np.outer(gt["R2"] + gr["R2"], alpha))
    s = (sysm.amp_const()[None, :] * (G / (Rt * Rr))[:, None] * loss * Seff *
         np.exp(-2j * np.pi * np.outer(tau, f)))
    info = dict(gt=gt, gr=gr, tau=tau, G=G, Rt=Rt, Rr=Rr, tTM=tTM, tTE=tTE)
    return s, info


def multi_target_echo(sysm, Ptx, Prx, targets):
    """targets: lista de (pos, scatter|None, peso complejo)."""
    s = 0
    for pos, sc, wgt in targets:
        si, _ = target_echo(sysm, Ptx, Prx, pos, sc)
        s = s + wgt * si
    return s


# ------------------------------------------------------- superficie: especular
def specular_echo(sysm, Ptx, Prx, sigma_h=0.0):
    """Reflexion especular coherente del plano z=0 (punto especular por pulso)."""
    img = Ptx * np.array([1, 1, -1.0])
    D = Prx - img
    L = np.linalg.norm(D, axis=1)
    tq = Ptx[:, 2] / (Ptx[:, 2] + Prx[:, 2])
    Q = img + tq[:, None] * D
    th = np.arccos(np.clip((Ptx[:, 2] + Prx[:, 2]) / L, -1, 1))
    wTM, wTE = sysm.pol_weights()
    Reff = wTM * em.fresnel_r(th, sysm.eps, "TM") + wTE * em.fresnel_r(th, sysm.eps, "TE")
    k0 = 2 * np.pi * sysm.fc / sysm.c
    chi = em.coherent_factor(k0, sigma_h, th)
    dT = Q - Ptx
    dT /= np.linalg.norm(dT, axis=1, keepdims=True)
    dR = Q - Prx
    dR /= np.linalg.norm(dR, axis=1, keepdims=True)
    G = np.sqrt(gain(Ptx, dT, sysm.bw_az, sysm.bw_el, sysm.G0, sysm.aim) *
                gain(Prx, dR, sysm.bw_az, sysm.bw_el, sysm.G0, sysm.aim))
    lam = sysm.c / sysm.freqs
    a = np.sqrt(sysm.Pt) * lam / (4 * np.pi)
    s = (a[None, :] * (G * Reff * chi / L)[:, None] * np.exp(-2j * np.pi * np.outer(L / sysm.c, sysm.freqs)))
    return s, dict(Q=Q, theta=th, L=L, G=G, Reff=Reff, chi=chi)


# ---------------------------------------------------------- superficie: difusa
@nb.njit(cache=True, inline="always")
def _kz(eps, s2):
    q = np.sqrt(eps - s2)
    if q.real < 0:
        q = -q
    return q


@nb.njit(parallel=True, cache=True)
def _diffuse_kernel(tx, rx, bt, ht, et, br, hr, er, fac, amp, eps, k0, sh, cl, wTM, wTE,
                    a_az, a_el, G0, f0, df, Nf, Pt, c, out):
    M = tx.shape[0]
    P = fac.shape[0]
    ln2 = np.log(2.0)
    for m in nb.prange(M):
        for p in range(P):
            xp, yp = fac[p, 0], fac[p, 1]
            dtx, dty, dtz = tx[m, 0] - xp, tx[m, 1] - yp, tx[m, 2]
            RT = np.sqrt(dtx * dtx + dty * dty + dtz * dtz)
            drx, dry, drz = rx[m, 0] - xp, rx[m, 1] - yp, rx[m, 2]
            RR = np.sqrt(drx * drx + dry * dry + drz * drz)
            ix, iy, iz = dtx / RT, dty / RT, dtz / RT
            sx, sy, sz = drx / RR, dry / RR, drz / RR
            ci, cs = iz, sz
            si2, ss2 = 1.0 - ci * ci, 1.0 - cs * cs
            si, ss = np.sqrt(max(si2, 0.0)), np.sqrt(max(ss2, 0.0))
            # cos(phi_s): s_h frente a la direccion directa -i_h
            nh = np.sqrt(ix * ix + iy * iy) * np.sqrt(sx * sx + sy * sy)
            cphi = -(ix * sx + iy * sy) / nh if nh > 0 else 1.0
            qi = _kz(eps, si2)
            qs = _kz(eps, ss2)
            a_hh = (eps - 1) * cphi / ((ci + qi) * (cs + qs))
            a_vv = (eps - 1) * (eps * si * ss - cphi * qi * qs) / ((eps * ci + qi) * (eps * cs + qs))
            aeff = wTM * a_vv + wTE * a_hh
            K2 = k0 * k0 * (ss2 + si2 - 2 * ss * si * cphi)
            W = cl * cl / (4 * np.pi) * np.exp(-K2 * cl * cl / 4)
            sig0 = 8 * k0 ** 4 * sh * sh * ci * ci * cs * cs * (aeff.real ** 2 + aeff.imag ** 2) * W
            # ganancias (direccion antena -> faceta = -i, -s)
            cb = -(ix * bt[m, 0] + iy * bt[m, 1] + iz * bt[m, 2])
            daz = np.arctan2(-(ix * ht[m, 0] + iy * ht[m, 1] + iz * ht[m, 2]), cb)
            dele = np.arctan2(-(ix * et[m, 0] + iy * et[m, 1] + iz * et[m, 2]), cb)
            gT = G0 * np.exp(-4 * ln2 * ((daz / a_az) ** 2 + (dele / a_el) ** 2))
            cb = -(sx * br[m, 0] + sy * br[m, 1] + sz * br[m, 2])
            daz = np.arctan2(-(sx * hr[m, 0] + sy * hr[m, 1] + sz * hr[m, 2]), cb)
            dele = np.arctan2(-(sx * er[m, 0] + sy * er[m, 1] + sz * er[m, 2]), cb)
            gR = G0 * np.exp(-4 * ln2 * ((daz / a_az) ** 2 + (dele / a_el) ** 2))
            A = np.sqrt(sig0 * gT * gR) / (RT * RR)
            tau = (RT + RR) / c
            ph0 = -2 * np.pi * f0 * tau
            z = A * amp[p] * (np.cos(ph0) + 1j * np.sin(ph0))
            st = np.cos(-2 * np.pi * df * tau) + 1j * np.sin(-2 * np.pi * df * tau)
            for k in range(Nf):
                lam = c / (f0 + k * df)
                out[m, k] += z * (np.sqrt(Pt) * lam / (4 * np.pi) ** 1.5)
                z *= st


def surface_facets(half, step, taper_frac=0.15, center=(0.0, 0.0)):
    """Facetas cuadradas en z=0 (lado 'step') sobre [-half, half]^2 con ventana de Tukey."""
    x = np.arange(-half, half + step / 2, step)
    X, Y = np.meshgrid(x + center[0], x + center[1], indexing="ij")
    fac = np.stack([X.ravel(), Y.ravel()], 1)

    def tukey(u):
        a = taper_frac
        w = np.ones_like(u)
        if a <= 0:
            return w
        r = np.abs(u) / half
        m = r > 1 - a
        w[m] = 0.5 * (1 + np.cos(np.pi * (r[m] - (1 - a)) / a))
        return w

    w = tukey(fac[:, 0] - center[0]) * tukey(fac[:, 1] - center[1])
    return fac, w, step * step


def diffuse_echo(sysm, Ptx, Prx, rough, rng, facets=None):
    """Eco difuso de la superficie rugosa (SPM) para una realizacion de la superficie."""
    if facets is None:
        facets = surface_facets(rough["patch_half"], rough["facet"], rough.get("taper", 0.15))
    fac, w, dA = facets
    a = (rng.standard_normal(len(fac)) + 1j * rng.standard_normal(len(fac))) / np.sqrt(2)
    amp = a * np.sqrt(dA) * w
    bt, ht, et = antenna_frames(Ptx, sysm.aim)
    br, hr, er = antenna_frames(Prx, sysm.aim)
    wTM, wTE = sysm.pol_weights()
    out = np.zeros((Ptx.shape[0], sysm.Nf), np.complex128)
    df = sysm.freqs[1] - sysm.freqs[0]
    _diffuse_kernel(np.ascontiguousarray(Ptx), np.ascontiguousarray(Prx), bt, ht, et, br, hr, er,
                    np.ascontiguousarray(fac), amp.astype(np.complex128), complex(sysm.eps),
                    2 * np.pi * sysm.fc / sysm.c, float(rough["sigma_h"]), float(rough["corr_len"]),
                    float(wTM), float(wTE), np.radians(sysm.bw_az), np.radians(sysm.bw_el),
                    10 ** (sysm.G0 / 10), float(sysm.freqs[0]), float(df), sysm.Nf, float(sysm.Pt),
                    float(sysm.c), out)
    return out
