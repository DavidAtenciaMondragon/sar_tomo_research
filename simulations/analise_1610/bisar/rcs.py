"""Mallas cerradas y Optica Fisica (PO) biestatica co-polar para PEC.

Puerto de verifica_ganancia_brewster.m (rcsPO, mallaObjeto, orientarFuera), con:
  * polarizaciones TM/TE generales fuera del plano xz (vectores v/h continuos en s_hat),
  * numero de onda por muestra (barrido en frecuencia),
  * nucleo numba paralelo.
Amplitud compleja de dispersion (fase referida al centro del objeto):
  S = (2 sqrt(pi)/lambda) sum_{n.i>0} [e_r.(n x (k_i x e_i))] I_T(k (i+s)),   sigma = |S|^2
"""
import numpy as np
import numba as nb


# ----------------------------------------------------------------------------- mallas
def _orient_outward(V, F):
    c = V.mean(0)
    p0, p1, p2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    n = np.cross(p1 - p0, p2 - p0)
    flip = np.einsum("ij,ij->i", n, (p0 + p1 + p2) / 3 - c) < 0
    F = F.copy()
    F[flip] = F[flip][:, [0, 2, 1]]
    return F


def _nseg(r, lam):
    return max(12, int(np.ceil(np.pi / np.arccos(max(1 - lam / (100 * r), -1)))))


def mesh_box(lx, ly, lz):
    V = 0.5 * np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                        [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float) * [lx, ly, lz]
    F = np.array([[1, 3, 2], [1, 4, 3], [5, 6, 7], [5, 7, 8], [1, 2, 6], [1, 6, 5],
                  [4, 8, 7], [4, 7, 3], [1, 5, 8], [1, 8, 4], [2, 3, 7], [2, 7, 6]]) - 1
    return V, _orient_outward(V, F)


def mesh_cylinder(r, L, lam):
    N = _nseg(r, lam)
    ph = 2 * np.pi * np.arange(N) / N
    ring = np.stack([r * np.cos(ph), r * np.sin(ph)], 1)
    V = np.vstack([np.c_[ring, -L / 2 * np.ones(N)], np.c_[ring, L / 2 * np.ones(N)],
                   [[0, 0, -L / 2], [0, 0, L / 2]]])
    i1 = np.arange(N)
    i2 = (i1 + 1) % N
    F = np.vstack([np.c_[i1, i2, N + i2], np.c_[i1, N + i2, N + i1],
                   np.c_[np.full(N, 2 * N), i2, i1], np.c_[np.full(N, 2 * N + 1), N + i1, N + i2]])
    return V, _orient_outward(V, F)


def mesh_cone(r, L, lam):
    N = _nseg(r, lam)
    ph = 2 * np.pi * np.arange(N) / N
    V = np.vstack([np.c_[r * np.cos(ph), r * np.sin(ph), -L / 2 * np.ones(N)], [[0, 0, L / 2], [0, 0, -L / 2]]])
    i1 = np.arange(N)
    i2 = (i1 + 1) % N
    F = np.vstack([np.c_[i1, i2, np.full(N, N)], np.c_[np.full(N, N + 1), i2, i1]])
    return V, _orient_outward(V, F)


def mesh_sphere(r, lam):
    N = _nseg(r, lam)
    M = int(np.ceil(N / 2))
    ph = 2 * np.pi * np.arange(N) / N
    th = np.pi * np.arange(1, M) / M
    X = r * np.outer(np.sin(th), np.cos(ph))
    Y = r * np.outer(np.sin(th), np.sin(ph))
    Z = r * np.outer(np.cos(th), np.ones(N))
    # vertice 0 = polo norte; luego meridiano b, anillo a -> 1 + b*(M-1) + a
    V = np.vstack([[0, 0, r], np.c_[X.T.ravel(), Y.T.ravel(), Z.T.ravel()], [0, 0, -r]])
    idx = lambda a, b: 1 + (b % N) * (M - 1) + a
    F = []
    for b in range(N):
        F.append([0, idx(0, b), idx(0, b + 1)])
        for a in range(M - 2):
            F.append([idx(a, b), idx(a + 1, b), idx(a + 1, b + 1)])
            F.append([idx(a, b), idx(a + 1, b + 1), idx(a, b + 1)])
        F.append([len(V) - 1, idx(M - 2, b + 1), idx(M - 2, b)])
    return V, _orient_outward(V, np.array(F))


def make_mesh(kind, dim, lam):
    if kind == "cubo":
        return mesh_box(dim, dim, dim)
    if kind == "caja":
        return mesh_box(*dim)
    if kind in ("cilindro", "tubo"):
        return mesh_cylinder(dim[0], dim[1], lam)
    if kind == "cono":
        return mesh_cone(dim[0], dim[1], lam)
    if kind == "esfera":
        return mesh_sphere(dim, lam)
    raise ValueError(kind)


def rot_xyz(rx_deg=0.0, ry_deg=0.0, rz_deg=0.0):
    """R = Rz Ry Rx (grados)."""
    a, b, c = np.radians([rx_deg, ry_deg, rz_deg])
    Rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    Ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    Rz = np.array([[np.cos(c), -np.sin(c), 0], [np.sin(c), np.cos(c), 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def random_so3(n, rng):
    """Rotaciones uniformes en SO(3) (cuaternion gaussiano normalizado), igual que
    orientacionesAleatorias(...,'so3') de verifica_ganancia_brewster.m."""
    q = rng.standard_normal((n, 4))
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    w, x, y, z = q.T
    R = np.empty((n, 3, 3))
    R[:, 0] = np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], 1)
    R[:, 1] = np.stack([2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], 1)
    R[:, 2] = np.stack([2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], 1)
    return R


# ------------------------------------------------------------------- polarizacion
def pol_vectors(kdir, pol):
    """Vector de polarizacion para una onda que se propaga en kdir (N,3).
    h = z x k / |z x k| (TE, horizontal);  v = h x k (TM, en el plano de incidencia).
    Ambos son continuos en kdir (salvo nadir): mono y bi convergen cuando dphi -> 0."""
    z = np.array([0.0, 0.0, 1.0])
    h = np.cross(z, kdir)
    nh = np.linalg.norm(h, axis=1, keepdims=True)
    bad = nh[:, 0] < 1e-12                      # nadir exacto: base arbitraria (y)
    h[bad] = [0.0, 1.0, 0.0]
    nh[bad] = 1.0
    h /= nh
    if pol.upper() == "TE":
        return h
    return np.cross(h, kdir)


def aligned_rx(ihat, shat):
    """Rx 'alineado' (polarizacion adaptada al Tx): Tx TM; el Rx recibe la componente de e_i
    transversal a s_hat, e_r = (e_i - (e_i.s) s)/|.|.  En mono (s = i) coincide con TM-TM.
    Devuelve (e_i, e_r, cos_psi) con cos_psi = e_r . v(s) (fraccion TM del Rx; el resto es TE),
    necesario para la transmitancia de subida.  Requiere un Rx de doble polarizacion."""
    ei = pol_vectors(-ihat, "TM")
    er = ei - np.einsum("ij,ij->i", ei, shat)[:, None] * shat
    er /= np.linalg.norm(er, axis=1, keepdims=True)
    cpsi = np.abs(np.einsum("ij,ij->i", er, pol_vectors(shat, "TM")))
    return ei, er, cpsi


# ------------------------------------------------------------------------- nucleo PO
@nb.njit(cache=True, inline="always")
def _phi1_s(t):
    if abs(t) < 1e-4:
        return 1 + 0.5j * t - t * t / 6
    return (np.exp(1j * t) - 1) / (1j * t)


@nb.njit(cache=True, inline="always")
def _g_s(x, y):
    if abs(x) > abs(y):
        x, y = y, x
    if abs(y) < 1e-3:
        return (0.5 + 1j * (x + y) / 6 - (x * x + x * y + y * y) / 24
                - 1j * (x ** 3 + x * x * y + x * y * y + y ** 3) / 120)
    return (np.exp(1j * y) * _phi1_s(x - y) - _phi1_s(x)) / (1j * y)


@nb.njit(parallel=True, cache=True)
def _po_kernel(I, S, EI, ER, kk, p0, p1, p2, n, area2):
    nd, m = I.shape[0], p0.shape[0]
    out = np.empty(nd, dtype=np.complex128)
    for d in nb.prange(nd):
        ix, iy, iz = I[d, 0], I[d, 1], I[d, 2]
        # k_i = -i ; h_i = k_i x e_i ; q = h_i x e_r ; factor = n . q = e_r . (n x h_i)
        kix, kiy, kiz = -ix, -iy, -iz
        ex, ey, ez = EI[d, 0], EI[d, 1], EI[d, 2]
        hx = kiy * ez - kiz * ey
        hy = kiz * ex - kix * ez
        hz = kix * ey - kiy * ex
        rx, ry, rz = ER[d, 0], ER[d, 1], ER[d, 2]
        qx = hy * rz - hz * ry
        qy = hz * rx - hx * rz
        qz = hx * ry - hy * rx
        k = kk[d]
        wx, wy, wz = k * (ix + S[d, 0]), k * (iy + S[d, 1]), k * (iz + S[d, 2])
        acc = 0j
        for j in range(m):
            if ix * n[j, 0] + iy * n[j, 1] + iz * n[j, 2] <= 0.0:
                continue
            fac = n[j, 0] * qx + n[j, 1] * qy + n[j, 2] * qz
            a = wx * p0[j, 0] + wy * p0[j, 1] + wz * p0[j, 2]
            b = wx * p1[j, 0] + wy * p1[j, 1] + wz * p1[j, 2]
            c = wx * p2[j, 0] + wy * p2[j, 1] + wz * p2[j, 2]
            acc += fac * area2[j] * np.exp(1j * a) * _g_s(b - a, c - a)
        # 2 sqrt(pi) / lambda = 2 sqrt(pi) k / (2 pi)
        out[d] = acc * (np.sqrt(np.pi) * k / np.pi)
    return out


class PO:
    """Objeto PEC facetado (opcionalmente rotado por R). amplitude(i, s, k, pol) -> S complejo."""

    def __init__(self, V, F, R=None):
        V = np.asarray(V, float)
        if R is not None:
            V = V @ np.asarray(R).T
        F = np.asarray(F, int)
        self.V, self.F = V, F
        self.p0 = np.ascontiguousarray(V[F[:, 0]])
        self.p1 = np.ascontiguousarray(V[F[:, 1]])
        self.p2 = np.ascontiguousarray(V[F[:, 2]])
        cr = np.cross(self.p1 - self.p0, self.p2 - self.p0)
        self.area2 = np.linalg.norm(cr, axis=1)
        self.n = np.ascontiguousarray(cr / self.area2[:, None])

    def amplitude(self, ihat, shat, k, pol="TM"):
        ihat = np.ascontiguousarray(np.atleast_2d(ihat), float)
        shat = np.ascontiguousarray(np.atleast_2d(shat), float)
        k = np.ascontiguousarray(np.broadcast_to(np.asarray(k, float), (ihat.shape[0],)))
        if pol.upper() == "TMAL":
            ei, er, _ = aligned_rx(ihat, shat)
        else:
            ei = pol_vectors(-ihat, pol)
            er = pol_vectors(shat, pol)
        ei, er = np.ascontiguousarray(ei), np.ascontiguousarray(er)
        return _po_kernel(ihat, shat, ei, er, k, self.p0, self.p1, self.p2, self.n, self.area2)

    def sigma(self, ihat, shat, k, pol="TM"):
        return np.abs(self.amplitude(ihat, shat, k, pol)) ** 2
