"""RCS monoestática "exacta" (full-wave) de un PEC con el Método de los Momentos (bempp-cl).

EFIE con funciones base RWG y test SNC. La matriz se factoriza una sola vez (LU)
y cada ángulo/polarización es sólo un lado derecho más -> barrer miles de
ángulos cuesta casi lo mismo que uno.

Coste: memoria ~ 16 * Ndof^2 bytes, con Ndof ~ 1.5 * n_triángulos y una malla de
~lambda/10.  Práctico hasta objetos de ~2-4 lambda con matriz densa en un PC.
Para frecuencias de resonancia interna de un cuerpo cerrado la EFIE pierde
precisión (conocido); se puede comprobar moviendo ligeramente la frecuencia.
"""

import numpy as np
import scipy.linalg as sla
import bempp_cl.api as bempp
from bempp_cl.api.operators.boundary.maxwell import electric_field
from bempp_cl.api.operators.far_field.maxwell import electric_field as ff_electric
from bempp_cl.api.integration.triangle_gauss import rule
from scipy.sparse import coo_matrix

from .po import C0, directions


def pol_vectors(az_deg, el_deg):
    """Vectores de polarización V (en elevación) y H (en azimut) para cada dirección."""
    az, el = np.meshgrid(np.radians(az_deg), np.radians(el_deg), indexing="ij")
    v = np.stack([-np.sin(el) * np.cos(az), -np.sin(el) * np.sin(az), np.cos(el)], axis=-1)
    h = np.stack([-np.sin(az), np.cos(az), np.zeros_like(az)], axis=-1)
    return v, h


class MoMSolver:
    def __init__(self, V, F, freq, quad_order=4, verbose=True):
        self.freq = freq
        self.k = 2 * np.pi * freq / C0
        grid = bempp.Grid(np.asarray(V, float).T, np.asarray(F, np.uint32).T)
        self.grid = grid
        self.rwg = bempp.function_space(grid, "RWG", 0)
        self.snc = bempp.function_space(grid, "SNC", 0)
        self.ndof = self.rwg.global_dof_count
        if verbose:
            print(f"  MoM: {grid.number_of_elements} triángulos, {self.ndof} incógnitas")

        A = electric_field(self.rwg, self.rwg, self.snc, self.k,
                           assembler="dense").weak_form().A
        self.lu = sla.lu_factor(A)

        # Cuadratura: puntos globales, pesos*jacobiano, normales y matrices base
        loc, w = rule(quad_order)
        nq = loc.shape[1]
        gd = grid.data()
        verts = gd.vertices[:, gd.elements].transpose(2, 1, 0)  # (Ne,3 vértices,3)
        s, t = loc
        self.xq = (verts[:, 0, None, :] * (1 - s - t)[None, :, None]
                   + verts[:, 1, None, :] * s[None, :, None]
                   + verts[:, 2, None, :] * t[None, :, None]).reshape(-1, 3)
        self.wq = (gd.integration_elements[:, None] * w[None, :]).ravel()
        self.nq = np.repeat(gd.normals, nq, axis=0)
        self.Mrwg = self._basis_map(self.rwg, loc)
        self.Msnc = self._basis_map(self.snc, loc)

    @staticmethod
    def _basis_map(space, loc):
        """Lista de 3 matrices dispersas (Nq_total x Ndof): componente c de la base en cada punto."""
        ne = space.grid.number_of_elements
        nq = loc.shape[1]
        rows, cols, vals = [], [], []
        for e in range(ne):
            val = space.evaluate(e, loc)  # (3, nshape, nq)
            dofs = space.local2global[e]
            for j, d in enumerate(dofs):
                rows.append(e * nq + np.arange(nq))
                cols.append(np.full(nq, d))
                vals.append(val[:, j, :])
        rows = np.concatenate(rows)
        cols = np.concatenate(cols)
        vals = np.concatenate(vals, axis=1)
        T = space.dof_transformation
        shape = (ne * nq, space.grid_dof_count)
        return [(coo_matrix((vals[c], (rows, cols)), shape=shape).tocsr() @ T) for c in range(3)]

    def _apply(self, M, vec):
        """sum_c M_c^T * vec[:, c, ...]."""
        return sum(M[c].T @ vec[:, c] for c in range(3))

    def _eval(self, M, coef):
        return np.stack([M[c] @ coef for c in range(3)], axis=1)  # (Nq,3,...)

    def solve_currents(self, d, p):
        """Coeficientes RWG para ondas planas que llegan desde las direcciones d (N,3) con polarización p (N,3)."""
        phase = np.exp(-1j * self.k * (self.xq @ d.T))  # (Nq,N): onda viajando hacia -d
        einc = p.T[None, :, :] * phase[:, None, :]  # (Nq,3,N)
        tang = np.cross(einc, self.nq[:, :, None], axis=1)  # E x n
        rhs = self._apply(self.Msnc, tang * self.wq[:, None, None])
        return sla.lu_solve(self.lu, -rhs)

    def far_field(self, coef, d):
        """Campo lejano (N,3) en la dirección d[i] producido por la solución coef[:, i]."""
        J = self._eval(self.Mrwg, coef)  # (Nq,3,N)
        phase = np.exp(-1j * self.k * (self.xq @ d.T)) * self.wq[:, None]
        F = np.einsum("qcn,qn->nc", J, phase)
        F -= np.einsum("nc,nc->n", F, d)[:, None] * d  # parte transversal
        return self.ff_scale * F

    ff_scale = 1.0  # se calibra con el operador de campo lejano de bempp (calibrate())

    def calibrate(self):
        """Fija la constante del campo lejano comparando con el operador far-field de bempp."""
        d = np.array([[0.3, -0.5, 0.81]])
        d /= np.linalg.norm(d)
        p = np.cross(d, [[0, 0, 1.0]])
        p /= np.linalg.norm(p)
        coef = self.solve_currents(d, p)
        self.ff_scale = 1.0
        mine = self.far_field(coef, d)[0]
        ref = ff_electric(self.rwg, d.T, self.k).evaluate(
            bempp.GridFunction(self.rwg, coefficients=coef[:, 0]))[:, 0]
        self.ff_scale = np.vdot(mine, ref) / np.vdot(mine, mine)
        err = np.linalg.norm(self.ff_scale * mine - ref) / np.linalg.norm(ref)
        return self.ff_scale, err

    def rcs(self, az_deg, el_deg, tx=("V", "H"), chunk=400):
        """Devuelve dict con sigma [m^2] para 'VV','HH','VH','HV', forma (len(az), len(el)).

        V = TM (campo E en el plano que contiene el eje z y la dirección de
        propagación), H = TE.  `tx` limita las polarizaciones transmitidas
        (cada una es un conjunto de resoluciones).
        """
        if self.ff_scale == 1.0:
            self.calibrate()
        az_deg, el_deg = np.atleast_1d(az_deg), np.atleast_1d(el_deg)
        D = directions(az_deg, el_deg).reshape(-1, 3)
        v, h = (x.reshape(-1, 3) for x in pol_vectors(az_deg, el_deg))
        out = {rx + t: np.empty(len(D)) for t in tx for rx in ("V", "H")}
        for i in range(0, len(D), chunk):
            sl = slice(i, i + chunk)
            for t, ptx in (("V", v[sl]), ("H", h[sl])):
                if t not in tx:
                    continue
                Fs = self.far_field(self.solve_currents(D[sl], ptx), D[sl])
                for rx, prx in (("V", v[sl]), ("H", h[sl])):
                    out[rx + t][sl] = 4 * np.pi * np.abs(np.einsum("nc,nc->n", Fs, prx)) ** 2
        shape = (len(az_deg), len(el_deg))
        return {key: val.reshape(shape) for key, val in out.items()}
