"""Mallas triangulares cerradas (normales hacia fuera) generadas con gmsh (OpenCASCADE).

Todas las funciones devuelven (V, F):
    V : (N, 3) float  -> vértices [m]
    F : (M, 3) int    -> triángulos (índices en V), orientados con normal saliente
"""

import numpy as np
import gmsh


def _mesh_current_model(h, curv=0):
    if curv:  # nº de elementos por cada 2*pi de curvatura
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", curv)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0 if curv else h)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay
    gmsh.model.mesh.generate(2)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    V = coords.reshape(-1, 3)
    index = {t: i for i, t in enumerate(tags)}
    F = []
    for dim, tag in gmsh.model.getEntities(2):
        etypes, _, enodes = gmsh.model.mesh.getElements(dim, tag)
        for et, en in zip(etypes, enodes):
            if et == 2:  # triángulo de 3 nodos
                F.append(np.vectorize(index.get)(en).reshape(-1, 3))
    F = np.vstack(F)
    # compactar vértices no usados
    used, F = np.unique(F, return_inverse=True)
    return V[used], F.reshape(-1, 3)


def _orient_outward(V, F):
    """Asegura normales salientes (válido para cuerpos estrellados respecto al centroide)."""
    c = V.mean(axis=0)
    p0, p1, p2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    n = np.cross(p1 - p0, p2 - p0)
    flip = np.einsum("ij,ij->i", n, (p0 + p1 + p2) / 3 - c) < 0
    F = F.copy()
    F[flip] = F[flip][:, [0, 2, 1]]
    return F


def _build(add_shape, h, curv=0):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("m")
    add_shape()
    gmsh.model.occ.synchronize()
    try:
        V, F = _mesh_current_model(h, curv)
    finally:
        gmsh.finalize()
    return V, _orient_outward(V, F)


def box(lx, ly, lz, h, curv=0):
    """Paralelepípedo centrado en el origen."""
    return _build(lambda: gmsh.model.occ.addBox(-lx / 2, -ly / 2, -lz / 2, lx, ly, lz), h, curv)


def cube(a, h, curv=0):
    return box(a, a, a, h, curv)


def cylinder(radius, length, h, curv=0):
    """Cilindro con tapas, eje z, centrado en el origen."""
    return _build(lambda: gmsh.model.occ.addCylinder(0, 0, -length / 2, 0, 0, length, radius), h, curv)


def sphere(radius, h, curv=0):
    return _build(lambda: gmsh.model.occ.addSphere(0, 0, 0, radius), h, curv)


def cone(radius, length, h, curv=0):
    """Cono con base en z=-L/2 y vértice en z=+L/2."""
    return _build(lambda: gmsh.model.occ.addCone(0, 0, -length / 2, 0, 0, length, radius, 0), h, curv)


def plate(lx, ly, h):
    """Placa fina (superficie abierta, normal +z). Útil para validar PO."""
    def add():
        gmsh.model.occ.addRectangle(-lx / 2, -ly / 2, 0, lx, ly)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("m")
    add()
    gmsh.model.occ.synchronize()
    try:
        V, F = _mesh_current_model(h)
    finally:
        gmsh.finalize()
    p0, p1, p2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    n = np.cross(p1 - p0, p2 - p0)
    F[n[:, 2] < 0] = F[n[:, 2] < 0][:, [0, 2, 1]]
    return V, F


def from_file(path):
    """Carga cualquier malla (STL, OBJ, PLY...) con trimesh."""
    import trimesh
    m = trimesh.load(path, force="mesh")
    m.fix_normals()
    return np.asarray(m.vertices, float), np.asarray(m.faces, int)
