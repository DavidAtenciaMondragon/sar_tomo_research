"""SCRIPT PRINCIPAL: simula la RCS monoestática de un objeto PEC y la muestra en 3D.

Genera una página HTML interactiva (se abre sola en el navegador) con:
  1. Escenario: objeto, radar, dirección de incidencia k y campo E (TM o TE).
  2. Malla 3D que se ha creado para el cálculo.
  3. Diagrama 3D de RCS (radio y color = dBsm).
  4. Mapa azimut x elevación.
  5-6. Cortes en azimut y en elevación.

Ejemplos:
    .venv\\Scripts\\python simular.py cubo --dim 1.0 --freq 400e6
    .venv\\Scripts\\python simular.py cilindro --dim 0.3 1.5 --freq 400e6 --vista 45 10
    .venv\\Scripts\\python simular.py esfera --dim 0.5 --freq 400e6 --pol TE
    .venv\\Scripts\\python simular.py stl --archivo modelo.stl --freq 10e9

Convenciones
  * Azimut en el plano xy desde +x; elevación desde el plano xy (+90° = +z).
  * TM: campo E en el plano que contiene el eje z y la dirección de
    propagación (vector unitario de elevación).  TE: campo E horizontal.
  * TM-TM = transmitir TM y recibir TM (co-polar, "VV").
"""

import argparse
import time
import webbrowser
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from rcs import geometry as g
from rcs.po import POSolver, C0, dbsm, directions
from rcs.mom import pol_vectors

MAX_DOF_MOM = 12000  # ~2.3 GB de matriz; por encima se recomienda PO


# ----------------------------------------------------------------------------- malla
def crear_malla(forma, dim, lam, metodo, archivo=None):
    if forma == "stl":
        return g.from_file(archivo)
    if metodo == "mom":
        h, curv = lam / 10, 0  # MoM: ~10 elementos por longitud de onda
    else:
        # PO integra la fase exactamente por triángulo: sólo hay que seguir la
        # curvatura (flecha < lambda/100, error < ~0.2 dB).
        h = max(dim) / 2
        curv = int(np.ceil(2 * np.pi / np.sqrt(0.08 * lam / dim[0])))
    if forma == "cubo":
        return g.cube(dim[0], h)
    if forma == "caja":
        return g.box(*dim, h)
    if forma == "cilindro":
        return g.cylinder(dim[0], dim[1], h, curv)
    if forma == "esfera":
        return g.sphere(dim[0], h, curv)
    if forma == "cono":
        return g.cone(dim[0], dim[1], h, curv)
    raise ValueError(forma)


def tamano_electrico(forma, dim, lam):
    if forma in ("cubo",):
        L = dim[0] * np.sqrt(3)
    elif forma == "caja":
        L = np.linalg.norm(dim)
    elif forma in ("cilindro", "cono"):
        L = np.hypot(2 * dim[0], dim[1])
    elif forma == "esfera":
        L = 2 * dim[0]
    else:
        L = np.nan
    return L / lam


# ----------------------------------------------------------------------------- figuras
def _aristas(V, F):
    E = np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    E = np.unique(E, axis=0)
    P = np.full((len(E) * 3, 3), np.nan)
    P[0::3], P[1::3] = V[E[:, 0]], V[E[:, 1]]
    return P


def _flecha(fig, p0, vec, color, nombre, row, col, texto=None, ancho=6):
    p1 = p0 + vec
    fig.add_trace(go.Scatter3d(x=[p0[0], p1[0]], y=[p0[1], p1[1]], z=[p0[2], p1[2]],
                               mode="lines", line=dict(color=color, width=ancho),
                               name=nombre, showlegend=False, hoverinfo="name"), row, col)
    n = np.linalg.norm(vec)
    fig.add_trace(go.Cone(x=[p1[0]], y=[p1[1]], z=[p1[2]], u=[vec[0] / n], v=[vec[1] / n],
                          w=[vec[2] / n], anchor="tip", sizemode="absolute", sizeref=0.18 * n,
                          colorscale=[[0, color], [1, color]], showscale=False,
                          name=nombre, hoverinfo="name"), row, col)
    if texto:
        q = p1 + 0.12 * vec
        fig.add_trace(go.Scatter3d(x=[q[0]], y=[q[1]], z=[q[2]], mode="text", text=[texto],
                                   textfont=dict(color=color, size=13), showlegend=False,
                                   hoverinfo="skip"), row, col)


def _objeto(fig, V, F, row, col, color="#9aa4b1", opacidad=1.0, bordes=False):
    fig.add_trace(go.Mesh3d(x=V[:, 0], y=V[:, 1], z=V[:, 2], i=F[:, 0], j=F[:, 1], k=F[:, 2],
                            color=color, opacity=opacidad, flatshading=True,
                            lighting=dict(ambient=0.45, diffuse=0.8, specular=0.2),
                            name="objeto", hoverinfo="skip", showlegend=False), row, col)
    if bordes:
        P = _aristas(V, F)
        fig.add_trace(go.Scatter3d(x=P[:, 0], y=P[:, 1], z=P[:, 2], mode="lines",
                                   line=dict(color="#1f2937", width=1.5), name="aristas",
                                   hoverinfo="skip", showlegend=False), row, col)


def figura(V, F, az, el, s_db, cross_db, info, vista, pol):
    L = np.ptp(V, axis=0).max()
    fig = make_subplots(
        rows=3, cols=2, vertical_spacing=0.06, horizontal_spacing=0.06,
        row_heights=[0.38, 0.38, 0.24],
        specs=[[{"type": "scene"}, {"type": "scene"}],
               [{"type": "scene"}, {"type": "xy"}],
               [{"type": "xy"}, {"type": "xy"}]],
        subplot_titles=[
            f"1. Escenario — radar en az {vista[0]:g}°, el {vista[1]:g}°, polarización {pol}",
            f"2. Malla: {len(F)} triángulos (arista media {info['arista_lambda']})",
            f"3. RCS {pol}-{pol} en 3D (radio y color en dBsm)",
            "4. Mapa RCS azimut × elevación [dBsm]",
            "5. Corte en azimut (el = 0°)",
            f"6. Corte en elevación (az = {vista[0]:g}°)"])

    # --- 1. Escenario
    _objeto(fig, V, F, 1, 1)
    for eje, vec, c in (("x", [1, 0, 0], "#d62728"), ("y", [0, 1, 0], "#2ca02c"),
                        ("z", [0, 0, 1], "#1f77b4")):
        _flecha(fig, np.zeros(3), 1.2 * L * np.array(vec, float), c, eje, 1, 1, eje, ancho=3)
    d = directions([vista[0]], [vista[1]])[0, 0]
    v, h = (p[0, 0] for p in pol_vectors([vista[0]], [vista[1]]))
    e = v if pol == "TM" else h
    radar = 2.2 * L * d
    fig.add_trace(go.Scatter3d(x=[radar[0]], y=[radar[1]], z=[radar[2]], mode="markers+text",
                               marker=dict(size=7, color="#111827", symbol="diamond"),
                               text=["  Radar"], textposition="middle right",
                               showlegend=False, hoverinfo="skip"), 1, 1)
    _flecha(fig, radar - 0.15 * L * d, -1.2 * L * d, "#7c3aed", "k incidente", 1, 1, "k")
    _flecha(fig, radar - 0.15 * L * d, 0.6 * L * e, "#ea580c", f"E ({pol})", 1, 1, f"E {pol}")
    # frente de onda plano
    c0 = radar - 0.5 * L * d
    sq = np.array([c0 + 0.35 * L * (a * v + b * h) for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1))])
    fig.add_trace(go.Mesh3d(x=sq[:, 0], y=sq[:, 1], z=sq[:, 2], i=[0, 0], j=[1, 2], k=[2, 3],
                            color="#a78bfa", opacity=0.25, hoverinfo="skip",
                            showlegend=False), 1, 1)

    # --- 2. Malla
    _objeto(fig, V, F, 1, 2, color="#cbd5e1", bordes=True)

    # --- 3. Diagrama 3D de RCS
    vmax = np.nanmax(s_db)
    vmin = vmax - 40
    azc = np.append(az, az[0] + 360)  # cerrar la superficie
    sc = np.vstack([s_db, s_db[:1]])
    A, E = np.meshgrid(np.radians(azc), np.radians(el), indexing="ij")
    r = np.clip(sc - vmin, 0, None)
    fig.add_trace(go.Surface(
        x=r * np.cos(E) * np.cos(A), y=r * np.cos(E) * np.sin(A), z=r * np.sin(E),
        surfacecolor=sc, cmin=vmin, cmax=vmax, colorscale="Turbo",
        customdata=np.dstack([np.degrees(A), np.degrees(E), sc]),
        hovertemplate="az %{customdata[0]:.0f}°<br>el %{customdata[1]:.0f}°"
                      "<br>%{customdata[2]:.1f} dBsm<extra></extra>",
        colorbar=dict(title="dBsm", x=1.01, y=0.5, len=0.33)), 2, 1)

    # --- 4. Mapa az x el
    fig.add_trace(go.Heatmap(x=az, y=el, z=s_db.T, zmin=vmin, zmax=vmax, colorscale="Turbo",
                             showscale=False,
                             hovertemplate="az %{x}°<br>el %{y}°<br>%{z:.1f} dBsm<extra></extra>"),
                  2, 2)
    fig.add_trace(go.Scatter(x=[vista[0]], y=[vista[1]], mode="markers",
                             marker=dict(symbol="x", size=12, color="white",
                                         line=dict(width=2, color="black")),
                             name="radar", showlegend=False), 2, 2)
    fig.update_xaxes(title_text="azimut [°]", row=2, col=2)
    fig.update_yaxes(title_text="elevación [°]", row=2, col=2)

    # --- 5/6. Cortes
    i_el0 = np.argmin(np.abs(el))
    i_az = np.argmin(np.abs(((az - vista[0] + 180) % 360) - 180))
    for col, x, y, yc, xt in ((1, az, s_db[:, i_el0], None if cross_db is None else cross_db[:, i_el0], "azimut [°]"),
                              (2, el, s_db[i_az, :], None if cross_db is None else cross_db[i_az, :], "elevación [°]")):
        fig.add_trace(go.Scatter(x=x, y=y, mode="lines", line=dict(color="#2563eb", width=2),
                                 name=f"{pol}-{pol} (co-polar)", showlegend=(col == 1)), 3, col)
        if yc is not None:
            otra = "TE" if pol == "TM" else "TM"
            fig.add_trace(go.Scatter(x=x, y=yc, mode="lines", line=dict(color="#9ca3af", width=1),
                                     name=f"{otra}-{pol} (cruzada)", showlegend=(col == 1)), 3, col)
        fig.update_xaxes(title_text=xt, row=3, col=col)
        fig.update_yaxes(title_text="RCS [dBsm]", range=[vmin, vmax + 3], row=3, col=col)

    escena = dict(aspectmode="data", xaxis_title="x [m]", yaxis_title="y [m]",
                  zaxis_title="z [m]", camera=dict(eye=dict(x=1.5, y=1.2, z=0.9)))
    # cámara del escenario: de lado respecto al radar para que no tape el objeto
    ca, ce = np.radians(vista[0] - 75), np.radians(25)
    ojo = dict(zip("xyz", 1.9 * np.array([np.cos(ce) * np.cos(ca), np.cos(ce) * np.sin(ca), np.sin(ce)])))
    fig.update_layout(
        scene=dict(escena, camera=dict(eye=ojo)), scene2=escena,
        scene3=dict(aspectmode="data", xaxis_title="", yaxis_title="", zaxis_title="",
                    camera=dict(eye=dict(x=1.5, y=1.2, z=0.9))),
        height=1650, template="plotly_white", title=dict(text=info["titulo"], x=0.5),
        legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.035),
        margin=dict(l=20, r=20, t=90, b=40))
    return fig


# ----------------------------------------------------------------------------- MATLAB
def guardar_mat(ruta, az, el, s, cross, V, F, a, lam, metodo):
    """Guarda los resultados en .mat (formato MATLAB 5, se lee con load()).

    Las matrices de RCS se guardan como (Nel x Naz): filas = elevación,
    columnas = azimut, para que en MATLAB funcione directamente
    imagesc(az_deg, el_deg, rcs_dbsm).
    """
    from scipy.io import savemat

    i, j = np.unravel_index(np.argmax(s), s.shape)
    otra = "TE" if a.pol == "TM" else "TM"
    cross_t = np.zeros((0, 0)) if cross is None else cross.T
    savemat(ruta, {
        "az_deg": az[None, :],
        "el_deg": el[None, :],
        "rcs_m2": s.T,
        "rcs_dbsm": dbsm(s).T,
        "cross_m2": cross_t,
        "cross_dbsm": dbsm(cross_t) if cross is not None else cross_t,
        "freq_Hz": float(a.freq),
        "lambda_m": float(lam),
        "pol": f"{a.pol}-{a.pol}",
        "pol_cross": f"{a.pol}-{otra}" if cross is not None else "",
        "metodo": metodo.upper(),
        "forma": a.forma,
        "dim_m": np.array(a.dim, float)[None, :],
        "rcs_max_dbsm": float(dbsm(s[i, j])),
        "az_max_deg": float(az[i]),
        "el_max_deg": float(el[j]),
        "mesh_vertices_m": V,
        "mesh_faces": F.astype(np.float64) + 1,  # índices base 1 para MATLAB
        "descripcion": (
            "RCS monoestatica de un PEC. rcs_m2/rcs_dbsm: matriz Nel x Naz "
            "(filas = el_deg, columnas = az_deg). Azimut en el plano xy desde +x; "
            "elevacion desde el plano xy (+90 = +z). TM: E en el plano que contiene "
            "z y la direccion de propagacion (VV); TE: E horizontal (HH). "
            "cross_*: polarizacion cruzada (vacia con PO). "
            "Malla: patch('Vertices',mesh_vertices_m,'Faces',mesh_faces)."),
    }, do_compression=True, oned_as="row")


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("forma", choices=["cubo", "caja", "cilindro", "esfera", "cono", "stl"])
    ap.add_argument("--dim", type=float, nargs="+", default=[1.0],
                    help="cubo: lado | caja: lx ly lz | cilindro/cono: radio largo | esfera: radio  [m]")
    ap.add_argument("--archivo", help="malla STL/OBJ (forma=stl)")
    ap.add_argument("--freq", type=float, default=400e6, help="frecuencia portadora [Hz]")
    ap.add_argument("--pol", choices=["TM", "TE"], default="TM")
    ap.add_argument("--metodo", choices=["auto", "po", "mom"], default="auto",
                    help="auto: MoM si el objeto mide <= 3 lambda, si no PO")
    ap.add_argument("--paso", type=float, help="paso angular [°] (defecto: 3 con MoM, 1 con PO)")
    ap.add_argument("--vista", type=float, nargs=2, default=[30, 20], metavar=("AZ", "EL"),
                    help="dirección del radar que se dibuja en el escenario [°]")
    ap.add_argument("--salida", help="nombre base de los archivos de salida (en resultados/)")
    ap.add_argument("--no-abrir", action="store_true", help="no abrir el navegador")
    a = ap.parse_args()

    lam = C0 / a.freq
    tam = tamano_electrico(a.forma, a.dim, lam)
    metodo = a.metodo
    if metodo == "auto":
        metodo = "mom" if tam <= 3 else "po"
    paso = a.paso or (3.0 if metodo == "mom" else 1.0)

    print(f"Objeto: {a.forma} {a.dim} m | f = {a.freq / 1e6:g} MHz | λ = {lam:.3f} m | "
          f"tamaño ≈ {tam:.2f} λ | método: {metodo.upper()} | pol: {a.pol}-{a.pol}")
    t0 = time.time()
    V, F = crear_malla(a.forma, a.dim, lam, metodo, a.archivo)
    E = np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    arista = np.linalg.norm(np.diff(V[np.unique(E, axis=0)], axis=1), axis=-1).mean()
    print(f"Malla: {len(V)} nodos, {len(F)} triángulos, arista media {arista * 100:.1f} cm (λ/{lam / arista:.0f})")

    if metodo == "mom" and 1.5 * len(F) > MAX_DOF_MOM:
        print(f"AVISO: ~{int(1.5 * len(F))} incógnitas (~{16 * (1.5 * len(F)) ** 2 / 1e9:.1f} GB). "
              "Puede ser muy lento; considera --metodo po.")

    az = np.arange(0, 360, paso)
    el = np.arange(-90, 90 + 1e-9, paso)
    cross = None
    if metodo == "po":
        s = POSolver(V, F).rcs(az, el, a.freq)
        if tam < 3:
            print("AVISO: el objeto mide menos de 3 λ; PO es sólo orientativa aquí (usa --metodo mom).")
    else:
        from rcs.mom import MoMSolver
        tx = "V" if a.pol == "TM" else "H"
        rx_cross = "H" if tx == "V" else "V"
        r = MoMSolver(V, F, a.freq).rcs(az, el, tx=(tx,))
        s, cross = r[tx + tx], r[rx_cross + tx]
    print(f"Cálculo: {len(az) * len(el)} direcciones en {time.time() - t0:.1f} s")

    s_db = dbsm(s)
    cross_db = None if cross is None else dbsm(cross)
    i, j = np.unravel_index(np.argmax(s), s.shape)
    iv = np.argmin(np.abs(((az - a.vista[0] + 180) % 360) - 180))
    jv = np.argmin(np.abs(el - a.vista[1]))
    donde = f"el {el[j]:g}° (polo)" if abs(el[j]) == 90 else f"az {az[i]:g}°, el {el[j]:g}°"
    print(f"RCS máxima: {s_db[i, j]:.2f} dBsm ({s[i, j]:.3g} m²) en {donde}")
    print(f"RCS en la vista del radar (az {az[iv]:g}°, el {el[jv]:g}°): {s_db[iv, jv]:.2f} dBsm")
    print(f"RCS media (esférica): {dbsm(np.average(s, axis=1, weights=np.cos(np.radians(el))).mean()):.2f} dBsm")

    salida = Path("resultados") / (a.salida or f"{a.forma}_{'x'.join(f'{d:g}' for d in a.dim)}m_"
                                              f"{a.freq / 1e6:g}MHz_{a.pol}{a.pol}_{metodo}")
    salida.parent.mkdir(exist_ok=True)
    np.savez(salida.with_suffix(".npz"), az=az, el=el, rcs_m2=s,
             cross_m2=np.array([]) if cross is None else cross,
             V=V, F=F, freq=a.freq, pol=a.pol, metodo=metodo)
    guardar_mat(salida.with_suffix(".mat"), az, el, s, cross, V, F, a, lam, metodo)

    info = dict(
        arista_lambda=f"≈ λ/{lam / arista:.0f}",
        titulo=(f"<b>{a.forma} {a.dim} m — f = {a.freq / 1e6:g} MHz (λ = {lam:.3f} m, "
                f"≈ {tam:.2f} λ)</b> — {metodo.upper()} — {a.pol}-{a.pol} — "
                f"máx {s_db[i, j]:.1f} dBsm en {donde}"))
    fig = figura(V, F, az, el, s_db, cross_db, info, a.vista, a.pol)
    html = salida.with_suffix(".html")
    fig.write_html(html, include_plotlyjs=True)
    print(f"Guardado: {html}, {salida.with_suffix('.mat')} y {salida.with_suffix('.npz')}")
    if not a.no_abrir:
        webbrowser.open(html.resolve().as_uri())


if __name__ == "__main__":
    main()
