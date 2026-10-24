"""Exp. C - Trayectorias completas en 3D: helice base, F1 (helice girada) y F2 (cono de
incidencia constante), con los conos de incidencia constante con vertice en el centro de la
superficie.  Usa geometry.familia_C (la misma funcion que exp_C_brewster.py).
Uso: python exp_C_fig_trayectorias3d.py
"""
import numpy as np
import matplotlib.pyplot as plt

from bisar import config, em, geometry, plotting as P

NAME = "expC"
p = config.load(NAME)
out = config.outdir(NAME)
t, P0 = geometry.tx_trajectory(p)
thB = em.pseudo_brewster_deg(em.eps_complex(p["suelos"]["S1"]["eps_r"]))[0]
th_base = geometry.incidence_deg(P0).mean()
CB, CG = P.C["s8"], P.C["ink2"]
LMAX = 200.0                              # longitud de las generatrices de los conos dibujados


def cone(ax, theta_deg, color, alpha=0.10, L=LMAX):
    """Superficie de incidencia constante theta: vertice en (0,0,0), eje vertical."""
    ph = np.linspace(0, 2 * np.pi, 73)
    ll = np.linspace(0, L, 12)
    PH, LL = np.meshgrid(ph, ll)
    s, c = np.sin(np.radians(theta_deg)), np.cos(np.radians(theta_deg))
    ax.plot_surface(LL * s * np.cos(PH), LL * s * np.sin(PH), LL * c, color=color, alpha=alpha,
                    linewidth=0, shade=False)
    ax.plot(L * s * np.cos(ph), L * s * np.sin(ph), L * c + 0 * ph, color=color, lw=0.6, alpha=0.6)


def rays(ax, Ptraj, idx, color):
    """Lineas de vista desde el centro de la superficie a algunos pulsos, con su incidencia."""
    for k in idx:
        q = Ptraj[k]
        ax.plot([0, q[0]], [0, q[1]], [0, q[2]], color=color, lw=0.8, ls="--")
        th = geometry.incidence_deg(Ptraj[k:k + 1])[0]
        ax.text(q[0] * 1.06, q[1] * 1.06, q[2] + 6, f"{th:.1f}°", fontsize=8, color=P.C["ink"])


def ground(ax):
    g = np.linspace(-190, 190, 2)
    G1, G2 = np.meshgrid(g, g)
    ax.plot_surface(G1, G2, 0 * G1, color="#c9bfae", alpha=0.25, linewidth=0, shade=False)
    ax.scatter([0], [0], [0], color=P.C["ink"], s=14)


def path(ax, Ptraj, color, label):
    ax.plot(*Ptraj.T, color=color, lw=2.0, label=label)
    ax.scatter(*Ptraj[0], color=color, s=30, marker="o", edgecolor="k", linewidth=0.5, zorder=5)
    ax.scatter(*Ptraj[-1], color=color, s=34, marker="s", edgecolor="k", linewidth=0.5, zorder=5)


def style(ax, title):
    ax.set_xlim(-190, 190); ax.set_ylim(-190, 190); ax.set_zlim(0, 140)
    ax.set_box_aspect((1, 1, 0.62), zoom=1.08)
    ax.set_xlabel("x [m]", labelpad=-2); ax.set_ylabel("y [m]", labelpad=-2); ax.set_zlabel("z [m]", labelpad=-4)
    ax.tick_params(labelsize=7, pad=-2)
    ax.view_init(elev=14, azim=-58)
    ax.set_title(title, fontsize=9, pad=0)


F1 = geometry.familia_C(P0, "F1", thB)
F2 = geometry.familia_C(P0, "F2", thB)
th0 = geometry.incidence_deg(P0)
th1 = geometry.incidence_deg(F1)
idx = [0, 800, len(t) - 1]

fig = plt.figure(figsize=(16, 6.2))
# (a) helice base
ax = fig.add_subplot(1, 3, 1, projection="3d")
ground(ax)
cone(ax, th0.min(), P.C["mono"], 0.07)
cone(ax, th0.max(), P.C["mono"], 0.07)
path(ax, P0, P.C["mono"], f"hélice base (= F1 en $\\theta_c$={th_base:.1f}°)")
rays(ax, P0, idx, CG)
style(ax, f"(a) Hélice base: 2 vueltas, baja de {P0[0, 2]:.0f} a {P0[-1, 2]:.0f} m\n"
          f"y se abre de {np.hypot(*P0[0, :2]):.0f} a {np.hypot(*P0[-1, :2]):.0f} m: "
          f"la incidencia barre {th0.min():.1f}°→{th0.max():.1f}°\n(conos de {th0.min():.1f}° y {th0.max():.1f}°)")
# (b) F1 en theta_pB
ax = fig.add_subplot(1, 3, 2, projection="3d")
ground(ax)
cone(ax, thB, CB, 0.10)
path(ax, F1, P.C["s7"], "F1 en $\\theta_c=\\theta_{pB}$")
rays(ax, F1, idx, CG)
style(ax, f"(b) F1 en $\\theta_c=\\theta_{{pB}}$: la misma hélice girada hacia abajo;\n"
          f"baja de {F1[0, 2]:.0f} a {F1[-1, 2]:.0f} m y barre {th1.min():.1f}°→{th1.max():.1f}°:\n"
          f"cruza el cono de Brewster (rojo) una sola vez")
# (c) F2 en theta_pB
ax = fig.add_subplot(1, 3, 3, projection="3d")
ground(ax)
cone(ax, thB, CB, 0.14)
path(ax, F2, P.C["bi"], "F2 en $\\theta_c=\\theta_{pB}$")
rays(ax, F2, [0, 267, 533], CG)
style(ax, f"(c) F2 en $\\theta_c=\\theta_{{pB}}$: todos los pulsos sobre el cono de {thB:.2f}°\n"
          f"(vértice en el centro de la superficie) → círculo a z≈{F2[:, 2].mean():.0f} m,\n"
          f"radio ≈{np.hypot(F2[:, 0], F2[:, 1]).mean():.0f} m, recorrido 2 veces (Δz < {np.ptp(F2[:, 2]):.1f} m)")
fig.text(0.5, 0.015, "Cono de incidencia constante = lugar geométrico de los puntos vistos desde el centro de la "
         "superficie (0,0,0) con el mismo ángulo θ respecto de la vertical. Círculo = inicio (pulso 0); "
         "cuadrado = fin; líneas discontinuas = línea de vista con su incidencia.", ha="center", fontsize=8.5, color=CG)
fig.subplots_adjust(left=0.0, right=0.97, top=0.84, bottom=0.06, wspace=0.06)
P.save(fig, "C_trayectorias_3d", out)
print("ok", thB, F2[:, 2].min(), F2[:, 2].max())
