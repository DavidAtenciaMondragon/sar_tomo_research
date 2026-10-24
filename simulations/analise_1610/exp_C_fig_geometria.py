"""Exp. C - Figura de geometria de las familias de trayectoria F1 y F2.

(a) plano vertical (rho, z): helice base, F1 y F2 para varias incidencias centrales theta_c,
    incluidos theta_pB + error; (b) incidencia del Tx en cada pulso, theta_Tx(t);
(c) planta de F2 en theta_pB con el Rx a dphi = 90, 135 y 180 y el punto especular de la
    superficie de cada par; (d) incidencia en el punto especular vs dphi.
Usa geometry.familia_C (la misma funcion que exp_C_brewster.py) y scene.specular_echo.
Uso: python exp_C_fig_geometria.py
"""
import numpy as np
import matplotlib.pyplot as plt

from bisar import config, em, geometry, plotting as P, scene

NAME = "expC"
p = config.load(NAME)
out = config.outdir(NAME)
t, P0 = geometry.tx_trajectory(p)
thB = em.pseudo_brewster_deg(em.eps_complex(p["suelos"]["S1"]["eps_r"]))[0]
th_base = geometry.incidence_deg(P0).mean()
sysm = scene.System(p)
C1, C2, CB, CG = P.C["mono"], P.C["bi"], P.C["s8"], P.C["ink2"]

fig, axs = plt.subplots(2, 2, figsize=(12.5, 10.4))
fig.subplots_adjust(hspace=0.62, wspace=0.28)

# ---------------------------------------------------------------- (a) plano vertical
a = axs[0, 0]
rr = np.linalg.norm(P0, axis=1).mean()
ang = np.radians(np.linspace(0, 90, 200))
a.plot(rr * np.sin(ang), rr * np.cos(ang), color=P.C["grid"], lw=1, zorder=0)
for thc in p["theta_fuera"]:
    if abs(thc - th_base) < 0.1:
        continue
    F2 = geometry.familia_C(P0, "F2", thc)
    a.plot(np.hypot(F2[:, 0], F2[:, 1]).mean(), F2[:, 2].mean(), "o", color=C2, ms=5, alpha=0.45)
    a.annotate(f"{thc:.0f}°", (rr * np.sin(np.radians(thc)), rr * np.cos(np.radians(thc))),
               xytext=(6, 4), textcoords="offset points", fontsize=7.5, color=CG)
a.plot(np.hypot(P0[:, 0], P0[:, 1]), P0[:, 2], color=C1, lw=3.0, label=f"hélice base = F1 en $\\theta_c$={th_base:.1f}°")
F1b = geometry.familia_C(P0, "F1", thB)
a.plot(np.hypot(F1b[:, 0], F1b[:, 1]), F1b[:, 2], color=P.C["s7"], lw=3.0,
       label=f"F1 en $\\theta_c=\\theta_{{pB}}$ (±7° de dispersión)")
for d in p["errores_deg"]:
    F2 = geometry.familia_C(P0, "F2", thB + d)
    a.plot(np.hypot(F2[:, 0], F2[:, 1]).mean(), F2[:, 2].mean(), "o", color=C2, ms=4 if d else 8,
           mfc=C2 if d == 0 else "white", label=None)
a.plot([], [], "o", color=C2, ms=8, label="F2 (cono): un único punto $(\\rho, z)$ por $\\theta_c$")
a.plot([], [], "o", color=C2, ms=4, mfc="white", label="F2 en $\\theta_{pB}+\\delta$, $\\delta$=±0.5, ±1, ±2, ±5°")
for thc, col, lab in [(th_base, CG, None), (thB, CB, f"$\\theta_{{pB}}$(S1) = {thB:.2f}°")]:
    s = np.radians(thc)
    a.plot([0, 200 * np.sin(s)], [0, 200 * np.cos(s)], color=col, lw=0.9, ls=":", label=lab)
a.axhline(0, color=P.C["ink"], lw=0.8)
a.plot(0, -2, "v", color=P.C["ink"], ms=6)
a.text(4, -12, "superficie z = 0, blanco en (0,0,−2) m", fontsize=7.5, color=CG)
a.set_xlim(0, 200); a.set_ylim(-20, 200); a.set_aspect("equal")
a.set_xlabel("radio horizontal $\\rho$ [m]"); a.set_ylabel("altura z [m]")
a.set_title("(a) Plano vertical: todas las trayectorias están a la misma distancia\n"
            f"del centro (≈{rr:.0f} m); sólo cambia la incidencia", fontsize=9)
a.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)

# ---------------------------------------------------------------- (b) incidencia por pulso
b = axs[0, 1]
b.plot(t, geometry.incidence_deg(P0), color=C1, lw=2.4, label=f"hélice base (F1, $\\theta_c$={th_base:.1f}°)")
b.plot(t, geometry.incidence_deg(F1b), color=P.C["s7"], lw=2.4, label="F1, $\\theta_c=\\theta_{pB}$")
for d in p["errores_deg"]:
    F2 = geometry.familia_C(P0, "F2", thB + d)
    b.plot(t, geometry.incidence_deg(F2), color=C2, lw=2.4 if d == 0 else 0.9, alpha=1 if d == 0 else 0.7,
           label="F2, $\\theta_c=\\theta_{pB}$" if d == 0 else ("F2, $\\theta_{pB}+\\delta$" if d == p["errores_deg"][0] else None))
    if d > 0:
        b.text(t[-1] + 0.2, thB + d, f"±{d:g}°", fontsize=6.5, va="center", color=CG)
b.axhline(thB, color=CB, lw=0.8, ls=":")
b.set_xlim(0, t[-1] + 1.6)
b.set_xlabel("t [s]"); b.set_ylabel("incidencia del Tx en el centro $\\theta_{Tx}(t)$ [°]")
b.set_title("(b) F1 barre ±7° alrededor de $\\theta_c$; F2 la mantiene constante.\n"
            "El error angular $\\delta$ desplaza todo el cono: $\\theta_c=\\theta_{pB}+\\delta$", fontsize=9)
b.legend(fontsize=7, loc="lower right")

# ---------------------------------------------------------------- (c) planta F2 + punto especular
c = axs[1, 0]
F2b = geometry.familia_C(P0, "F2", thB)
k = 0
c.plot(F2b[:, 0], F2b[:, 1], color=C1, lw=1.8, label="Tx (F2 en $\\theta_{pB}$)")
c.plot(*F2b[k, :2], "o", color=C1, ms=8)
for dphi, mk in [(90, "s"), (135, "D"), (180, "^")]:
    Prx = geometry.rx_from_tx(F2b, dphi)
    _, info = scene.specular_echo(sysm, F2b, Prx)
    Q = info["Q"][k]
    c.plot(*Prx[k, :2], mk, color=C2, ms=7)
    c.plot([F2b[k, 0], Q[0], Prx[k, 0]], [F2b[k, 1], Q[1], Prx[k, 1]], color=CG, lw=0.7, ls="--")
    c.plot(*Q[:2], "x", color=CB, ms=8, mew=2)
    c.annotate(f"Rx {dphi}°", Prx[k, :2], xytext=(6, -10), textcoords="offset points", fontsize=7.5)
    c.annotate(f"$\\theta_{{esp}}$={np.degrees(info['theta'][k]):.1f}°", Q[:2], xytext=(6, 4),
               textcoords="offset points", fontsize=7.5, color=CB)
c.plot(*F2b[k, :2], "x", color=CB, ms=0)
c.plot([], [], "x", color=CB, ms=8, mew=2, label="punto especular (Tx–superficie–Rx)")
c.plot([], [], "s", color=C2, ms=6, label="Rx a $\\Delta\\phi$ = 90°, 135°, 180° (Ref1)")
h = p["rugosidad"]["patch_half"]
c.add_patch(plt.Rectangle((-h, -h), 2 * h, 2 * h, fill=False, ec=P.C["s6"], lw=1.2))
c.text(h + 4, -h - 2, "parche difuso\n±12 m", fontsize=7, color=P.C["s6"], va="top")
c.annotate("mono: especular en el nadir\ndel Tx ($\\theta_{esp}$=0)", F2b[k, :2], xytext=(-125, -45),
           textcoords="offset points", fontsize=7.5, arrowprops=dict(arrowstyle="->", color=CG, lw=0.7))
c.set_aspect("equal"); c.set_xlim(-200, 200); c.set_ylim(-200, 200)
c.set_xlabel("x [m]"); c.set_ylabel("y [m]")
c.set_title(f"(c) Planta de F2 en $\\theta_{{pB}}$ (altura {F2b[:, 2].mean():.1f} m, radio "
            f"{np.hypot(F2b[:, 0], F2b[:, 1]).mean():.0f} m; pulso 0)\n"
            "y punto especular de la superficie para cada separación", fontsize=9)
c.legend(fontsize=7, loc="lower left")

# ---------------------------------------------------------------- (d) incidencia en el punto especular
d_ = axs[1, 1]
dph = np.linspace(0, 180, 181)
for thc, col, ls in [(th_base, C1, "-"), (thB, C2, "-")]:
    th_sp = np.degrees(np.arctan(np.tan(np.radians(thc)) * np.sin(np.radians(dph) / 2)))
    d_.plot(dph, th_sp, color=col, ls=ls, lw=2, label=f"F2 con $\\theta_c$ = {thc:.1f}°: "
            "$\\arctan(\\tan\\theta_c\\,\\sin\\frac{\\Delta\\phi}{2})$")
    Pc = geometry.familia_C(P0, "F2", thc)
    for dphi in [0, 45, 90, 135, 180]:
        Prx = Pc.copy() if dphi == 0 else geometry.rx_from_tx(Pc, dphi)
        _, info = scene.specular_echo(sysm, Pc, Prx)
        d_.plot(dphi, np.degrees(info["theta"]).mean(), "o", color=col, mfc="white", ms=6)
d_.plot([], [], "o", color=CG, mfc="white", ms=6, label="simulación (scene.specular_echo)")
d_.axhline(thB, color=CB, lw=0.9, ls=":", label=f"$\\theta_{{pB}}$(S1) = {thB:.2f}°")
d_.set_xlabel("separación $\\Delta\\phi$ [°]"); d_.set_ylabel("incidencia en el punto especular $\\theta_{esp}$ [°]")
d_.set_xticks(range(0, 181, 30)); d_.set_ylim(0, 80)
d_.set_title("(d) La incidencia del Tx ($\\theta_c$) sólo coincide con la del punto\n"
             "especular en $\\Delta\\phi$=180° (Ref1); a menor $\\Delta\\phi$ la especular es más vertical", fontsize=9)
d_.legend(fontsize=7, loc="upper left")

fig.suptitle("Exp. C: geometría de las familias F1 (hélice girada) y F2 (cono de incidencia constante)", y=0.965)
P.save(fig, "C_geometria", out)
print("ok", round(thB, 3), round(th_base, 2))
