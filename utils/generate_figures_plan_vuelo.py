"""
generate_figures_plan_vuelo.py
Genera todas las figuras para explicacion_plan_de_vuelo.tex
Guarda en: hypotheses/figures/pv_fig*.pdf

Ejecutar desde la raíz del proyecto:
    python utils/generate_figures_plan_vuelo.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch, Arc, FancyBboxPatch, Wedge
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import os

# ── directorios ──────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(SCRIPT_DIR, '..', 'hypotheses', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# ── estilo global ─────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':      'DejaVu Serif',
    'font.size':        11,
    'axes.labelsize':   12,
    'axes.titlesize':   13,
    'xtick.labelsize':  10,
    'ytick.labelsize':  10,
    'lines.linewidth':  2.0,
    'grid.alpha':       0.35,
    'figure.dpi':       150,
    'savefig.dpi':      200,
    'savefig.bbox':     'tight',
    'savefig.pad_inches': 0.1,
})

C = {
    'sky':    '#2c7bb6',
    'soil':   '#8c510a',
    'target': '#d7191c',
    'ray':    '#fdae61',
    'green':  '#1a9641',
    'purple': '#762a83',
    'gray':   '#636363',
    'orange': '#e08214',
    'teal':   '#018571',
    'tx':     '#2166ac',
    'rx':     '#d6604d',
    'brew':   '#f4a582',
}

def save(name):
    path = os.path.join(FIG_DIR, name)
    plt.savefig(path, bbox_inches='tight', pad_inches=0.12)
    plt.close()
    print(f'  Saved: {name}')

# ─────────────────────────────────────────────────────────────────────────────
# Física de Fresnel TM
# ─────────────────────────────────────────────────────────────────────────────
def fresnel_tm(theta1_deg, n1, n2):
    """Coeficientes de Fresnel TM. Devuelve r, t, R, T (potencia)."""
    th1 = np.radians(theta1_deg)
    sin_th2 = n1 / n2 * np.sin(th1)
    sin_th2 = np.clip(sin_th2, -1, 1)
    th2 = np.arcsin(sin_th2)
    cos1, cos2 = np.cos(th1), np.cos(th2)
    denom = n2 * cos1 + n1 * cos2
    r = np.where(np.abs(denom) > 1e-12,
                 (n2 * cos1 - n1 * cos2) / denom, 1.0)
    t = np.where(np.abs(denom) > 1e-12,
                 (2 * n1 * cos1) / denom, 0.0)
    R_pow = r ** 2
    T_pow = np.where(np.abs(cos1) > 1e-12,
                     (n2 * cos2) / (n1 * cos1) * t ** 2, 0.0)
    T_pow = np.clip(T_pow, 0, 1)
    return r, t, R_pow, T_pow

# ─────────────────────────────────────────────────────────────────────────────
# FIG 1: Coeficientes de Fresnel TM
# ─────────────────────────────────────────────────────────────────────────────
def fig01_fresnel_tm():
    n1, n2 = 1.0, 2.0
    theta = np.linspace(0, 90, 1800)
    theta_crit = np.degrees(np.arcsin(1.0))  # no hay ángulo crítico n1<n2
    theta_B = np.degrees(np.arctan(n2 / n1))  # ≈63.43°

    r, t, R, T = fresnel_tm(theta, n1, n2)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # ── Panel izquierdo: transmitancia y reflectancia de POTENCIA ──
    ax = axes[0]
    ax.plot(theta, T, color=C['sky'], lw=2.5, label='$T$ (transmitancia de potencia)')
    ax.plot(theta, R, color=C['target'], lw=2.5, label='$R$ (reflectancia de potencia)')
    ax.plot(theta, T + R, 'k--', lw=1.2, alpha=0.5, label='$T+R$ (conservación)')
    ax.axvline(theta_B, color=C['orange'], lw=1.8, ls=':', label=f'Brewster $\\theta_B={theta_B:.1f}°$')
    ax.fill_betweenx([0, 1.05], theta_B - 3, theta_B + 3, alpha=0.12, color=C['orange'])
    ax.annotate('$T=1,\\ R=0$\nen Brewster',
                xy=(theta_B, 1.00), xytext=(theta_B + 12, 0.82),
                fontsize=10, color=C['orange'],
                arrowprops=dict(arrowstyle='->', color=C['orange'], lw=1.5))
    ax.set_xlabel('Ángulo de incidencia $\\theta_1$ [°]')
    ax.set_ylabel('Potencia transmitida / reflejada')
    ax.set_title('(a) Transmitancia y Reflectancia de Potencia — TM')
    ax.set_xlim(0, 90)
    ax.set_ylim(-0.02, 1.10)
    ax.legend(loc='lower left', fontsize=9)
    ax.grid(True)

    # ── Panel derecho: coeficientes de CAMPO ──
    ax = axes[1]
    ax.plot(theta, np.abs(t), color=C['sky'], lw=2.5, label='$|t_{TM}|$ (transmisión campo)')
    ax.plot(theta, np.abs(r), color=C['target'], lw=2.5, label='$|r_{TM}|$ (reflexión campo)')
    ax.axvline(theta_B, color=C['orange'], lw=1.8, ls=':',
               label=f'Brewster $\\theta_B={theta_B:.1f}°$')
    ax.annotate('$r_{TM}=0$', xy=(theta_B, 0.0), xytext=(theta_B + 10, 0.18),
                fontsize=10, color=C['orange'],
                arrowprops=dict(arrowstyle='->', color=C['orange'], lw=1.5))
    ax.set_xlabel('Ángulo de incidencia $\\theta_1$ [°]')
    ax.set_ylabel('Coeficiente de campo eléctrico')
    ax.set_title('(b) Coeficientes de Campo — TM')
    ax.set_xlim(0, 90)
    ax.set_ylim(-0.02, 1.55)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True)

    plt.suptitle(f'Coeficientes de Fresnel para polarización TM ($n_1={n1}$, $n_2={n2}$, $\\varepsilon_r=4$)',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    save('pv_fig01_fresnel_tm.pdf')
    print(f'    theta_B = {theta_B:.2f} deg,  T(theta_B) ~= {np.interp(theta_B, theta, T):.4f}')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 2: NO GENERADA POR CÓDIGO — Geometría Fermat 3D (placeholder)
# pv_fig02 se omite; se usa tikz en el .tex directamente
# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
# FIG 3: Función de Snell 1D y convergencia Newton-Raphson
# ─────────────────────────────────────────────────────────────────────────────
def snell_f_and_deriv(u, D, sz, tz, n1, n2):
    v = D - u
    d1 = np.sqrt(u**2 + sz**2)
    d2 = np.sqrt(v**2  + tz**2)
    f  = n1 * u / d1 - n2 * v / d2
    fp = n1 * sz**2 / d1**3 + n2 * tz**2 / d2**3
    return f, fp

def fermat_newton(sz, tz, D, n1, n2, u0=None, max_iter=12):
    if u0 is None:
        u0 = D * (sz / n1) / (sz / n1 + tz / n2)
    u = np.clip(u0, 0, D)
    traj = [u]
    for _ in range(max_iter):
        f, fp = snell_f_and_deriv(u, D, sz, tz, n1, n2)
        du = -f / fp
        u = np.clip(u + du, 0, D)
        traj.append(u)
        if abs(du) < 1e-10 * D:
            break
    return u, traj

def fig03_snell_1d():
    n1, n2 = 1.0, 2.0
    sz, tz, D = 100.0, 5.0, 80.0

    u_arr = np.linspace(0, D, 800)
    f_arr = np.array([snell_f_and_deriv(u, D, sz, tz, n1, n2)[0] for u in u_arr])

    u_star, traj = fermat_newton(sz, tz, D, n1, n2)
    traj = np.array(traj)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # ── Panel izquierdo: f(u) ──
    ax = axes[0]
    ax.plot(u_arr, f_arr, color=C['sky'], lw=2.5)
    ax.axhline(0, color='k', lw=0.8, ls='--')
    ax.axvline(u_star, color=C['green'], lw=1.8, ls=':', label=f'$u^*={u_star:.2f}$ m')
    ax.fill_between(u_arr, 0, f_arr, where=(f_arr < 0), alpha=0.12, color=C['target'], label='$f<0$')
    ax.fill_between(u_arr, 0, f_arr, where=(f_arr > 0), alpha=0.12, color=C['sky'], label='$f>0$')
    ax.plot(u_star, 0, 'o', ms=9, color=C['green'], zorder=5, label='Raíz única $f(u^*)=0$')
    ax.set_xlabel('$u$ [m]  (posición del punto de refracción en la dirección horizontal)')
    ax.set_ylabel("$f(u) = n_1 u/d_1 - n_2(D{-}u)/d_2$")
    ax.set_title('(a) Función de Snell 1D — monótona creciente')
    ax.legend(fontsize=9)
    ax.grid(True)
    params_str = (f'$n_1={n1}$, $n_2={n2}$\n'
                  f'$s_z={sz}$ m, $t_z={tz}$ m, $D={D}$ m')
    ax.text(0.04, 0.96, params_str, transform=ax.transAxes,
            fontsize=9, va='top', bbox=dict(boxstyle='round', fc='white', ec='gray', alpha=0.8))

    # ── Panel derecho: iteraciones Newton-Raphson ──
    ax = axes[1]
    ax.plot(u_arr, f_arr, color=C['sky'], lw=2.5, label='$f(u)$', zorder=2)
    ax.axhline(0, color='k', lw=0.8, ls='--')

    colors_nr = plt.cm.Reds(np.linspace(0.4, 0.9, len(traj) - 1))
    for k in range(min(5, len(traj) - 1)):
        u_k = traj[k]
        f_k, fp_k = snell_f_and_deriv(u_k, D, sz, tz, n1, n2)
        # Tangente en u_k
        u_tan = np.array([u_k - 8, u_k + 8])
        f_tan = f_k + fp_k * (u_tan - u_k)
        u_next = u_k - f_k / fp_k
        # Vertical de bajada
        ax.plot([u_k, u_k], [0, f_k], ':', color=colors_nr[k], lw=1.5, zorder=3)
        # Tangente
        ax.plot(u_tan, f_tan, '-', color=colors_nr[k], lw=1.5, zorder=3)
        # Punto k
        ax.plot(u_k, f_k, 'o', ms=7, color=colors_nr[k], zorder=5,
                label=f'$u_{k}$' if k < 3 else '')
        # Flecha hacia nuevo punto
        ax.annotate('', xy=(u_next, 0), xytext=(u_k, 0),
                    arrowprops=dict(arrowstyle='->', color=colors_nr[k], lw=1.5))

    ax.plot(u_star, 0, '*', ms=14, color=C['green'], zorder=6, label=f'$u^*={u_star:.2f}$ m')

    err = np.abs(np.array(traj) - u_star)
    err_str = '\n'.join([f'  $|e_{k}|={err[k]:.4f}$ m' for k in range(min(4, len(traj)))])
    ax.text(0.04, 0.96, 'Errores:\n' + err_str, transform=ax.transAxes,
            fontsize=9, va='top', family='monospace',
            bbox=dict(boxstyle='round', fc='lightyellow', ec='gray', alpha=0.9))

    ax.set_xlabel('$u$ [m]')
    ax.set_ylabel("$f(u)$")
    ax.set_title('(b) Convergencia Newton-Raphson (4 iteraciones)')
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(True)

    plt.suptitle('Solución de Snell 1D por Newton-Raphson para interfaz plana horizontal',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    save('pv_fig03_snell_1d.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 4: Ecuación de radar biestática — diagrama geométrico
# ─────────────────────────────────────────────────────────────────────────────
def fig04_bistatic_radar():
    fig, ax = plt.subplots(figsize=(11, 6))

    # Interfaz
    ax.axhline(0, color=C['soil'], lw=2, ls='--', label='Interfaz $z=0$')
    ax.fill_between([-0.5, 10.5], [-3, -3], [0, 0], color=C['soil'], alpha=0.10)
    ax.fill_between([-0.5, 10.5], [0, 0], [5.5, 5.5], color='lightcyan', alpha=0.25)
    ax.text(0.15, 0.15, 'Aire  $n_1=1$', fontsize=11, color=C['sky'], transform=ax.transAxes)
    ax.text(0.15, 0.08, 'Suelo  $n_2=2$', fontsize=11, color=C['soil'], transform=ax.transAxes)

    # Posiciones
    TX = np.array([1.5, 4.5])
    Q1 = np.array([3.5, 0.0])
    P  = np.array([5.0, -2.2])
    Q2 = np.array([6.8, 0.0])
    RX = np.array([8.5, 3.5])

    def arrow(ax, A, B, col, lw=2, style='->', label=None, **kw):
        dx, dy = B - A
        ax.annotate('', xy=B, xytext=A,
                    arrowprops=dict(arrowstyle=style, color=col, lw=lw, **kw))
        if label:
            mid = (A + B) / 2
            ax.text(mid[0], mid[1] + 0.2, label, ha='center', fontsize=10, color=col,
                    fontweight='bold')

    # Rayos con etiquetas de distancia
    arrow(ax, TX, Q1, C['tx'], lw=2.2)
    arrow(ax, Q1, P,  C['tx'], lw=2.2)
    arrow(ax, P,  Q2, C['rx'], lw=2.2)
    arrow(ax, Q2, RX, C['rx'], lw=2.2)

    # Etiquetas de segmentos
    mid_TQ1 = (TX + Q1) / 2
    ax.text(mid_TQ1[0] - 0.3, mid_TQ1[1] + 0.15, '$R_1$', fontsize=12, color=C['tx'], fontweight='bold')
    mid_Q1P = (Q1 + P) / 2
    ax.text(mid_Q1P[0] - 0.4, mid_Q1P[1], '$R_2$', fontsize=12, color=C['tx'], fontweight='bold')
    mid_PQ2 = (P + Q2) / 2
    ax.text(mid_PQ2[0] + 0.12, mid_PQ2[1], '$R_3$', fontsize=12, color=C['rx'], fontweight='bold')
    mid_Q2R = (Q2 + RX) / 2
    ax.text(mid_Q2R[0] + 0.1, mid_Q2R[1] + 0.15, '$R_4$', fontsize=12, color=C['rx'], fontweight='bold')

    # Puntos
    for pt, lab, col, mk in [(TX, 'TX', C['tx'], 's'), (RX, 'RX', C['rx'], 's'),
                               (Q1, '$Q_1^*$', C['orange'], 'D'), (Q2, '$Q_2^*$', C['orange'], 'D'),
                               (P,  'Target $P$', C['target'], '*')]:
        ax.plot(*pt, mk, ms=10 if mk == '*' else 8, color=col, zorder=5)
        offset = [0.12, 0.15]
        if lab in ('RX', '$Q_2^*$'):
            offset[0] = 0.12
        ax.text(pt[0] + offset[0], pt[1] + offset[1], lab, fontsize=11, color=col, fontweight='bold')

    # Llaves con R_T y R_R
    ax.annotate('', xy=(Q1[0], -0.6), xytext=(TX[0], -0.6),
                arrowprops=dict(arrowstyle='<->', color=C['tx'], lw=1.5))
    ax.annotate('', xy=(P[0], -0.6), xytext=(Q1[0], -0.6),
                arrowprops=dict(arrowstyle='<->', color=C['tx'], lw=1.5))
    ax.text((TX[0] + P[0]) / 2, -1.0,
            r'$R_T = R_1 + R_2$', ha='center', fontsize=11, color=C['tx'], fontweight='bold')

    ax.annotate('', xy=(Q2[0], -0.6), xytext=(P[0], -0.6),
                arrowprops=dict(arrowstyle='<->', color=C['rx'], lw=1.5))
    ax.annotate('', xy=(RX[0], -0.6), xytext=(Q2[0], -0.6),
                arrowprops=dict(arrowstyle='<->', color=C['rx'], lw=1.5))
    ax.text((P[0] + RX[0]) / 2, -1.0,
            r'$R_R = R_3 + R_4$', ha='center', fontsize=11, color=C['rx'], fontweight='bold')

    # Transmitancias en interfaz
    ax.text(Q1[0] + 0.1, 0.25, '$T_1$', fontsize=12, color=C['orange'], fontweight='bold')
    ax.text(Q2[0] + 0.1, 0.25, '$T_2$', fontsize=12, color=C['orange'], fontweight='bold')

    # Fórmula radar
    formula = (r'$P_r = \dfrac{P_t\,G_T\,G_R\,\sigma\,T_1\,T_2\,\lambda^2}'
               r'{(4\pi)^3\,R_T^2\,R_R^2}$')
    ax.text(0.5, 0.96, formula, transform=ax.transAxes, ha='center', va='top',
            fontsize=13, bbox=dict(boxstyle='round', fc='lightyellow', ec='orange', lw=1.5))

    ax.set_xlim(-0.3, 10.5)
    ax.set_ylim(-1.6, 5.8)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Ecuación de Radar Biestático con Refracción en Dos Medios',
                 fontsize=13, fontweight='bold', pad=12)
    save('pv_fig04_bistatic_radar.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 5: Mapa de costo J(xR, yR) y perfil 1D
# ─────────────────────────────────────────────────────────────────────────────
def compute_cost_map():
    """Versión simplificada de J para visualización (sin solver Snell completo)."""
    n1, n2 = 1.0, 2.0
    Pt  = 1e5
    lam = 0.75
    sig = 1.0
    G   = 4.0
    Rx_z = 100.0
    alpha_res = 0.3

    # Tx fijo
    TX = np.array([160.0, 0.0, 100.0])
    # Target centroid
    tg_c = np.array([0.0, 0.0, -5.0])
    theta_B = np.arctan(n2 / n1)

    def approx_cost(xR, yR):
        RX = np.array([xR, yR, Rx_z])
        # distancias aproximadas (sin refracción exacta para el mapa rápido)
        R_T = np.linalg.norm(TX - tg_c)
        R_R = np.linalg.norm(RX - tg_c)

        # Ángulo de incidencia aprox. para Tx (desde vertical)
        dxy_TX = np.linalg.norm(TX[:2] - tg_c[:2])
        th_TX = np.arctan2(dxy_TX, TX[2] + abs(tg_c[2]))
        _, _, _, T1 = fresnel_tm(np.degrees(th_TX), n1, n2)

        # Ángulo de incidencia aprox. para Rx
        dxy_RX = np.linalg.norm(RX[:2] - tg_c[:2])
        th_RX = np.arctan2(dxy_RX, Rx_z + abs(tg_c[2]))
        _, _, _, T2 = fresnel_tm(np.degrees(th_RX), n1, n2)

        Pr = (Pt * G * G * sig * T1 * T2 * lam**2) / ((4*np.pi)**3 * R_T**2 * R_R**2)
        Pr = max(Pr, 1e-300)

        # Resolución aproximada
        rho_RX = np.sqrt(xR**2 + yR**2)
        psi0 = np.arctan2(rho_RX, Rx_z)
        psi0 = max(psi0, 0.01)

        dxy = 0.60 * lam / (np.pi * np.sin(psi0) + 1e-9)

        # kspace vertical nominal
        B_hz = 50e6
        c  = 3e8
        cos_tht = np.sqrt(max(1 - (np.sin(psi0)/n2)**2, 0))
        B_perp = 47.17  # m
        R0 = np.sqrt(rho_RX**2 + Rx_z**2) + 1
        Wz = n2*B_hz*cos_tht + c*B_perp*np.sin(psi0)*np.cos(psi0)/(lam*R0*n2*cos_tht+1e-9)
        dz = c / (2*Wz) if Wz > 0 else 100.0
        dz = min(dz, 100.0)

        log_Pr  = np.log10(Pr)
        log_res = np.log10(dxy * dz)

        J = -(1-alpha_res)*log_Pr + alpha_res*log_res
        return J, Pr, dxy, dz

    return approx_cost, TX, tg_c, theta_B

def fig05_cost_function():
    approx_cost, TX, tg_c, theta_B = compute_cost_map()

    # Grid 2D
    xs = np.linspace(-200, 200, 120)
    ys = np.linspace(-200, 200, 120)
    XX, YY = np.meshgrid(xs, ys)
    JJ = np.zeros_like(XX)
    for i in range(XX.shape[0]):
        for j in range(XX.shape[1]):
            JJ[i, j] = approx_cost(XX[i, j], YY[i, j])[0]

    # Perfil 1D a lo largo de la dirección opuesta al TX
    tx_dir = np.array([TX[0], TX[1]])
    opp_dir = -tx_dir / np.linalg.norm(tx_dir)
    rs = np.linspace(10, 200, 400)
    J_profile = []
    for r in rs:
        xR, yR = r * opp_dir
        J_profile.append(approx_cost(xR, yR)[0])
    J_profile = np.array(J_profile)

    # Brewster distance (approximate)
    d_brew = (100.0 + 5.0) / np.tan(theta_B)
    xB, yB = d_brew * opp_dir

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # Mapa 2D
    ax = axes[0]
    vmin, vmax = np.percentile(JJ, [5, 95])
    cm = ax.contourf(XX, YY, JJ, levels=30, cmap='RdYlBu_r', vmin=vmin, vmax=vmax)
    plt.colorbar(cm, ax=ax, label='Costo $J$')
    ax.contour(XX, YY, JJ, levels=15, colors='k', linewidths=0.3, alpha=0.4)

    # Tx y candidatos de inicio
    ax.plot(TX[0], TX[1], 's', ms=12, color=C['tx'], zorder=8, label='TX')
    ax.plot(tg_c[0], tg_c[1], '*', ms=14, color=C['target'], zorder=8, label='Centroide targets')

    # 5 candidatos
    d = np.linalg.norm(TX[:2] - tg_c[:2])
    opp = np.array([tg_c[0] - TX[0], tg_c[1] - TX[1]])
    opp = opp / np.linalg.norm(opp)
    perp = np.array([-opp[1], opp[0]])
    cands = [
        tg_c[:2] + opp * d * 0.8,
        tg_c[:2] + perp * d * 0.8,
        tg_c[:2] - perp * d * 0.8,
        tg_c[:2] + d_brew * np.array([np.cos(np.pi + np.arctan2(TX[1]-tg_c[1], TX[0]-tg_c[0])),
                                        np.sin(np.pi + np.arctan2(TX[1]-tg_c[1], TX[0]-tg_c[0]))]) * 0.9,
        tg_c[:2] + d_brew * np.array([np.cos(np.pi + np.arctan2(TX[1]-tg_c[1], TX[0]-tg_c[0])),
                                        np.sin(np.pi + np.arctan2(TX[1]-tg_c[1], TX[0]-tg_c[0]))]) * 1.1,
    ]
    labels = ['Cand. 1\n(opuesto)', 'Cand. 2\n(perp+)', 'Cand. 3\n(perp−)',
              'Cand. 4\n(Brewster)', 'Cand. 5\n(Brewster)']
    for k, (c_, lab) in enumerate(zip(cands, labels)):
        ax.plot(*c_, 'x', ms=10, mew=2.5, color='white', zorder=9)
        ax.text(c_[0] + 5, c_[1] + 5, f'{k+1}', fontsize=9, color='white',
                fontweight='bold', zorder=10)

    # Mínimo aproximado
    imin = np.unravel_index(np.argmin(JJ), JJ.shape)
    ax.plot(XX[imin], YY[imin], 'D', ms=10, color=C['green'], zorder=10, label='Mínimo global')
    ax.set_xlabel('$x_R$ [m]')
    ax.set_ylabel('$y_R$ [m]')
    ax.set_title('(a) Mapa de costo $J(x_R, y_R)$ — $\\alpha=0.3$')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_xlim(-200, 200)
    ax.set_ylim(-200, 200)
    ax.set_aspect('equal')

    # Perfil 1D
    ax = axes[1]
    ax.plot(rs, J_profile, color=C['sky'], lw=2.5)
    idx_min = np.argmin(J_profile)
    ax.plot(rs[idx_min], J_profile[idx_min], 'D', ms=10, color=C['green'],
            zorder=5, label=f'Mínimo $r^*={rs[idx_min]:.0f}$ m')
    ax.axvline(d_brew, color=C['orange'], ls=':', lw=1.8,
               label=f'$d_B={d_brew:.0f}$ m (Brewster)')
    # Candidatos sobre el perfil
    for k, r_c in enumerate([d * 0.8, d_brew * 0.9, d_brew * 1.1]):
        j_c = approx_cost(r_c * opp[0], r_c * opp[1])[0]
        ax.plot(r_c, j_c, 'x', ms=10, mew=2.5, color=C['target'], zorder=7)
    ax.set_xlabel('Radio del RX desde el origen [m] (dirección opuesta al TX)')
    ax.set_ylabel('Costo $J$')
    ax.set_title('(b) Perfil de $J$ en la dirección opuesta al TX')
    ax.legend(fontsize=9)
    ax.grid(True)
    ax.text(0.04, 0.96,
            f'$\\alpha={0.3}$, TX en $({TX[0]:.0f}, {TX[1]:.0f}, {TX[2]:.0f})$ m',
            transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(boxstyle='round', fc='white', ec='gray', alpha=0.8))

    plt.suptitle('Función de costo multi-objetivo $J = -(1-\\alpha)\\log_{10}P_r + \\alpha\\log_{10}(\\delta_{xy}\\delta_z)$',
                 fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    save('pv_fig05_cost_function.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 6: Candidatos multi-inicio (vista cenital)
# ─────────────────────────────────────────────────────────────────────────────
def fig06_multistart():
    fig, ax = plt.subplots(figsize=(8, 8))

    TX = np.array([160.0, 0.0])
    tg_c = np.array([0.0, 0.0])
    n1, n2 = 1.0, 2.0
    theta_B = np.arctan(n2 / n1)
    d = np.linalg.norm(TX - tg_c)
    d_brew = (100.0 + 5.0) / np.tan(theta_B)

    opp = (tg_c - TX) / np.linalg.norm(tg_c - TX)
    perp = np.array([-opp[1], opp[0]])

    az_opp = np.arctan2(opp[1], opp[0])

    cands = {
        '1\nOpuesto': tg_c + opp * d * 0.8,
        '2\nPerp. +': tg_c + perp * d * 0.8,
        '3\nPerp. −': tg_c - perp * d * 0.8,
        '4\nBrewster': tg_c + d_brew * np.array([np.cos(az_opp), np.sin(az_opp)]),
        '5\nBrewster\n(jitter)': tg_c + d_brew * 1.18 * np.array(
            [np.cos(az_opp + 0.2), np.sin(az_opp + 0.2)]),
    }

    # Círculo de búsqueda
    theta_circ = np.linspace(0, 2*np.pi, 300)
    ax.plot(200*np.cos(theta_circ), 200*np.sin(theta_circ), '--',
            color='gray', lw=1.0, alpha=0.5, label='Región de búsqueda ($r=200$ m)')

    # Círculo de Brewster
    ax.plot(d_brew*np.cos(theta_circ), d_brew*np.sin(theta_circ),
            ':', color=C['orange'], lw=1.5, alpha=0.7, label=f'Radio Brewster $d_B={d_brew:.0f}$ m')

    # Flechas desde target hacia candidatos
    colors_cands = [C['sky'], C['purple'], C['teal'], C['orange'], C['brew']]
    for (lab, pt), col in zip(cands.items(), colors_cands):
        ax.annotate('', xy=pt, xytext=tg_c,
                    arrowprops=dict(arrowstyle='->', color=col, lw=1.5, alpha=0.6))
        ax.plot(*pt, 'x', ms=12, mew=3.0, color=col, zorder=8)
        ax.text(pt[0] + 6, pt[1] + 6, lab, fontsize=9.5, color=col,
                fontweight='bold', ha='left', va='bottom')

    # TX y centroide
    ax.plot(*TX, 's', ms=13, color=C['tx'], zorder=9, label='TX')
    ax.text(TX[0]+6, TX[1]+6, 'TX', fontsize=11, color=C['tx'], fontweight='bold')

    # Targets (pequeña nube)
    np.random.seed(42)
    for _ in range(25):
        pt = tg_c + np.random.uniform(-20, 20, 2)
        ax.plot(*pt, 'x', ms=5, mew=1.0, color=C['target'], alpha=0.5, zorder=4)
    ax.plot(*tg_c, '*', ms=16, color=C['target'], zorder=9, label='Centroide targets')

    # Línea TX → centroide
    ax.plot([TX[0], tg_c[0]], [TX[1], tg_c[1]], ':', color=C['tx'], lw=1.2, alpha=0.5)

    ax.set_xlim(-230, 230)
    ax.set_ylim(-230, 230)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('$x$ [m]', fontsize=12)
    ax.set_ylabel('$y$ [m]', fontsize=12)
    ax.set_title('Vista Cenital: 5 Candidatos de Inicio para la Búsqueda del RX\n'
                 f'($n_1={n1}$, $n_2={n2}$, $\\theta_B={np.degrees(theta_B):.1f}°$)',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='lower left')
    save('pv_fig06_multistart.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 7: Frontera de Pareto
# ─────────────────────────────────────────────────────────────────────────────
def fig07_pareto():
    approx_cost, TX, tg_c, theta_B = compute_cost_map()
    n1, n2 = 1.0, 2.0
    c_light = 3e8
    lam = 0.75
    B_hz = 50e6
    d_brew = (100.0 + abs(tg_c[2])) / np.tan(theta_B)
    opp_dir = -(TX[:2] / np.linalg.norm(TX[:2]))

    alphas = np.linspace(0, 1, 41)
    Pr_pareto = []
    res_pareto = []

    # Para cada alpha, evaluar en el punto de Brewster (aproximación a la solución óptima)
    for alpha in alphas:
        # Barrido en radio a lo largo de la dirección opuesta
        rs = np.linspace(20, 200, 200)
        best_J = np.inf
        best_r = d_brew
        for r in rs:
            xR, yR = r * opp_dir
            J, Pr, dxy, dz = approx_cost(xR, yR)
            J_alpha = -(1-alpha)*np.log10(Pr) + alpha*np.log10(dxy*dz)
            if J_alpha < best_J:
                best_J = J_alpha
                best_r = r
                best_Pr = Pr
                best_res = dxy * dz
        Pr_pareto.append(best_Pr)
        res_pareto.append(best_res)

    Pr_dBm = 10 * np.log10(np.array(Pr_pareto) * 1000)
    res_cm2 = np.array(res_pareto) * 1e4  # m² -> cm²

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # ── Panel izquierdo: Pareto en (Pr, resolución) ──
    ax = axes[0]
    sc = ax.scatter(Pr_dBm, res_cm2, c=alphas, cmap='plasma', s=60, zorder=5)
    ax.plot(Pr_dBm, res_cm2, '-', color='gray', lw=1.2, alpha=0.5, zorder=4)
    plt.colorbar(sc, ax=ax, label='$\\alpha$ (peso resolución)')

    # Marcadores alpha=0, 0.3, 0.87, 1
    for alpha_mark, lab in [(0, '$\\alpha=0$\n(solo potencia)'),
                             (0.3, '$\\alpha=0.3$\n(nominal)'),
                             (0.87, '$\\alpha^*\\approx0.87$'),
                             (1, '$\\alpha=1$\n(solo resolución)')]:
        idx = np.argmin(np.abs(alphas - alpha_mark))
        ax.plot(Pr_dBm[idx], res_cm2[idx], 'D', ms=10, zorder=8,
                color=plt.cm.plasma(alpha_mark))
        ax.annotate(lab, xy=(Pr_dBm[idx], res_cm2[idx]),
                    xytext=(Pr_dBm[idx] - 1.5, res_cm2[idx] + 0.5*(1+alpha_mark)),
                    fontsize=8.5, ha='right',
                    arrowprops=dict(arrowstyle='->', lw=1.0, color='gray'))

    ax.set_xlabel('Potencia recibida $P_r$ [dBm]', fontsize=11)
    ax.set_ylabel('Área de resolución $\\delta_{xy}\\cdot\\delta_z$ [cm$^2$]', fontsize=11)
    ax.set_title('(a) Frontera de Pareto: Potencia vs. Resolución', fontsize=12)
    ax.grid(True)
    ax.invert_xaxis()  # mayor dBm = mejor potencia = izquierda

    # ── Panel derecho: J* vs alpha ──
    ax = axes[1]
    J_power = [-np.log10(p) for p in Pr_pareto]
    J_res   = [np.log10(r) for r in res_pareto]
    J_total = [-(1-a)*jp + a*jr for a, jp, jr in zip(alphas, J_power, J_res)]

    ax.plot(alphas, J_power, color=C['sky'], lw=2, label='Contrib. potencia $-(1-\\alpha)\\log_{10}P_r$')
    ax.plot(alphas, J_res,   color=C['target'], lw=2, label='Contrib. resolución $\\alpha\\log_{10}(\\delta\\delta_z)$')
    ax.plot(alphas, J_total, 'k-', lw=2.5, label='Costo total $J^*$')
    ax.axvline(0.87, color=C['orange'], ls=':', lw=1.8, label='$\\alpha^*\\approx0.87$')
    ax.set_xlabel('Peso de resolución $\\alpha$', fontsize=11)
    ax.set_ylabel('Valor del costo en el óptimo', fontsize=11)
    ax.set_title('(b) Descomposición de $J^*$ vs. $\\alpha$', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True)

    plt.suptitle('Frontera de Pareto: equilibrio Potencia — Resolución 3D',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    save('pv_fig07_pareto.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 8: Geometría Brewster 2D (sección transversal)
# ─────────────────────────────────────────────────────────────────────────────
def fig08_brewster_geometry():
    fig, ax = plt.subplots(figsize=(10, 7))

    n1, n2 = 1.0, 2.0
    theta_B = np.arctan(n2 / n1)
    z_TX = 100.0
    z_tg = -5.0
    x_tg = 0.0

    # Posición Tx para que el rayo incida en Brewster
    d_B = (z_TX + abs(z_tg)) / np.tan(theta_B)  # distancia horizontal TX a Q1
    x_TX = -d_B
    x_Q1 = 0.0  # punto de refracción Tx→Target en x=0

    # Rx simétrico (Target→Rx incide también en Brewster)
    x_Q2 = x_tg  # coincide con x_tg cuando target en x=0
    x_RX = d_B
    z_RX = z_TX

    # Interfaz
    ax.axhline(0, color=C['soil'], lw=2.5, ls='--')
    ax.fill_between([-200, 200], [-30, -30], [0, 0], color=C['soil'], alpha=0.08)
    ax.fill_between([-200, 200], [0, 0], [140, 140], color='lightcyan', alpha=0.20)
    ax.text(-180, 120, 'Aire  ($n_1=1$)', fontsize=11, color=C['sky'])
    ax.text(-180, -22, 'Suelo  ($n_2=2$)', fontsize=11, color=C['soil'])

    # Rayos
    def draw_ray(ax, A, B, col, lw=2.5, ls='-', arrowfrac=0.55):
        mid = A + arrowfrac * (B - A)
        ax.annotate('', xy=mid + 0.001*(B-A), xytext=mid,
                    arrowprops=dict(arrowstyle='->', color=col, lw=lw))
        ax.plot([A[0], B[0]], [A[1], B[1]], ls, color=col, lw=lw)

    TX  = np.array([x_TX, z_TX])
    Q1  = np.array([x_Q1, 0.0])
    TG  = np.array([x_tg, z_tg])
    Q2  = np.array([x_Q2 + 0.0, 0.0])   # ligeramente desplazado para visibilidad
    RX  = np.array([x_RX, z_RX])

    draw_ray(ax, TX, Q1, C['tx'])
    draw_ray(ax, Q1, TG, C['tx'])
    draw_ray(ax, TG, Q2, C['rx'])
    draw_ray(ax, Q2, RX, C['rx'])

    # Normales
    for xn in [x_Q1, x_Q2]:
        ax.plot([xn, xn], [-20, 20], ':', color='gray', lw=1.0, alpha=0.7)

    # Ángulos de Brewster (arcos)
    for (xq, xA, zA), col in [((x_Q1, x_TX, z_TX), C['tx']), ((x_Q2, x_RX, z_RX), C['rx'])]:
        angle_deg = np.degrees(np.arctan2(abs(zA), abs(xA - xq)))
        arc = Arc((xq, 0), 30, 30, angle=90, theta1=0, theta2=angle_deg,
                  color=col, lw=1.5, ls='-')
        ax.add_patch(arc)
        mid_angle = np.radians(angle_deg / 2)
        ax.text(xq + 18*np.sin(mid_angle), 16*np.cos(mid_angle),
                f'$\\theta_B={np.degrees(theta_B):.1f}°$',
                fontsize=10, color=col, ha='center')

    # Puntos y etiquetas
    for pt, lab, col, mk in [(TX, 'TX', C['tx'], 's'), (RX, 'RX', C['rx'], 's'),
                               (Q1, '$Q_1^*$', C['orange'], 'D'), (Q2, '$Q_2^*$', C['orange'], 'D'),
                               (TG, 'Target $P$', C['target'], '*')]:
        ax.plot(*pt, mk, ms=11 if mk=='*' else 9, color=col, zorder=8)
        offset = np.array([5, 4])
        if lab in ('TX', '$Q_1^*$'):
            offset[0] = -45
        ax.text(pt[0]+offset[0], pt[1]+offset[1], lab, fontsize=11, color=col, fontweight='bold')

    # Cota de distancia Brewster
    ax.annotate('', xy=(x_Q1, -28), xytext=(x_TX, -28),
                arrowprops=dict(arrowstyle='<->', color=C['tx'], lw=1.5))
    ax.text((x_TX + x_Q1)/2, -35, f'$d_B = {int(d_B)}$ m', ha='center',
            fontsize=10, color=C['tx'])

    ax.set_xlim(-200, 200)
    ax.set_ylim(-45, 145)
    ax.set_xlabel('Posición horizontal [m]', fontsize=12)
    ax.set_ylabel('Altura / Profundidad [m]', fontsize=12)
    ax.set_title('Geometría de Brewster Biestática: Configuración Óptima en Sección Transversal\n'
                 f'($n_1={n1}$, $n_2={n2}$, $\\theta_B={np.degrees(theta_B):.1f}°$)',
                 fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Caja con la distancia de Brewster
    ax.text(0.02, 0.97,
            f'$d_B = (z_{{TX}}+|z_P|)/\\tan\\theta_B$\n'
            f'$= ({int(z_TX)}+{int(abs(z_tg))})/\\tan({np.degrees(theta_B):.1f}°)$\n'
            f'$\\approx {d_B:.0f}$ m',
            transform=ax.transAxes, fontsize=10, va='top',
            bbox=dict(boxstyle='round', fc='lightyellow', ec='orange', alpha=0.9))
    save('pv_fig08_brewster_geometry.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 9: Atenuación en suelo P-band vs X-band + sensibilidad de alpha
# ─────────────────────────────────────────────────────────────────────────────
def fig09_alpha_sensitivity():
    """Dos paneles: (a) atenuación en suelo vs freq, (b) sensibilidad de J a alpha."""

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # ── Panel izquierdo: Atenuación en suelo ──
    ax = axes[0]

    # Datos aproximados de atenuación en suelo (dB/m) para distintos suelos
    # Fuente: Ulaby et al., tablas de suelos típicos
    freqs_ghz = np.array([0.1, 0.2, 0.4, 0.8, 1.5, 3.0, 6.0, 10.0])
    # Suelo seco (low-loss): σ≈0.001 S/m, ε_r≈4
    atten_dry = np.array([0.02, 0.04, 0.08, 0.16, 0.30, 0.58, 1.1, 1.8])
    # Suelo húmedo (medium): σ≈0.01-0.1 S/m, ε_r≈8-15
    atten_wet = np.array([0.5, 0.8, 1.5, 3.0, 6.0, 12.0, 22.0, 38.0])
    # Suelo saturado (high-loss)
    atten_sat = np.array([2.0, 3.5, 7.0, 14.0, 28.0, 55.0, 100.0, 160.0])

    ax.semilogy(freqs_ghz, atten_dry, 'o-', color=C['green'], lw=2.2, ms=6,
                label='Suelo seco ($\\sigma\\approx 0.001$ S/m)')
    ax.semilogy(freqs_ghz, atten_wet, 's-', color=C['orange'], lw=2.2, ms=6,
                label='Suelo húmedo ($\\sigma\\approx 0.01$ S/m)')
    ax.semilogy(freqs_ghz, atten_sat, 'D-', color=C['target'], lw=2.2, ms=6,
                label='Suelo saturado ($\\sigma\\approx 0.1$ S/m)')

    # Límites de profundidad para detección a 5 m (necesita < 60 dB total = 12 dB/m @ 5m)
    ax.axhline(12.0, color='gray', ls='--', lw=1.5, alpha=0.7, label='Límite ($<12$ dB/m para 5 m)')

    # Marcadores de bandas
    ax.axvspan(0.3, 1.0, alpha=0.10, color='blue', label='Banda P/L (0.3–1 GHz)')
    ax.axvspan(8, 12, alpha=0.10, color='red', label='Banda X (8–12 GHz)')
    ax.axvline(0.4, color=C['sky'], ls=':', lw=2.0)
    ax.text(0.4, 0.3, '  400 MHz\n  (sistema)', fontsize=9.5, color=C['sky'], fontweight='bold')
    ax.axvline(10.0, color=C['target'], ls=':', lw=2.0)
    ax.text(10.2, 0.3, '10 GHz', fontsize=9.5, color=C['target'], fontweight='bold')

    ax.set_xlabel('Frecuencia [GHz]', fontsize=12)
    ax.set_ylabel('Atenuación en suelo [dB/m]', fontsize=12)
    ax.set_title('(a) Atenuación EM en suelo vs. frecuencia\n(razón para elegir banda P)', fontsize=11)
    ax.legend(fontsize=8.5, loc='lower right')
    ax.grid(True, which='both', alpha=0.35)
    ax.set_xlim(0.08, 12)

    # Anotación de pérdida total
    ax.text(0.02, 0.97,
            'Pérdida total ida-vuelta @ 5 m:\n'
            '  Banda P: $\\approx 2\\times 5\\times 0.08 = 0.8$ dB\n'
            '  Banda X: $\\approx 2\\times 5\\times 38 = 380$ dB (!)',
            transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(boxstyle='round', fc='lightyellow', ec='orange', alpha=0.9))

    # ── Panel derecho: sensibilidad de J a alpha ──
    ax = axes[1]
    approx_cost, TX, tg_c, theta_B = compute_cost_map()
    d_brew = (100.0 + 5.0) / np.tan(theta_B)
    opp_dir = -(TX[:2] / np.linalg.norm(TX[:2]))

    alphas = np.linspace(0, 1, 61)
    Pr_opt, res_opt = [], []
    r_opt_arr = []

    for alpha in alphas:
        rs = np.linspace(20, 200, 300)
        best_J = np.inf
        best_Pr = 1e-15
        best_res = 100.0
        best_r = d_brew
        for r in rs:
            xR, yR = r * opp_dir
            J, Pr, dxy, dz = approx_cost(xR, yR)
            J_a = -(1-alpha)*np.log10(Pr) + alpha*np.log10(dxy*dz)
            if J_a < best_J:
                best_J, best_Pr, best_res, best_r = J_a, Pr, dxy*dz, r
        Pr_opt.append(best_Pr)
        res_opt.append(best_res)
        r_opt_arr.append(best_r)

    Pr_dBm = 10*np.log10(np.array(Pr_opt)*1000)
    res_cm2 = np.array(res_opt) * 1e4

    ax2_r = ax.twinx()
    l1, = ax.plot(alphas, Pr_dBm,  color=C['sky'],    lw=2.5, label='$P_r^*$ [dBm] (eje izq.)')
    l2, = ax2_r.plot(alphas, res_cm2, color=C['target'], lw=2.5, ls='--',
                     label='$\\delta_{xy}\\delta_z^*$ [cm²] (eje der.)')
    ax.axvline(0.87, color=C['orange'], ls=':', lw=1.8, label='$\\alpha^*\\approx0.87$')
    ax.axvline(0.3,  color=C['green'],  ls=':', lw=1.5, label='$\\alpha=0.3$ (nominal)')

    ax.set_xlabel('Parámetro de peso $\\alpha$', fontsize=12)
    ax.set_ylabel('Potencia óptima $P_r^*$ [dBm]', fontsize=12, color=C['sky'])
    ax2_r.set_ylabel('Área resolución $\\delta_{xy}\\delta_z^*$ [cm²]', fontsize=12, color=C['target'])
    ax.set_title('(b) Sensibilidad de la solución óptima a $\\alpha$', fontsize=11)
    ax.tick_params(axis='y', labelcolor=C['sky'])
    ax2_r.tick_params(axis='y', labelcolor=C['target'])
    lines = [l1, l2] + [ax.axvline(0.87, color=C['orange'], ls=':', lw=1.8),
                         ax.axvline(0.3, color=C['green'], ls=':', lw=1.5)]
    labels = ['$P_r^*$', '$\\delta_{xy}\\delta_z^*$', '$\\alpha^*\\approx0.87$', '$\\alpha=0.3$']
    ax.legend(lines, labels, fontsize=9, loc='center left')
    ax.grid(True, alpha=0.35)

    plt.suptitle('Atenuación en suelo (motivación banda P) y sensibilidad al peso $\\alpha$',
                 fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    save('pv_fig09_alpha_sensitivity.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 10: Vista 3D del sistema completo
# ─────────────────────────────────────────────────────────────────────────────
def fig10_sistema_3d():
    fig = plt.figure(figsize=(13, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Parámetros de la hélice cónica del TX
    rho_min, rho_max = 147.5, 172.5
    z_min, z_max = 80.0, 120.0
    N_turns = 2
    n_pts = 600
    t = np.linspace(0, 2*np.pi*N_turns, n_pts)
    rho_t = rho_min + (rho_max - rho_min) * t / (2*np.pi*N_turns)
    z_t   = z_max  - (z_max - z_min)   * t / (2*np.pi*N_turns)
    PxT = rho_t * np.cos(t)
    PyT = rho_t * np.sin(t)
    PzT = z_t

    # Interfaz z=0
    xg = np.linspace(-200, 200, 30)
    yg = np.linspace(-200, 200, 30)
    XG, YG = np.meshgrid(xg, yg)
    ZG = np.zeros_like(XG)
    ax.plot_surface(XG, YG, ZG, color=C['soil'], alpha=0.12, zorder=1)
    ax.plot_wireframe(XG[::3, ::3], YG[::3, ::3], ZG[::3, ::3],
                      color=C['soil'], lw=0.3, alpha=0.25, zorder=2)

    # Hélice TX
    ax.plot(PxT, PyT, PzT, '-', color=C['tx'], lw=2.0, label='Hélice TX', zorder=5)

    # Grid de targets (4×4×4 = 64 puntos, simplificado 3D)
    n1, n2, theta_B = 1.0, 2.0, np.arctan(2.0)
    tg_x = np.linspace(-40, 40, 4)
    tg_y = np.linspace(-40, 40, 4)
    tg_z = np.linspace(-10, -2, 4)
    for tx_ in tg_x:
        for ty_ in tg_y:
            for tz_ in tg_z:
                ax.plot([tx_], [ty_], [tz_], 'x', ms=4, mew=1.2,
                        color=C['target'], alpha=0.6, zorder=4)
    ax.plot([], [], 'x', ms=6, color=C['target'], label='Grid de targets')

    # Posiciones RX decimadas (simulación: misma hélice, NorthOffset=90°)
    decim = 25
    idx_dec = range(0, n_pts, decim)
    rx_x_list, rx_y_list, rx_z_list = [], [], []
    for k in idx_dec:
        tx_pos = np.array([PxT[k], PyT[k], PzT[k]])
        # RX guiado por Brewster (opuesto + distancia Brewster aproximada)
        tx_dir = tx_pos[:2] / np.linalg.norm(tx_pos[:2])
        d_brew = (PzT[k] + 5.0) / np.tan(theta_B)
        rx_xy = -tx_dir * min(d_brew, 190.0)
        rx_x_list.append(rx_xy[0])
        rx_y_list.append(rx_xy[1])
        rx_z_list.append(PzT[k])

    ax.scatter(rx_x_list, rx_y_list, rx_z_list, s=25, color=C['rx'],
               zorder=6, label='RX óptimo (decimado)')

    # Caminos refractados para una posición TX seleccionada (k=100)
    k_show = 100
    TX_show = np.array([PxT[k_show], PyT[k_show], PzT[k_show]])
    RX_show = np.array([rx_x_list[k_show//decim], rx_y_list[k_show//decim], rx_z_list[k_show//decim]])
    tg_show = np.array([0.0, 0.0, -5.0])  # target central

    # Refracción TX→Target (Newton-Raphson exacto)
    D_TX = np.linalg.norm(TX_show[:2] - tg_show[:2])
    sz_TX, tz_TX = TX_show[2], abs(tg_show[2])
    u_TX, _ = fermat_newton(sz_TX, tz_TX, D_TX, n1, n2)
    dir_TX = (tg_show[:2] - TX_show[:2]) / D_TX
    Q1_show = np.array([TX_show[0] + u_TX*dir_TX[0],
                        TX_show[1] + u_TX*dir_TX[1], 0.0])

    # Refracción Target→RX
    D_RX = np.linalg.norm(tg_show[:2] - RX_show[:2])
    sz_RX, tz_RX = abs(tg_show[2]), RX_show[2]
    u_RX, _ = fermat_newton(sz_RX, tz_RX, D_RX, n2, n1)
    dir_RX = (RX_show[:2] - tg_show[:2]) / D_RX
    Q2_show = np.array([tg_show[0] + u_RX*dir_RX[0],
                        tg_show[1] + u_RX*dir_RX[1], 0.0])

    for seg, col in [([TX_show, Q1_show, tg_show], C['ray']),
                     ([tg_show, Q2_show, RX_show], C['rx'])]:
        xs_ = [p[0] for p in seg]
        ys_ = [p[1] for p in seg]
        zs_ = [p[2] for p in seg]
        ax.plot(xs_, ys_, zs_, '-', color=col, lw=2.5, zorder=8)

    # Puntos
    ax.plot(*TX_show, 's', ms=9, color=C['tx'], zorder=10)
    ax.plot(*RX_show, 's', ms=9, color=C['rx'], zorder=10)
    ax.plot(*tg_show, '*', ms=12, color=C['target'], zorder=10)
    ax.plot(*Q1_show, 'D', ms=7, color=C['orange'], zorder=10)
    ax.plot(*Q2_show, 'D', ms=7, color=C['orange'], zorder=10)

    ax.set_xlabel('$x$ [m]')
    ax.set_ylabel('$y$ [m]')
    ax.set_zlabel('$z$ [m]')
    ax.set_title('Sistema SAR Biestático Helicoidal — Plan de Vuelo Optimizado',
                 fontsize=12, fontweight='bold', pad=15)
    ax.legend(fontsize=9, loc='upper left')
    ax.set_xlim(-200, 200)
    ax.set_ylim(-200, 200)
    ax.set_zlim(-15, 130)
    ax.view_init(elev=22, azim=45)
    save('pv_fig10_sistema_3d.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print('Generando figuras para explicacion_plan_de_vuelo.tex ...')
    print()

    print('[1/8] Fig 01 — Coeficientes de Fresnel TM ...')
    fig01_fresnel_tm()

    print('[2/8] Fig 03 — Snell 1D y Newton-Raphson ...')
    fig03_snell_1d()

    print('[3/8] Fig 04 — Ecuación de radar biestático ...')
    fig04_bistatic_radar()

    print('[4/8] Fig 05 — Función de costo J(xR,yR) ...')
    fig05_cost_function()

    print('[5/8] Fig 06 — Candidatos multi-inicio ...')
    fig06_multistart()

    print('[6/8] Fig 07 — Frontera de Pareto ...')
    fig07_pareto()

    print('[7/8] Fig 08 — Geometría de Brewster ...')
    fig08_brewster_geometry()

    print('[8/8] Fig 09+10 — Atenuación/sensibilidad alpha y 3D sistema ...')
    fig09_alpha_sensitivity()
    fig10_sistema_3d()

    print()
    print(f'Todas las figuras guardadas en: {FIG_DIR}')
    print('Para compilar el documento:')
    print('  cd hypotheses')
    print('  pdflatex explicacion_plan_de_vuelo.tex')
    print('  pdflatex explicacion_plan_de_vuelo.tex  # segunda pasada para referencias')
