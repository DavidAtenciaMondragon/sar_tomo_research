"""
generate_figures.py
Genera todas las figuras para explicacion_resolucion.tex
Guarda en: hypotheses/figures/
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch, Arc, FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.special import j0
import os

# ── directorios ──────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(SCRIPT_DIR, '..', 'hypotheses', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# ── estilo global ─────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'DejaVu Serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'lines.linewidth': 2.0,
    'grid.alpha': 0.35,
    'figure.dpi': 150,
    'savefig.dpi': 200,
    'savefig.bbox': 'tight',
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
}

def save(name):
    path = os.path.join(FIG_DIR, name)
    plt.savefig(path, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    print(f'  Saved: {name}')

# ─────────────────────────────────────────────────────────────────────────────
# FIG 1: Geometría básica de un radar
# ─────────────────────────────────────────────────────────────────────────────
def fig01_radar_basics():
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_xlim(-1, 9); ax.set_ylim(-0.5, 6)
    ax.axis('off')

    # suelo
    ax.axhline(0, color='#8c510a', lw=2.5)
    ax.fill_between([-1, 9], [-0.5, -0.5], [0, 0], color='#d5a97a', alpha=0.4)
    ax.text(8.5, -0.25, 'Suelo\n($z=0$)', ha='right', va='center', fontsize=9, color='#5a3a1a')

    # sensor
    H = 4.5
    sx = 1.5
    ax.plot(sx, H, 's', ms=14, color=C['sky'], zorder=5)
    ax.text(sx - 0.25, H + 0.2, 'Sensor\n(radar)', ha='center', va='bottom', fontsize=10, color=C['sky'])

    # target
    tx = 6.5
    ax.plot(tx, 0, '*', ms=14, color=C['target'], zorder=5)
    ax.text(tx + 0.2, 0.15, 'Objetivo\n(target)', ha='left', va='bottom', fontsize=10, color=C['target'])

    # slant range
    ax.annotate('', xy=(tx, 0), xytext=(sx, H),
                arrowprops=dict(arrowstyle='->', color=C['ray'], lw=2.5))
    mid_x = (sx + tx) / 2 + 0.1
    mid_y = H / 2
    ax.text(mid_x + 0.3, mid_y, r'$R$ (slant range)', rotation=-45, fontsize=11,
            color=C['ray'], ha='left')

    # altura H
    ax.annotate('', xy=(sx, 0), xytext=(sx, H),
                arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
    ax.text(sx - 0.15, H / 2, '$H$', ha='right', va='center', fontsize=12)

    # ground range
    ax.annotate('', xy=(tx, -0.35), xytext=(sx, -0.35),
                arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
    ax.text((sx + tx) / 2, -0.45, r'$r_g$ (ground range)', ha='center', va='top', fontsize=10)

    # ángulo de incidencia ψ
    ang_rad = np.arctan2(tx - sx, H)
    arc_ang = Arc((sx, H), 1.2, 1.2, angle=0, theta1=270 - np.degrees(ang_rad), theta2=270,
                  color='black', lw=1.5)
    ax.add_patch(arc_ang)
    ax.text(sx + 0.7, H - 0.75, r'$\psi$', fontsize=13)

    # pulso chirp (esquema)
    bx = np.linspace(7, 8.5, 200)
    chirp_t = np.linspace(0, 1, 200)
    chirp_y = 0.4 * np.sin(2 * np.pi * (chirp_t + 2 * chirp_t**2)) * np.exp(-3 * (chirp_t - 0.5)**2)
    ax.plot(bx, chirp_y + 3.5, color=C['green'], lw=1.5)
    ax.text(7.75, 3.1, 'Pulso LFM\n(bandwidth $B$)', ha='center', va='top', fontsize=9, color=C['green'])
    ax.annotate('', xy=(8.5, 2.8), xytext=(7.0, 2.8),
                arrowprops=dict(arrowstyle='<->', color=C['green'], lw=1.5))
    ax.text(7.75, 2.65, r'$T_p$', ha='center', fontsize=10, color=C['green'])

    # resolución en rango
    ax.text(4.5, 5.4, r'Resolución en rango: $\delta_r = \dfrac{c}{2B}$',
            ha='center', fontsize=12,
            bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', ec='#aaa', alpha=0.9))

    ax.set_title('Geometría básica de un radar pulsado', fontsize=13, pad=8)
    save('fig01_radar_basics.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 2: Apertura real vs. sintética
# ─────────────────────────────────────────────────────────────────────────────
def fig02_aperture_real_vs_synth():
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    for ax in axes:
        ax.set_xlim(-5, 5); ax.set_ylim(-0.5, 5.5)
        ax.axis('off')
        ax.axhline(0, color='#8c510a', lw=2)
        ax.fill_between([-5, 5], [-0.5, -0.5], [0, 0], color='#d5a97a', alpha=0.35)

    # ── apertura real ──────────────────────────────────────────────────────
    ax = axes[0]
    ax.set_title('(a) Apertura Real (RAR)', fontsize=12, pad=6)
    H = 4.0
    # antena física
    ax.plot(0, H, 'sb', ms=16, zorder=5)
    ant_w = 1.2
    ax.fill_between([-ant_w/2, ant_w/2], [H-0.15, H-0.15], [H+0.15, H+0.15],
                    color=C['sky'], alpha=0.6)
    ax.text(0, H + 0.35, r'$d_a$ (antena física)', ha='center', fontsize=10, color=C['sky'])

    # haz
    half_bw = np.degrees(np.arctan(0.5 * 0.05 / 0.3))  # λ/da
    ang = np.radians(8)
    for s in [-1, 1]:
        ax.plot([0, s * H * np.tan(ang)], [H, 0], '--', color=C['ray'], alpha=0.7)
    ax.fill_between([- H * np.tan(ang), H * np.tan(ang)], [0, 0], [0, 0],
                    alpha=0)
    ax.fill([0, -H * np.tan(ang), H * np.tan(ang)], [H, 0, 0],
            color=C['ray'], alpha=0.15)

    target_x = 0
    ax.plot(target_x, 0, '*', ms=14, color=C['target'], zorder=5)
    bw = 2 * H * np.tan(ang)
    ax.annotate('', xy=(bw/2, -0.35), xytext=(-bw/2, -0.35),
                arrowprops=dict(arrowstyle='<->', color=C['target'], lw=1.5))
    ax.text(0, -0.45, r'$\delta_a^{RAR} = \frac{\lambda}{d_a}\,R \propto R$',
            ha='center', va='top', fontsize=10, color=C['target'])
    ax.text(0, 1.5, r'$\Theta_{3dB} = \lambda/d_a$', ha='center', fontsize=10, color=C['ray'])

    # ── apertura sintética ─────────────────────────────────────────────────
    ax = axes[1]
    ax.set_title('(b) Apertura Sintética (SAR)', fontsize=12, pad=6)
    H = 4.0
    # trayectoria
    xs = np.linspace(-3, 3, 7)
    ax.plot(xs, [H] * len(xs), 'D', ms=7, color=C['sky'], zorder=4)
    ax.plot([-3.5, 3.5], [H, H], '-', color=C['sky'], lw=1.5, alpha=0.5)
    # apertura sintética
    ax.annotate('', xy=(3, H + 0.4), xytext=(-3, H + 0.4),
                arrowprops=dict(arrowstyle='<->', color=C['sky'], lw=2))
    ax.text(0, H + 0.7, r'$L_{SA}$ (apertura sintética)', ha='center', fontsize=10, color=C['sky'])

    target_x = 0
    ax.plot(target_x, 0, '*', ms=14, color=C['target'], zorder=5)
    # rayos desde extremos
    for xs_e in [-3, 3]:
        ax.plot([xs_e, target_x], [H, 0], '--', color=C['ray'], alpha=0.5, lw=1.2)
    ax.fill([-3, target_x, 3], [H, 0, H], color=C['ray'], alpha=0.1)

    bw_sar = 0.5
    ax.annotate('', xy=(bw_sar/2, -0.35), xytext=(-bw_sar/2, -0.35),
                arrowprops=dict(arrowstyle='<->', color=C['target'], lw=1.5))
    ax.text(0, -0.45, r'$\delta_a^{SAR} = \frac{d_a}{2}$ (independiente de $R$)',
            ha='center', va='top', fontsize=10, color=C['target'])
    ax.text(0, 2.5, 'Coherente en\ntodos los pulsos', ha='center', fontsize=9, color=C['gray'],
            style='italic')

    fig.suptitle('Resolución en Azimut: Apertura Real vs. Apertura Sintética', fontsize=13, y=1.01)
    plt.tight_layout()
    save('fig02_real_vs_synth.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 3: SAR circular — trayectoria + k-espacio + PSF
# ─────────────────────────────────────────────────────────────────────────────
def fig03_circular_sar():
    fig = plt.figure(figsize=(13, 4.5))

    # ── (a) trayectoria 3D ──────────────────────────────────────────────────
    ax1 = fig.add_subplot(131, projection='3d')
    rho = 1.5; H = 2.0; alpha = np.linspace(0, 2 * np.pi, 200)
    ax1.plot(rho * np.cos(alpha), rho * np.sin(alpha), H, color=C['sky'], lw=2.5,
             label='Trayectoria')
    ax1.scatter([0], [0], [0], c=C['target'], s=80, zorder=5, marker='*',
                label='Objetivo')
    # líneas de rango desde algunos puntos
    for a in np.linspace(0, 2 * np.pi, 7, endpoint=False):
        ax1.plot([rho * np.cos(a), 0], [rho * np.sin(a), 0], [H, 0],
                 '--', color=C['ray'], alpha=0.4, lw=1.0)
    # plano interfaz
    u = np.linspace(-2, 2, 2)
    v = np.linspace(-2, 2, 2)
    U, V = np.meshgrid(u, v)
    ax1.plot_surface(U, V, np.zeros_like(U), alpha=0.12, color='#8c510a')
    ax1.set_xlabel('x'); ax1.set_ylabel('y'); ax1.set_zlabel('z')
    ax1.set_title('(a) Trayectoria\nCircular', fontsize=11)
    ax1.legend(fontsize=8, loc='upper left')

    # ángulos
    psi0 = np.arctan2(rho, H)
    ax1.quiver(0, 0, 0, 0, 0, H, length=H * 0.2, arrow_length_ratio=0.3, color='black')
    ax1.text(0.1, 0, H * 0.85, r'$\psi_0$', fontsize=11)

    # ── (b) k-espacio ──────────────────────────────────────────────────────
    ax2 = fig.add_subplot(132)
    # anillo del k-espacio
    f0 = 1.0; B = 0.3; lam = 1.0 / f0
    kc = 4 * np.pi * np.sin(psi0) / lam
    k_in = 4 * np.pi * (f0 - B / 2) * np.sin(psi0)
    k_out = 4 * np.pi * (f0 + B / 2) * np.sin(psi0)
    theta = np.linspace(0, 2 * np.pi, 300)
    for r, ls in [(k_in, '--'), (k_out, '-')]:
        ax2.plot(r * np.cos(theta), r * np.sin(theta), ls, color=C['green'], lw=2)
    ax2.fill_between(k_out * np.cos(theta), k_out * np.sin(theta),
                     k_in * np.cos(theta), alpha=0.2, color=C['green'])
    ax2.annotate('', xy=(k_out, 0), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
    ax2.text(k_out * 0.5, 0.15, r'$R_c = \frac{4\pi\sin\psi_0}{\lambda}$', fontsize=10)
    ax2.annotate('', xy=(k_in, -0.2), xytext=(k_out, -0.2),
                arrowprops=dict(arrowstyle='<->', color='darkblue', lw=1.5))
    ax2.text((k_in + k_out) / 2, -0.35, r'$\Delta k_r = \frac{4\pi B}{c}$',
             ha='center', fontsize=9, color='darkblue')
    ax2.set_xlabel(r'$k_x$ [m$^{-1}$]'); ax2.set_ylabel(r'$k_y$ [m$^{-1}$]')
    ax2.set_title('(b) Cobertura\nen k-espacio', fontsize=11)
    ax2.set_aspect('equal'); ax2.grid(True)

    # ── (c) PSF (perfil J0) ─────────────────────────────────────────────────
    ax3 = fig.add_subplot(133)
    r = np.linspace(-3, 3, 500)
    psf = j0(kc * r)
    ax3.plot(r, psf, color=C['purple'], lw=2.5)
    ax3.axhline(0, color='black', lw=0.8)
    ax3.axhline(1 / np.sqrt(2), color='gray', ls='--', lw=1.5, label=r'$-3\,\mathrm{dB}$')
    # marcar los puntos −3dB
    u_3db = 1.20  # J0(u)=1/sqrt(2) en u≈1.20
    r_3db = u_3db / kc
    ax3.axvline(-r_3db, color=C['orange'], ls=':', lw=1.5)
    ax3.axvline(r_3db, color=C['orange'], ls=':', lw=1.5)
    ax3.annotate('', xy=(r_3db, 0.5), xytext=(-r_3db, 0.5),
                arrowprops=dict(arrowstyle='<->', color=C['orange'], lw=2))
    ax3.text(0, 0.55, r'$\delta_{xy}$', ha='center', fontsize=13, color=C['orange'])
    ax3.text(0, -0.25, r'$\delta_{xy} = \dfrac{1.20\lambda}{4\pi\sin\psi_0}$',
             ha='center', fontsize=10,
             bbox=dict(boxstyle='round', fc='lightyellow', ec='#ccc', alpha=0.9))
    ax3.set_xlabel(r'Desplazamiento $\delta$ [m]')
    ax3.set_ylabel(r'PSF$_{xy}$ (normalizada)')
    ax3.set_title('(c) PSF Horizontal\n($J_0$ de Bessel)', fontsize=11)
    ax3.legend(fontsize=9, loc='upper right')
    ax3.grid(True); ax3.set_ylim(-0.45, 1.15)

    plt.tight_layout()
    save('fig03_circular_sar.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 4: Hélice cónica con todos los parámetros geométricos
# ─────────────────────────────────────────────────────────────────────────────
def fig04_helical_params():
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # parámetros
    rho_top = 1.0; rho_base = 1.6; z_top = 2.5; z_base = 0.8
    N = 3
    t = np.linspace(0, 1, 600)
    rho_t = rho_top + (rho_base - rho_top) * t
    z_t = z_top + (z_base - z_top) * t
    alpha_t = 2 * np.pi * N * t

    x_t = rho_t * np.cos(alpha_t)
    y_t = rho_t * np.sin(alpha_t)

    # trayectoria TX (azul)
    ax.plot(x_t, y_t, z_t, color=C['sky'], lw=2.5, label='TX ($\phi_{TX}=0°$)')
    # trayectoria RX (rojo, offset 90°)
    offset = np.pi / 2
    ax.plot(rho_t * np.cos(alpha_t + offset), rho_t * np.sin(alpha_t + offset), z_t,
            color=C['target'], lw=2, ls='--', label='RX ($\Delta\phi=90°$)')

    # target
    ax.scatter([0], [0], [0], c=C['target'], s=120, marker='*', zorder=5)
    ax.text(0.15, 0.15, 0, 'Target\n$(0,0,z_P)$', fontsize=10)

    # plano interfaz
    u = np.linspace(-2, 2, 2)
    v = np.linspace(-2, 2, 2)
    U, V = np.meshgrid(u, v)
    ax.plot_surface(U, V, np.zeros_like(U), alpha=0.1, color='#8c510a')

    # suelo
    ax.plot_surface(U, V, -0.15 * np.ones_like(U), alpha=0.25, color='#d5a97a')

    # parámetros radiales
    i_top = 0; i_base = -1
    # radio en cima
    ax.plot([0, rho_top], [0, 0], [z_top, z_top], 'k--', lw=1.5)
    ax.text(rho_top / 2, 0.1, z_top + 0.1, r'$\rho_{top}$', fontsize=11, color='black')
    # radio en base
    ax.plot([0, rho_base], [0, 0], [z_base, z_base], 'b--', lw=1.5)
    ax.text(rho_base / 2, 0.1, z_base - 0.25, r'$\rho_{base}$', fontsize=11, color='blue')

    # alturas
    ax.plot([1.8, 1.8], [0, 0], [z_base, z_top], 'g-', lw=2)
    ax.text(2.0, 0, (z_top + z_base) / 2, r'$\Delta z$', fontsize=11, color='green')
    # Δρ
    ax.plot([0, rho_base - rho_top], [1.5, 1.5], [z_base, z_base], 'm-', lw=2)
    ax.text((rho_base - rho_top) / 2, 1.6, z_base - 0.3, r'$\Delta\rho$', fontsize=11, color='m')

    # B_helix (vector)
    ax.quiver(rho_top, 0, z_top, rho_base - rho_top, 0, z_base - z_top,
              color='darkgreen', arrow_length_ratio=0.15, lw=2)
    ax.text((rho_top + rho_base) / 2 + 0.1, 0.2, (z_top + z_base) / 2 + 0.15,
            r'$B_{helix}$', fontsize=12, color='darkgreen')

    # eje óptico (LOS medio)
    rho0 = (rho_top + rho_base) / 2; z0 = (z_top + z_base) / 2
    ax.quiver(0, 0, 0, rho0, 0, z0, color='darkorange', arrow_length_ratio=0.1, lw=1.5)
    ax.text(rho0 / 2 + 0.1, 0.2, z0 / 2 - 0.2, r'$R_0$, $\psi_0$', fontsize=11, color='darkorange')

    # ángulo β
    ax.text(rho_top + 0.1, 0.1, z_top - 0.3, r'$\beta$', fontsize=13, color='darkgreen')

    # B_perp (perpendicular a LOS)
    psi0 = np.arctan2(rho0, z0)
    beta = np.arctan2(z_top - z_base, rho_base - rho_top)
    Bperp_len = np.sqrt((z_top - z_base)**2 + (rho_base - rho_top)**2) * abs(np.cos(beta - psi0))
    ax.text(-1.5, -1.5, 1.2, r'$B_\perp = B_{helix}|\cos(\beta-\psi_0)|$',
            fontsize=10, color='darkblue',
            bbox=dict(boxstyle='round', fc='lightyellow', ec='#aaa', alpha=0.85))

    ax.set_xlabel('x [m]'); ax.set_ylabel('y [m]'); ax.set_zlabel('z [m]')
    ax.set_title('Parámetros de la Hélice Cónica Biestática', fontsize=13, pad=10)
    ax.legend(loc='upper left', fontsize=9)
    ax.view_init(elev=20, azim=-55)

    save('fig04_helical_params.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 5: Ley de Snell — geometría en dos medios
# ─────────────────────────────────────────────────────────────────────────────
def fig05_snell_geometry():
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.set_xlim(-4, 6); ax.set_ylim(-4, 5)
    ax.axis('off')

    # medios
    ax.fill_between([-4, 6], [0, 0], [5, 5], color='#deebf7', alpha=0.7,
                    label='Aire ($n_1=1$)')
    ax.fill_between([-4, 6], [-4, -4], [0, 0], color='#d5c4a1', alpha=0.5,
                    label='Suelo ($n_2=\sqrt{\\varepsilon_r}$)')
    ax.axhline(0, color='#5a3a1a', lw=2.5)

    # etiquetas de medio
    ax.text(5.5, 3.5, 'Medio 1\n(Aire)\n$n_1=1$\n$v_1=c$',
            ha='right', va='top', fontsize=10, color=C['sky'],
            bbox=dict(boxstyle='round', fc='white', ec=C['sky'], alpha=0.8))
    ax.text(5.5, -3.5, 'Medio 2\n(Suelo)\n$n_2=\sqrt{\\varepsilon_r}$\n$v_2=c/n_2$',
            ha='right', va='bottom', fontsize=10, color=C['soil'],
            bbox=dict(boxstyle='round', fc='white', ec=C['soil'], alpha=0.8))

    # sensor TX
    sensor_x = -1.5; sensor_z = 3.5
    ax.plot(sensor_x, sensor_z, 's', ms=14, color=C['sky'], zorder=5)
    ax.text(sensor_x - 0.2, sensor_z + 0.2, 'TX', ha='right', fontsize=11, color=C['sky'])

    # sensor RX
    rx_x = 4.5; rx_z = 2.8
    ax.plot(rx_x, rx_z, 's', ms=14, color=C['green'], zorder=5)
    ax.text(rx_x + 0.2, rx_z + 0.2, 'RX', ha='left', fontsize=11, color=C['green'])

    # puntos de refracción
    qTx = 0.5; qRx = 2.8
    ax.plot(qTx, 0, 'o', ms=10, color='black', zorder=5)
    ax.text(qTx, 0.2, r'$\mathbf{Q}_{TX}^*$', ha='center', fontsize=10)
    ax.plot(qRx, 0, 'o', ms=10, color='black', zorder=5)
    ax.text(qRx, 0.2, r'$\mathbf{Q}_{RX}^*$', ha='center', fontsize=10)

    # target
    target_x = 1.5; target_z = -2.8
    ax.plot(target_x, target_z, '*', ms=18, color=C['target'], zorder=5)
    ax.text(target_x + 0.25, target_z - 0.1, r'$\mathbf{P}=(x_P,y_P,z_P)$',
            ha='left', fontsize=10, color=C['target'])

    # rayos TX → Q → P
    ax.annotate('', xy=(qTx, 0), xytext=(sensor_x, sensor_z),
                arrowprops=dict(arrowstyle='->', color=C['sky'], lw=2.5))
    ax.annotate('', xy=(target_x, target_z), xytext=(qTx, 0),
                arrowprops=dict(arrowstyle='->', color=C['sky'], lw=2.5))

    # rayos P → Q → RX
    ax.annotate('', xy=(qRx, 0), xytext=(target_x, target_z),
                arrowprops=dict(arrowstyle='->', color=C['green'], lw=2.5, ls='--'))
    ax.annotate('', xy=(rx_x, rx_z), xytext=(qRx, 0),
                arrowprops=dict(arrowstyle='->', color=C['green'], lw=2.5, ls='--'))

    # normal a la interfaz (en Q_TX)
    ax.plot([qTx, qTx], [-1.0, 1.5], 'k:', lw=1.5)
    ax.text(qTx + 0.1, 1.2, r'$\hat{n}$', fontsize=12)

    # ángulo de incidencia θ_i (TX)
    dx_inc = qTx - sensor_x; dz_inc = -(0 - sensor_z)
    ang_i = np.degrees(np.arctan2(abs(dx_inc), dz_inc))
    arc1 = Arc((qTx, 0), 1.4, 1.4, angle=0, theta1=90, theta2=90 + ang_i,
               color=C['sky'], lw=1.5)
    ax.add_patch(arc1)
    ax.text(qTx - 0.9, 0.7, r'$\theta_i$', fontsize=13, color=C['sky'])

    # ángulo de transmisión θ_t (TX)
    dx_tr = target_x - qTx; dz_tr = -(target_z - 0)
    ang_t = np.degrees(np.arctan2(abs(dx_tr), dz_tr))
    arc2 = Arc((qTx, 0), 1.2, 1.2, angle=0, theta1=270, theta2=270 + ang_t,
               color=C['orange'], lw=1.5)
    ax.add_patch(arc2)
    ax.text(qTx + 0.5, -0.8, r'$\theta_t$', fontsize=13, color=C['orange'])

    # distancias
    mid_air_x = (sensor_x + qTx) / 2; mid_air_z = (sensor_z + 0) / 2
    ax.text(mid_air_x - 0.5, mid_air_z, r'$d_1^{TX}$', fontsize=11, color=C['sky'], rotation=50)
    mid_soil_x = (qTx + target_x) / 2; mid_soil_z = (0 + target_z) / 2
    ax.text(mid_soil_x - 0.7, mid_soil_z, r'$d_2^{TX}$', fontsize=11, color=C['sky'], rotation=55)
    mid_air2_x = (qRx + rx_x) / 2; mid_air2_z = (0 + rx_z) / 2
    ax.text(mid_air2_x + 0.1, mid_air2_z, r'$d_1^{RX}$', fontsize=11, color=C['green'], rotation=-30)
    mid_soil2_x = (target_x + qRx) / 2; mid_soil2_z = (target_z + 0) / 2
    ax.text(mid_soil2_x + 0.1, mid_soil2_z, r'$d_2^{RX}$', fontsize=11, color=C['green'], rotation=-30)

    # ley de Snell
    ax.text(0, -3.5,
            r'Ley de Snell: $n_1\sin\theta_i = n_2\sin\theta_t$' + '\n' +
            r'Camino óptico: $R_{OP} = d_1^{TX} + n_2\,d_2^{TX} + n_2\,d_2^{RX} + d_1^{RX}$',
            ha='left', fontsize=10.5,
            bbox=dict(boxstyle='round', fc='lightyellow', ec='#999', alpha=0.9))

    ax.set_title('Geometría de Propagación Biestática en Dos Medios\n'
                 '(Principio de Fermat determina los puntos de refracción)',
                 fontsize=12, pad=8)
    save('fig05_snell_geometry.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 6: k-espacio vertical — dos contribuciones a ΔkZ
# ─────────────────────────────────────────────────────────────────────────────
def fig06_kspace_vertical():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # parámetros ilustrativos
    f0 = 10e9; B = 0.5e9; c = 3e8; n2 = 2.0; lam = c / f0
    psi0 = np.radians(58)
    psi_top = np.radians(51); psi_bot = np.radians(65)
    cos_t_top = np.sqrt(1 - np.sin(psi_top)**2 / n2**2)
    cos_t_bot = np.sqrt(1 - np.sin(psi_bot)**2 / n2**2)
    cos_t_mean = (cos_t_top + cos_t_bot) / 2

    # ── (a) variación de k_z con frecuencia y tiempo ─────────────────────
    ax = axes[0]
    f_arr = np.array([f0 - B / 2, f0, f0 + B / 2])
    # k_z para dos geometrías (cima y base)
    for cos_t, label, col, ls in [(cos_t_top, 'Cima ($\\psi_{top}$)', C['sky'], '-'),
                                    (cos_t_bot, 'Base ($\\psi_{bot}$)', C['soil'], '--')]:
        kz = 4 * np.pi * n2 * f_arr * cos_t / c
        ax.plot(f_arr / 1e9, kz, ls=ls, color=col, lw=2.5, label=label)

    # sombrear la región cubierta
    kz_max = 4 * np.pi * n2 * (f0 + B / 2) * cos_t_top / c
    kz_min = 4 * np.pi * n2 * (f0 - B / 2) * cos_t_bot / c
    ax.fill_between([f0 / 1e9 - B / 2e9, f0 / 1e9 + B / 2e9],
                    [kz_min, kz_min], [kz_max, kz_max], alpha=0.2, color=C['purple'])

    ax.annotate('', xy=(f0 / 1e9 + B / 2e9 + 0.05, kz_min),
                xytext=(f0 / 1e9 + B / 2e9 + 0.05, kz_max),
                arrowprops=dict(arrowstyle='<->', color=C['purple'], lw=2))
    ax.text(f0 / 1e9 + B / 2e9 + 0.12, (kz_min + kz_max) / 2,
            r'$\Delta k_z$', fontsize=13, color=C['purple'], va='center')

    ax.set_xlabel('Frecuencia [GHz]')
    ax.set_ylabel(r'$|k_z|$ [rad/m]')
    ax.set_title('(a) Variación de $k_z$ con\nfrecuencia y geometría', fontsize=11)
    ax.legend(fontsize=9); ax.grid(True)

    # ── (b) descomposición de W_z ─────────────────────────────────────────
    ax = axes[1]
    # calcular los dos términos
    rho0 = 160; z0 = 100; R0 = np.sqrt(rho0**2 + z0**2)
    B_helix = 47.17; beta = np.radians(58)
    Bperp = B_helix * abs(np.cos(beta - psi0))
    cos_t0 = np.sqrt(1 - np.sin(psi0)**2 / n2**2)
    B_hz = 50e6

    term_bw = n2 * B_hz * cos_t0
    term_tomo = c * Bperp * np.sin(psi0) * np.cos(psi0) / (lam * R0 * n2 * cos_t0)
    Wz = term_bw + term_tomo

    bars = ax.bar(['Ancho de\nbanda\n$n_2 B\cos\\theta_{t,0}$',
                   'Apertura\ntomográfica\n$\\frac{cB_\\perp\\sin\\psi_0\\cos\\psi_0}{\\lambda_0 R_0 n_2\\cos\\theta_{t,0}}$',
                   'Total\n$W_z$'],
                  [term_bw / 1e6, term_tomo / 1e6, Wz / 1e6],
                  color=[C['sky'], C['green'], C['purple']],
                  edgecolor='white', linewidth=1.5, width=0.5)

    for bar, val in zip(bars, [term_bw, term_tomo, Wz]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 5,
                f'{val/1e6:.0f} MHz', ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.set_ylabel(r'$W_z$ [MHz]')
    ax.set_title('(b) Descomposición de $W_z$\n($f_0=10$ GHz, $n_2=2$, $B=50$ MHz)', fontsize=11)
    ax.grid(axis='y', alpha=0.4)

    # anotación de δz
    dz = c / (2 * Wz)
    ax.text(1, Wz / 1e6 * 0.5,
            f'$\\delta_z = c/(2W_z)$\n$= {dz*100:.1f}$ cm',
            ha='center', fontsize=11,
            bbox=dict(boxstyle='round', fc='lightyellow', ec='#aaa', alpha=0.9))

    plt.tight_layout()
    save('fig06_kspace_vertical.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 7: Efecto del NorthOffset en el k-espacio horizontal
# ─────────────────────────────────────────────────────────────────────────────
def fig07_northoffset_kspace():
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    offsets = [0, 90, 180]
    titles = ['(a) NorthOffset = 0°\n(monoestático equivalente)',
              '(b) NorthOffset = 90°',
              '(c) NorthOffset = 180°\n(cancelación total)']
    colors_off = [C['sky'], C['green'], C['target']]

    psi0 = np.radians(58); n2 = 2.0; lam = 0.03
    f0 = 10e9; c = 3e8; B = 0.5e9
    sin_t0 = np.sin(psi0) / n2
    A = 2 * np.pi * f0 * n2 * sin_t0 / c  # radio base para Δφ=0

    theta = np.linspace(0, 2 * np.pi, 300)

    for ax, dphi_deg, title, col in zip(axes, offsets, titles, colors_off):
        dphi = np.radians(dphi_deg)
        Rc = A * 2 * abs(np.cos(dphi / 2))  # radio del círculo

        if Rc > 0:
            k_in = A * 2 * abs(np.cos(dphi / 2)) * (1 - B / (2 * f0))
            k_out = A * 2 * abs(np.cos(dphi / 2)) * (1 + B / (2 * f0))
            ax.fill_between(k_out * np.cos(theta), -k_out * np.sin(theta),
                            k_in * np.cos(theta) * np.sign(np.cos(theta)),
                            alpha=0.15, color=col)
            ax.plot(Rc * np.cos(theta), Rc * np.sin(theta), color=col, lw=2.5)
            # radio
            ax.annotate('', xy=(Rc, 0), xytext=(0, 0),
                        arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
            ax.text(Rc * 0.5, 0.2 * Rc, r'$R_c$', fontsize=12)
            # resolución
            delta_xy = 1.20 / Rc * 1e3  # mm
            ax.text(0, -Rc * 1.3,
                    rf'$\delta_{{xy}} \approx {delta_xy:.1f}$ mm',
                    ha='center', fontsize=10, color=col,
                    bbox=dict(boxstyle='round', fc='white', ec=col, alpha=0.8))
        else:
            ax.text(0, 0, 'Sin cobertura\nhorizontal\n($k_x=k_y=0$)',
                    ha='center', va='center', fontsize=11, color=col,
                    bbox=dict(boxstyle='round', fc='lightyellow', ec=col, alpha=0.9))
            ax.text(0, -0.25,
                    r'$\delta_{xy}\to\infty$',
                    ha='center', fontsize=12, color=col)

        ax.axhline(0, color='gray', lw=0.8); ax.axvline(0, color='gray', lw=0.8)
        ax.set_aspect('equal')
        ax.set_xlabel(r'$k_x$ [rad/m]'); ax.set_ylabel(r'$k_y$ [rad/m]')
        ax.set_title(title, fontsize=10)
        ax.grid(True)
        lim = max(A * 2.2, 0.5)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)

    fig.suptitle(r'Cobertura del k-espacio horizontal según NorthOffset $\Delta\phi$' + '\n' +
                 r'($R_c(\Delta\phi) = \frac{4\pi\sin\psi_0}{\lambda}|\cos(\Delta\phi/2)|$)',
                 fontsize=12, y=1.02)
    plt.tight_layout()
    save('fig07_northoffset_kspace.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 8: Sistema completo 3D — hélice cónica biestática + interfaz + target
# ─────────────────────────────────────────────────────────────────────────────
def fig08_full_system():
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # hélice TX
    N = 2; rho_top = 1.0; rho_base = 1.5; z_top = 2.2; z_base = 0.7
    t = np.linspace(0, 1, 400)
    rho_t = rho_top + (rho_base - rho_top) * t
    z_t = z_top + (z_base - z_top) * t
    alpha_t = 2 * np.pi * N * t
    xTX = rho_t * np.cos(alpha_t); yTX = rho_t * np.sin(alpha_t)
    xRX = rho_t * np.cos(alpha_t + np.pi/2); yRX = rho_t * np.sin(alpha_t + np.pi/2)

    ax.plot(xTX, yTX, z_t, color=C['sky'], lw=2.5, label='TX (NorthOffset=0°)')
    ax.plot(xRX, yRX, z_t, color=C['target'], lw=2, ls='--', label='RX (NorthOffset=90°)')

    # puntos TX y RX actuales (instante medio)
    i_mid = len(t) // 2
    ax.scatter([xTX[i_mid]], [yTX[i_mid]], [z_t[i_mid]], c=C['sky'], s=100, zorder=6, marker='s')
    ax.scatter([xRX[i_mid]], [yRX[i_mid]], [z_t[i_mid]], c=C['target'], s=100, zorder=6, marker='s')

    # interfaz (suelo superior)
    u = np.linspace(-2, 2, 3)
    v = np.linspace(-2, 2, 3)
    U, V = np.meshgrid(u, v)
    ax.plot_surface(U, V, np.zeros_like(U), alpha=0.12, color='#8c510a')
    ax.text(1.8, 1.8, 0.05, '$z=0$\n(Interfaz)', fontsize=9, color='#5a3a1a')

    # suelo
    ax.plot_surface(U, V, -0.25 * np.ones_like(U), alpha=0.3, color='#d5a97a')
    ax.text(-1.8, -1.8, -0.2, 'Suelo\n$n_2=\sqrt{\\varepsilon_r}$', fontsize=9, color='#5a3a1a')

    # aire label
    ax.text(1.5, 1.5, 1.8, 'Aire\n$n_1=1$', fontsize=9, color=C['sky'])

    # target enterrado
    ax.scatter([0], [0], [-1.0], c=C['target'], s=160, marker='*', zorder=6)
    ax.text(0.2, 0.1, -1.1, r'$\mathbf{P}=(0,0,z_P)$', fontsize=10, color=C['target'])

    # rayos de propagación (algunos pulsos)
    for idx in [50, 200, 350]:
        # rayo TX → target (aproximado, a través de interfaz)
        qx_TX = xTX[idx] * abs(z_t[idx]) / (z_t[idx] + 0)  # refracción simplificada
        qy_TX = yTX[idx] * abs(z_t[idx]) / (z_t[idx] + 0)
        qx_TX = 0.3 * xTX[idx]
        qy_TX = 0.3 * yTX[idx]
        ax.plot([xTX[idx], qx_TX, 0], [yTX[idx], qy_TX, 0], [z_t[idx], 0, -1.0],
                '--', color=C['ray'], alpha=0.35, lw=1.0)

    # dimensiones
    ax.quiver(0, 0, -1.0, 0, 0, 1.0, color='gray', arrow_length_ratio=0.15, lw=1.5)
    ax.text(0.1, 0.1, -0.5, r'$|z_P|$', fontsize=11, color='gray')

    ax.set_xlabel('x [m]'); ax.set_ylabel('y [m]'); ax.set_zlabel('z [m]')
    ax.set_title('Sistema Completo: SAR Biestático Helicoidal\n'
                 'con Dos Medios e Interfaz Plana', fontsize=12, pad=10)
    ax.legend(loc='upper left', fontsize=9)
    ax.view_init(elev=22, azim=-50)

    save('fig08_full_system.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 9: Comparación de PSFs (pancake vs elipsoide)
# ─────────────────────────────────────────────────────────────────────────────
def fig09_psf_comparison():
    fig = plt.figure(figsize=(12, 5))

    for col_idx, (dphi_deg, title, dxy_m, dz_m) in enumerate([
        (180, 'NorthOffset=180°\n(sin resolución xy)', 10.0, 0.21),
        (90,  'NorthOffset=90°\n($\\delta_{xy}\\approx 9$ mm)', 0.009, 0.21),
        (0,   'NorthOffset=0°\n(monoestático)', 0.0045, 0.21),
    ]):
        ax = fig.add_subplot(1, 3, col_idx + 1, projection='3d')

        # isosuperficie simplificada como elipsoide
        u = np.linspace(0, np.pi, 30)
        v = np.linspace(0, 2 * np.pi, 30)
        rx = dxy_m * 10 * np.outer(np.sin(u), np.cos(v))
        ry = dxy_m * 10 * np.outer(np.sin(u), np.sin(v))
        rz = dz_m * 0.5 * np.outer(np.cos(u), np.ones_like(v))

        col = [C['sky'], C['green'], C['purple']][col_idx]
        ax.plot_surface(rx, ry, rz, alpha=0.6, color=col)

        # ejes
        ax.quiver(0, 0, -dz_m * 0.55, 0, 0, dz_m * 1.1, color='black',
                  arrow_length_ratio=0.15, lw=1.5)
        ax.text(0, 0, dz_m * 0.62, 'z', fontsize=10, ha='center')

        ax.set_title(title, fontsize=10)
        ax.set_xlabel('x'); ax.set_ylabel('y')
        lim = max(dxy_m * 12, dz_m * 0.6)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.view_init(elev=20, azim=-60)

    fig.suptitle('Forma de la PSF 3D según NorthOffset\n'
                 '(misma resolución $\\delta_z$ en todos los casos)',
                 fontsize=12, y=1.01)
    plt.tight_layout()
    save('fig09_psf_comparison.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 10: Curvas de resolución vs parámetros del sistema
# ─────────────────────────────────────────────────────────────────────────────
def fig10_resolution_curves():
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    c = 3e8; n2 = 2.0; f0 = 10e9; lam0 = c / f0
    psi0 = np.radians(58); cos_t0 = np.sqrt(1 - np.sin(psi0)**2 / n2**2)
    R0 = 188.7; B_perp = 47.17

    # ── δz vs B_perp ─────────────────────────────────────────────────────
    ax = axes[0]
    Bperp_arr = np.linspace(0, 100, 200)
    for B_hz, col, lab in [(25e6, C['sky'], '25 MHz'), (50e6, C['green'], '50 MHz'),
                            (100e6, C['orange'], '100 MHz')]:
        term_bw = n2 * B_hz * cos_t0
        term_tomo = c * Bperp_arr * np.sin(psi0) * np.cos(psi0) / (lam0 * R0 * n2 * cos_t0)
        Wz = term_bw + term_tomo
        dz = c / (2 * Wz)
        ax.plot(Bperp_arr, dz * 100, color=col, lw=2, label=f'$B={lab}$')
    ax.set_xlabel(r'$B_\perp$ [m]')
    ax.set_ylabel(r'$\delta_z$ [cm]')
    ax.set_title(r'$\delta_z$ vs. $B_\perp$', fontsize=11)
    ax.legend(fontsize=9); ax.grid(True)
    ax.set_ylim(0, 50)

    # ── δxy vs NorthOffset ───────────────────────────────────────────────
    ax = axes[1]
    dphi_arr = np.linspace(0, 175, 300)
    sin_t0 = np.sin(psi0) / n2
    for f0_hz, col, lab in [(5e9, C['sky'], '5 GHz'), (10e9, C['green'], '10 GHz'),
                              (20e9, C['orange'], '20 GHz')]:
        lam = c / f0_hz
        Rc = (4 * np.pi / lam) * np.sin(psi0) * np.abs(np.cos(np.radians(dphi_arr / 2)))
        Rc = np.where(Rc < 1, np.nan, Rc)
        dxy = 2 * 1.20 / Rc * 100  # cm
        ax.semilogy(dphi_arr, dxy, color=col, lw=2, label=f'$f_0={lab}$')
    ax.set_xlabel(r'NorthOffset $\Delta\phi$ [°]')
    ax.set_ylabel(r'$\delta_{xy}$ [cm]')
    ax.set_title(r'$\delta_{xy}$ vs. NorthOffset', fontsize=11)
    ax.legend(fontsize=9); ax.grid(True, which='both')
    ax.set_xlim(0, 175)

    # ── δz vs β (ángulo de inclinación) ────────────────────────────────
    ax = axes[2]
    B_hz = 50e6; B_helix = 47.17
    psi0_arr = np.radians(58)
    beta_arr = np.linspace(0, 90, 200)
    Bperp_b = B_helix * np.abs(np.cos(np.radians(beta_arr) - psi0_arr))
    term_bw = n2 * B_hz * cos_t0
    term_tomo = c * Bperp_b * np.sin(psi0) * np.cos(psi0) / (lam0 * R0 * n2 * cos_t0)
    Wz_b = term_bw + term_tomo
    dz_b = c / (2 * Wz_b)
    ax.plot(beta_arr, dz_b * 100, color=C['purple'], lw=2.5)
    ax.axvline(np.degrees(psi0_arr), color='red', ls='--', lw=2,
               label=f'Óptimo: $\\beta=\\psi_0={np.degrees(psi0_arr):.0f}°$')
    ax.scatter([np.degrees(psi0_arr)], [c / (2 * (term_bw + c * B_helix * np.sin(psi0) * np.cos(psi0) / (lam0 * R0 * n2 * cos_t0))) * 100],
               c='red', s=80, zorder=5)
    ax.set_xlabel(r'Ángulo de inclinación $\beta$ [°]')
    ax.set_ylabel(r'$\delta_z$ [cm]')
    ax.set_title(r'$\delta_z$ vs. $\beta$ ($B_{helix}=47$ m)', fontsize=11)
    ax.legend(fontsize=9); ax.grid(True)

    plt.tight_layout()
    save('fig10_resolution_curves.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 11: Off-axis — ángulo equivalente de Góes
# ─────────────────────────────────────────────────────────────────────────────
def fig11_offaxis_correction():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    # ── (a) geometría off-axis (vista superior) ───────────────────────────
    ax = axes[0]
    ax.set_aspect('equal')
    rho0 = 3.0; xP = 0.8

    # hélice (círculo en 2D)
    theta = np.linspace(0, 2 * np.pi, 200)
    ax.plot(rho0 * np.cos(theta), rho0 * np.sin(theta), color=C['sky'], lw=2, label='Hélice TX')

    # target off-axis
    ax.plot(xP, 0, '*', ms=16, color=C['target'], zorder=5)
    ax.text(xP + 0.1, 0.15, r'$P=(x_P,0,z_P)$', fontsize=10, color=C['target'])

    # near-range y far-range
    nr_x = rho0; nr_y = 0  # near-range (α=0)
    fr_x = -rho0; fr_y = 0  # far-range (α=π)
    ax.plot(nr_x, nr_y, 'o', ms=10, color=C['green'], zorder=5)
    ax.text(nr_x + 0.1, 0.2, 'Near-range\n($\\alpha=0$)', fontsize=9, color=C['green'])
    ax.plot(fr_x, fr_y, 'o', ms=10, color=C['orange'], zorder=5)
    ax.text(fr_x - 0.1, 0.2, 'Far-range\n($\\alpha=\\pi$)', ha='right', fontsize=9, color=C['orange'])

    # distancias
    ax.annotate('', xy=(xP, 0), xytext=(nr_x, 0),
                arrowprops=dict(arrowstyle='<->', color=C['green'], lw=2))
    ax.text((nr_x + xP) / 2, 0.2, r'$\rho_0-x_P$', ha='center', fontsize=10, color=C['green'])

    ax.annotate('', xy=(xP, 0), xytext=(fr_x, 0),
                arrowprops=dict(arrowstyle='<->', color=C['orange'], lw=2))
    ax.text((fr_x + xP) / 2, -0.25, r'$\rho_0+x_P$', ha='center', fontsize=10, color=C['orange'])

    # eje
    ax.plot(0, 0, '+', ms=12, color='black')
    ax.text(0.1, 0.1, 'Centro\nhélice', fontsize=9)

    ax.set_xlabel('x [m]'); ax.set_ylabel('y [m]')
    ax.set_title('(a) Vista superior: target off-axis\n(asimetría near/far range)', fontsize=11)
    ax.legend(fontsize=9); ax.grid(True)
    ax.set_xlim(-rho0 * 1.3, rho0 * 1.3)

    # ── (b) corrección Góes: ψ̃₀ < ψ₀ ─────────────────────────────────────
    ax = axes[1]
    z0 = 100; rho0_m = 160; xP_m = 20

    psi0 = np.degrees(np.arctan2(rho0_m, z0))
    psi_tilde = np.degrees(np.arctan2(rho0_m - xP_m, z0))

    xP_vals = np.linspace(0, 80, 200)
    psi0_arr = np.full_like(xP_vals, psi0)
    psi_tilde_arr = np.degrees(np.arctan2(rho0_m - xP_vals, z0))

    ax.plot(xP_vals, psi0_arr, '--', color=C['gray'], lw=2, label=r'$\psi_0$ (on-axis)')
    ax.plot(xP_vals, psi_tilde_arr, color=C['purple'], lw=2.5, label=r'$\tilde{\psi}_0$ (Góes)')
    ax.fill_between(xP_vals, psi_tilde_arr, psi0_arr, alpha=0.2, color=C['purple'])

    ax.axvline(xP_m, color=C['target'], ls=':', lw=2)
    ax.scatter([xP_m], [psi_tilde], c=C['purple'], s=80, zorder=5)
    ax.annotate(f'$x_P=20$ m\n$\\tilde{{\\psi}}_0={psi_tilde:.1f}°$\n(vs $\\psi_0={psi0:.1f}°$)',
                xy=(xP_m, psi_tilde), xytext=(xP_m + 12, psi_tilde - 3),
                arrowprops=dict(arrowstyle='->', color='black'),
                fontsize=10,
                bbox=dict(boxstyle='round', fc='lightyellow', ec='gray', alpha=0.85))

    ax.set_xlabel(r'Desplazamiento off-axis $x_P$ [m]')
    ax.set_ylabel(r'Look angle efectivo [°]')
    ax.set_title('(b) Corrección de Góes: ángulo equivalente\n'
                 r'$\tilde{\psi}_0 = \arctan\!\left(\frac{\rho_0-x_P}{z_0}\right)$',
                 fontsize=11)
    ax.legend(fontsize=10); ax.grid(True)

    plt.tight_layout()
    save('fig11_offaxis_correction.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 12: Resumen — curva de validación
# ─────────────────────────────────────────────────────────────────────────────
def fig12_validation():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # datos experimentales
    experiments = [
        ('ΔΦ=180°\non-axis', 0.211, 0.2112, 0.240, 0.240),
        ('ΔΦ=90°\non-axis', 0.211, 0.2112, 9.54e-3, 8.8e-3),
        ('ΔΦ=90°\noff-axis\n$x_P=20$ m', 0.1872, 0.2005, 9.96e-3, 8.8e-3),
    ]

    labels = [e[0] for e in experiments]
    dz_pred = [e[1] for e in experiments]
    dz_meas = [e[2] for e in experiments]
    dxy_pred = [e[3] for e in experiments]
    dxy_meas = [e[4] for e in experiments]

    x = np.arange(len(labels))
    w = 0.35

    # ── δz ───────────────────────────────────────────────────────────────
    ax = axes[0]
    b1 = ax.bar(x - w/2, [v * 100 for v in dz_pred], w, label='Predicción', color=C['sky'], alpha=0.85)
    b2 = ax.bar(x + w/2, [v * 100 for v in dz_meas], w, label='Simulación', color=C['green'], alpha=0.85)

    for bar, pred, meas in zip(x, dz_pred, dz_meas):
        err = abs(pred - meas) / meas * 100
        ax.text(bar, max(pred, meas) * 100 + 0.5, f'{err:.1f}%\nerror', ha='center', fontsize=8.5)

    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel(r'$\delta_z$ [cm]')
    ax.set_title(r'Resolución Vertical $\delta_z$', fontsize=12)
    ax.legend(fontsize=10); ax.grid(axis='y', alpha=0.4)
    ax.set_ylim(0, 30)

    # ── δxy ──────────────────────────────────────────────────────────────
    ax = axes[1]
    # solo los casos donde δxy es medible (índices 1 y 2)
    idx_valid = [1, 2]
    lv = [labels[i] for i in idx_valid]
    xv = np.arange(len(lv))
    b1 = ax.bar(xv - w/2, [dxy_pred[i] * 1e3 for i in idx_valid], w,
                label='Predicción', color=C['sky'], alpha=0.85)
    b2 = ax.bar(xv + w/2, [dxy_meas[i] * 1e3 for i in idx_valid], w,
                label='Simulación', color=C['green'], alpha=0.85)

    for bar_idx, i in enumerate(idx_valid):
        err = abs(dxy_pred[i] - dxy_meas[i]) / dxy_meas[i] * 100
        ax.text(bar_idx, max(dxy_pred[i], dxy_meas[i]) * 1e3 + 0.3,
                f'{err:.1f}%\nerror', ha='center', fontsize=8.5)

    ax.set_xticks(xv); ax.set_xticklabels(lv, fontsize=9)
    ax.set_ylabel(r'$\delta_{xy}$ [mm]')
    ax.set_title(r'Resolución Horizontal $\delta_{xy}$ (grilla fina, 2 mm)', fontsize=12)
    ax.legend(fontsize=10); ax.grid(axis='y', alpha=0.4)

    plt.tight_layout()
    save('fig12_validation.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print('Generando figuras...')
    fig01_radar_basics()
    fig02_aperture_real_vs_synth()
    fig03_circular_sar()
    fig04_helical_params()
    fig05_snell_geometry()
    fig06_kspace_vertical()
    fig07_northoffset_kspace()
    fig08_full_system()
    fig09_psf_comparison()
    fig10_resolution_curves()
    fig11_offaxis_correction()
    fig12_validation()
    print(f'\nFiguras guardadas en: {FIG_DIR}')
