"""
generate_figures_tx.py
Genera las figuras para doc/otimiza_tx.tex (estrategia de optimizacion Tx -> Rx).
Guarda en: hypotheses/figures/otx_fig*.pdf

Ejecutar desde la raiz del proyecto:
    python utils/generate_figures_tx.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(SCRIPT_DIR, '..', 'hypotheses', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams.update({
    'font.family':      'DejaVu Serif',
    'font.size':        11,
    'axes.labelsize':   12,
    'axes.titlesize':   12,
    'xtick.labelsize':  10,
    'ytick.labelsize':  10,
    'lines.linewidth':  2.0,
    'grid.alpha':       0.35,
    'figure.dpi':       150,
    'savefig.dpi':       200,
    'savefig.bbox':     'tight',
    'savefig.pad_inches': 0.1,
})

C = {
    'sky':    '#2c7bb6',
    'target': '#d7191c',
    'green':  '#1a9641',
    'purple': '#762a83',
    'orange': '#e08214',
    'gray':   '#636363',
}


def fresnel_tm_T(theta_i_deg, n1, n2):
    th_i = np.radians(theta_i_deg)
    sin_t = np.clip(n1 / n2 * np.sin(th_i), -1, 1)
    th_t = np.arcsin(sin_t)
    cos_i, cos_t = np.cos(th_i), np.cos(th_t)
    r = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t)
    return 1 - r ** 2


# ── Parametros reais do sistema (radarTx_espiral_plano_voo.json / Tabela 6.1) ──
n1, n2 = 1.0, 2.0
theta_B = np.degrees(np.arctan(n2 / n1))
f0 = 400e6
c = 3e8
lam0 = c / f0
B = 50e6
RaioMenor, RaioMaior = 147.5, 172.5
AltMenor, AltMaior = 80.0, 120.0
rho0_nom = (RaioMenor + RaioMaior) / 2
z0_nom = (AltMenor + AltMaior) / 2
B_helix = np.hypot(RaioMaior - RaioMenor, AltMaior - AltMenor)
psi0_nom = np.degrees(np.arctan(rho0_nom / z0_nom))
R0_nom = np.hypot(rho0_nom, z0_nom)


def Wz_R0fixed(psi0_deg, R0):
    psi0_r = np.radians(psi0_deg)
    cos_t0 = np.sqrt(np.clip(1 - (np.sin(psi0_r) / n2) ** 2, 1e-9, None))
    return n2 * B * cos_t0 + c * B_helix * np.sin(psi0_r) * np.cos(psi0_r) / (lam0 * R0 * n2 * cos_t0)


def Wz_z0fixed(psi0_deg, z0):
    psi0_r = np.radians(psi0_deg)
    R0 = z0 / np.cos(psi0_r)
    cos_t0 = np.sqrt(np.clip(1 - (np.sin(psi0_r) / n2) ** 2, 1e-9, None))
    return n2 * B * cos_t0 + c * B_helix * np.sin(psi0_r) * np.cos(psi0_r) / (lam0 * R0 * n2 * cos_t0)


def save(name):
    path = os.path.join(FIG_DIR, name)
    plt.savefig(path, bbox_inches='tight', pad_inches=0.12)
    plt.close()
    print(f'  Saved: {name}')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 1: The psi0 conflict -- T1(psi0) vs Wz(psi0) (delta_z), two constraint
#         scenarios (R0 fixed / z0 fixed), with the combined-cost knee.
# ─────────────────────────────────────────────────────────────────────────────
def fig01_psi0_tradeoff():
    psi0 = np.linspace(5, 85, 2000)
    T1 = fresnel_tm_T(psi0, n1, n2)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))

    for ax, Wz_fn, fixed_val, fixed_label, tag in [
            (axes[0], lambda p: Wz_R0fixed(p, R0_nom), R0_nom, f'$R_0={R0_nom:.0f}$ m fixed', 'R0'),
            (axes[1], lambda p: Wz_z0fixed(p, z0_nom), z0_nom, f'$z_0={z0_nom:.0f}$ m fixed', 'z0')]:
        Wz = Wz_fn(psi0)
        dz = c / (2 * Wz) * 100  # cm

        ax2 = ax.twinx()
        l1, = ax.plot(psi0, T1, color=C['sky'], lw=2.3, label='$T_1(\\psi_0)$ (Tx power transmittance)')
        l2, = ax2.plot(psi0, dz, color=C['target'], lw=2.3, ls='--', label='$\\delta_z(\\psi_0)$')

        v1 = ax.axvline(45.0, color=C['gray'], lw=1.4, ls=':', alpha=0.85,
                         label='$45°$ (tomographic optimum)')
        v2 = ax.axvline(theta_B, color=C['orange'], lw=1.6, ls=':', alpha=0.9,
                         label=f'$\\theta_B={theta_B:.1f}°$ (Brewster)')
        v3 = ax.axvline(psi0_nom, color=C['purple'], lw=1.6, ls='-.', alpha=0.9,
                         label=f'current system ($\\psi_0={psi0_nom:.0f}°$)')

        ax.set_xlabel('Look angle $\\psi_0$ [deg]')
        ax.set_ylabel('$T_1$ (power transmittance)', color=C['sky'])
        ax2.set_ylabel('$\\delta_z$ [cm]', color=C['target'])
        ax.tick_params(axis='y', labelcolor=C['sky'])
        ax2.tick_params(axis='y', labelcolor=C['target'])
        ax.set_ylim(0.5, 1.03)
        ax.set_title(fixed_label, fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend([l1, l2, v1, v2, v3],
                  [l1.get_label(), l2.get_label(), v1.get_label(), v2.get_label(), v3.get_label()],
                  fontsize=7.8, loc='lower left', framealpha=0.92)

    plt.suptitle('The $\\psi_0$ trade-off: Tx power (Brewster) vs. vertical resolution',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    save('otx_fig01_psi0_tradeoff.pdf')


# ─────────────────────────────────────────────────────────────────────────────
# FIG 2: Tx log-cost J_Tx(psi0;gamma) family and resulting psi0*(gamma) curve,
#         mirroring the existing Rx alpha-Pareto analysis (R0-fixed scenario).
# ─────────────────────────────────────────────────────────────────────────────
def fig02_gamma_sweep():
    psi0 = np.linspace(5, 85, 2000)
    T1 = fresnel_tm_T(psi0, n1, n2)
    Wz = Wz_R0fixed(psi0, R0_nom)
    dz = c / (2 * Wz) * 100

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))

    ax = axes[0]
    gammas_plot = [0.0, 0.1, 0.3, 0.5, 0.7, 1.0]
    cmap = plt.cm.plasma
    for i, gamma in enumerate(gammas_plot):
        J = -(1 - gamma) * np.log10(T1) + gamma * np.log10(dz)
        J = (J - J.min())  # shift for display only
        ax.plot(psi0, J, color=cmap(i / (len(gammas_plot) - 1)), lw=2,
                label=f'$\\gamma={gamma:.1f}$')
        idx = np.argmin(J)
        ax.plot(psi0[idx], J[idx], 'o', ms=6, color=cmap(i / (len(gammas_plot) - 1)))
    ax.set_xlabel('Look angle $\\psi_0$ [deg]')
    ax.set_ylabel('$J_{Tx}(\\psi_0;\\gamma)$ (shifted, a.u.)')
    ax.set_title('(a) $J_{Tx}$ vs. $\\psi_0$ for several $\\gamma$', fontsize=11)
    ax.legend(fontsize=8.5, ncol=2)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    gammas = np.linspace(0, 1, 61)
    psi0_star = []
    for gamma in gammas:
        J = -(1 - gamma) * np.log10(T1) + gamma * np.log10(dz)
        psi0_star.append(psi0[np.argmin(J)])
    psi0_star = np.array(psi0_star)
    ax.plot(gammas, psi0_star, color=C['sky'], lw=2.3)
    ax.axhline(psi0_nom, color=C['purple'], ls='-.', lw=1.5,
               label=f'current system ($\\psi_0={psi0_nom:.0f}°$)')
    ax.axhline(theta_B, color=C['orange'], ls=':', lw=1.5,
               label=f'Brewster ($\\theta_B={theta_B:.1f}°$)')
    ax.axhline(45.0, color=C['gray'], ls=':', lw=1.5, label='tomographic opt. ($45°$)')
    gamma_ref = 0.3
    J_ref = -(1 - gamma_ref) * np.log10(T1) + gamma_ref * np.log10(dz)
    psi0_ref = psi0[np.argmin(J_ref)]
    ax.plot([gamma_ref], [psi0_ref], 'D', ms=9, color=C['target'], zorder=5,
            label=f'$\\gamma=0.3 \\Rightarrow \\psi_0^*={psi0_ref:.1f}°$')
    ax.set_xlabel('Resolution weight $\\gamma$')
    ax.set_ylabel('Optimal look angle $\\psi_0^*$ [deg]')
    ax.set_title('(b) $\\psi_0^*(\\gamma)$ — $R_0$ fixed scenario', fontsize=11)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.suptitle('Stage-1 (Tx) log-cost minimization: knee search over $\\psi_0$',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    save('otx_fig02_gamma_sweep.pdf')


if __name__ == '__main__':
    print('Generando figuras (otimiza_tx)...')
    fig01_psi0_tradeoff()
    fig02_gamma_sweep()
    print(f'\nFiguras guardadas en: {FIG_DIR}')
