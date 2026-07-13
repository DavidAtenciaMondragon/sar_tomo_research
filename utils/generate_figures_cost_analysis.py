"""
generate_figures_cost_analysis.py
Análisis comparativo de funciones de costo para el optimizador
energía/resolución del SAR biestático helicoidal.

Genera 5 figuras en utils/:
  fig_cf01_landscapes.pdf      — Paisajes de costo en el plano (x_Rx, y_Rx)
  fig_cf02_pareto_fronts.pdf   — Fronteras de Pareto en espacio de objetivos
  fig_cf03_alpha_sensitivity.pdf — Posición óptima del Rx vs. parámetro α
  fig_cf04_gradient_norm.pdf   — Norma del gradiente para cada función
  fig_cf05_brewster_penalty.pdf — Efecto de la regularización hacia θ_B

Convenciones:
  - Alvo en (0, 0, -5 m); Tx fijo en (160, 0, 100 m)
  - Rx varía en grilla cuadrada a altura z_Rx = 100 m
  - n1 = 1 (aire), n2 = 2 (suelo seco, ε_r = 4)
  - f0 = 10 GHz, B = 50 MHz (banda X típica)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import os

# ── directorios ───────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = SCRIPT_DIR          # guardar directamente en utils/
os.makedirs(FIG_DIR, exist_ok=True)

# ── estilo global (académico) ─────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':       'DejaVu Serif',
    'font.size':         10,
    'axes.labelsize':    11,
    'axes.titlesize':    11,
    'xtick.labelsize':   9,
    'ytick.labelsize':   9,
    'lines.linewidth':   1.8,
    'grid.alpha':        0.3,
    'figure.dpi':        150,
    'savefig.dpi':       200,
    'savefig.bbox':      'tight',
    'savefig.pad_inches': 0.1,
    'axes.spines.top':   False,
    'axes.spines.right': False,
})

PALETTE = ['#1f4e79', '#c0392b', '#27ae60', '#8e44ad', '#d35400']

def save(name):
    path = os.path.join(FIG_DIR, name)
    plt.savefig(path, bbox_inches='tight', pad_inches=0.12)
    plt.close()
    print(f'  Guardada: {name}')

# ══════════════════════════════════════════════════════════════════════════════
# PARÁMETROS DEL SISTEMA
# ══════════════════════════════════════════════════════════════════════════════
f0      = 10e9          # Hz — frecuencia portadora
c_light = 299_792_458.0 # m/s
lam     = c_light / f0  # m — longitud de onda

B       = 50e6          # Hz — ancho de banda chirp
n1      = 1.0           # índice aire
n2      = 2.0           # índice suelo (eps_r = 4)
theta_B = np.arctan(n2 / n1)   # ángulo de Brewster (rad)

z_rx    = 100.0         # m — altitud del Rx (y del Tx)
z_tg    = 5.0           # m — profundidad del blanco (positiva)

# Tx en posición fija sobre la hélice
rho_Tx  = 160.0         # m
x_Tx    = rho_Tx        # Tx en (160, 0, 100)
y_Tx    = 0.0
z_Tx    = z_rx

# Parámetros de la hélice (condición óptima β = ψ₀)
psi_Tx   = np.arctan(rho_Tx / z_Tx)  # ≈ 58°
B_helix  = 47.17                       # m
beta_hel = psi_Tx                      # condición óptima

# Ecuación de radar (parámetros de referencia)
Pt, sigma = 1.0, 1.0    # W, m²

# ══════════════════════════════════════════════════════════════════════════════
# FUNCIONES AUXILIARES
# ══════════════════════════════════════════════════════════════════════════════
def snell_bisect(d, z_air, z_soil, n_air=1.0, n_soil=2.0, n_iter=60):
    """Resolve el parámetro de rayo p = sin(θ_aire) por bisección vectorizada."""
    p_lo = np.zeros_like(d, dtype=float)
    p_hi = np.full_like(d, 1.0 - 1e-10, dtype=float)
    for _ in range(n_iter):
        p = 0.5 * (p_lo + p_hi)
        c1 = np.sqrt(np.clip(1.0 - p**2,           1e-14, None))
        c2 = np.sqrt(np.clip(1.0 - (p/n_soil)**2,  1e-14, None))
        d_p = z_air * p / c1 + z_soil * (p / n_soil) / c2
        go_up = d_p < d
        p_lo  = np.where(go_up, p, p_lo)
        p_hi  = np.where(go_up, p_hi, p)
    return 0.5 * (p_lo + p_hi)


def ray_geometry(xr, yr, zr, xt, yt, zt_depth):
    """
    Geometría del rayo refractado desde el receptor (xr,yr,zr)
    hasta el blanco (xt, yt, -zt_depth) pasando por la interface z=0.
    Retorna: R_total (m), theta_inc_air (rad) — vectorizados.
    """
    d = np.sqrt((xr - xt)**2 + (yr - yt)**2)
    p = snell_bisect(d, zr, zt_depth)
    c1 = np.sqrt(np.clip(1.0 - p**2,           1e-14, None))
    c2 = np.sqrt(np.clip(1.0 - (p/n2)**2,      1e-14, None))
    R_air  = zr      / c1
    R_soil = zt_depth / c2
    theta  = np.arcsin(np.clip(p, 0.0, 1.0 - 1e-10))
    return R_air + R_soil, theta


def fresnel_T(theta_air):
    """Transmitancia TM de potencia — vectorizada."""
    sin_t = np.clip(np.sin(theta_air) / n2, -1.0, 1.0)
    cos_t = np.sqrt(np.clip(1.0 - sin_t**2, 1e-14, None))
    cos_i = np.cos(theta_air)
    denom = n2 * cos_i + n1 * cos_t
    t_TM  = np.where(np.abs(denom) > 1e-14,
                     2.0 * n1 * cos_i / denom, 0.0)
    T     = (n2 * cos_t / (n1 * cos_i + 1e-14)) * t_TM**2
    return np.clip(T, 0.0, 1.0)


def bistatic_resolution(xr, yr, zr, xt_pos, yt_pos):
    """
    Resolución biestática analítica (δ_xy, δ_z) para Rx en (xr,yr,zr).
    Tx en posición global (x_Tx, y_Tx, z_Tx). Blanco en (0,0).
    """
    rho_r = np.sqrt((xr - xt_pos)**2 + (yr - yt_pos)**2)
    rho_r = np.maximum(rho_r, 1.0)          # evitar ψ₀=0
    psi0  = np.arctan(rho_r / zr)
    R0    = np.sqrt(rho_r**2 + zr**2)

    sin_t = np.clip(np.sin(psi0) / n2, 0.0, 1.0)
    cos_t = np.sqrt(np.clip(1.0 - sin_t**2, 1e-12, None))
    cos_t = np.maximum(cos_t, 1e-9)

    B_perp = B_helix * np.abs(np.cos(beta_hel - psi0))
    Wz     = (n2 * B * cos_t
              + c_light * B_perp * np.sin(psi0) * np.cos(psi0)
              / (lam * R0 * n2 * cos_t + 1e-14))
    Wz     = np.maximum(Wz, 1.0)
    dz     = c_light / (2.0 * Wz)

    az_Tx  = np.arctan2(y_Tx - yt_pos, x_Tx - xt_pos)
    az_Rx  = np.arctan2(yr   - yt_pos, xr   - xt_pos)
    DPhi   = (az_Rx - az_Tx + np.pi) % (2*np.pi) - np.pi
    cos_dp = np.abs(np.cos(DPhi / 2.0))
    cos_dp = np.maximum(cos_dp, 1e-9)

    sin_p0 = np.maximum(np.sin(psi0), 1e-9)
    dxy    = 0.60 * lam / (np.pi * sin_p0 * cos_dp)
    return dxy, dz

# ══════════════════════════════════════════════════════════════════════════════
# GRILLA DE POSICIONES DEL Rx
# ══════════════════════════════════════════════════════════════════════════════
N  = 90
XG = np.linspace(-270, 270, N)
YG = np.linspace(-270, 270, N)
XX, YY = np.meshgrid(XG, YG)

# Precomputa Tx → blanco (constante para Tx fijo)
R_T, theta_T = ray_geometry(x_Tx, y_Tx, z_Tx, 0.0, 0.0, z_tg)
T1 = float(fresnel_T(np.array([theta_T]))[0])

# Rx → blanco
R_R, theta_R_grid = ray_geometry(XX, YY, z_rx, 0.0, 0.0, z_tg)
T2_grid = fresnel_T(theta_R_grid)

# Potencia recibida (ecuación del radar biestático con refracción)
Pr_grid = (Pt * sigma * T1 * T2_grid * lam**2) / \
          ((4*np.pi)**3 * R_T**2 * R_R**2)
Pr_grid = np.maximum(Pr_grid, 1e-40)

# Resolución
dxy_grid, dz_grid = bistatic_resolution(XX, YY, z_rx, 0.0, 0.0)
# Enmascarar puntos degenerados
bad = ~np.isfinite(dxy_grid) | ~np.isfinite(dz_grid) | \
      (dxy_grid <= 0) | (dz_grid <= 0)
dxy_grid[bad] = np.nan
dz_grid[bad]  = np.nan

# Términos de objetivo
f1 = -np.log10(Pr_grid)                          # minimizar → maximizar Pr
f2 = np.log10(dxy_grid * dz_grid)                # minimizar → minimizar producto
f2 = np.where(bad, np.nan, f2)

# Rangos para normalización (percentil robusto)
f1_min, f1_max = np.nanpercentile(f1, 2), np.nanpercentile(f1, 98)
f2_min, f2_max = np.nanpercentile(f2, 2), np.nanpercentile(f2, 98)
f1_rng = f1_max - f1_min
f2_rng = f2_max - f2_min

f1n = np.clip((f1 - f1_min) / f1_rng, 0.0, 1.0)   # normalizado [0,1]
f2n = np.clip((f2 - f2_min) / f2_rng, 0.0, 1.0)

# ── Definición de las 5 funciones de costo ───────────────────────────────────
def J1(a): return (1-a)*f1 + a*f2
def J2(a): return (1-a)*f1n + a*f2n          # lineal normalizada
def J3(a):                                    # Tchebyshev
    with np.errstate(invalid='ignore'):
        return np.maximum((1-a)*f1n, a*f2n)
def J4():  return np.sqrt(f1n**2 + f2n**2)   # distancia euclídea (sin α)
def J5(a, gamma):                             # J1 + penalización Brewster
    pen = (theta_R_grid - theta_B)**2
    return (1-a)*f1 + a*f2 + gamma * pen

alpha_ref = 0.30   # valor de referencia para paisajes

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 1 — Paisajes de optimización
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 3, figsize=(13, 8.5), sharex=True, sharey=True)
fig.suptitle(r'Paisajes de las funciones de costo en el plano $(x_{Rx},\,y_{Rx})$'
             f'\n'
             r'($\alpha = 0{,}30$, $\gamma = 5\,\mathrm{rad}^{-2}$, Tx en estrella, blanco en cruz)',
             fontsize=10.5, y=1.01)

cost_maps = [
    (f1,              r'$f_1 = -\log_{10}\,P_r$',             'viridis'),
    (f2,              r'$f_2 = \log_{10}(\delta_{xy}\delta_z)$', 'viridis'),
    (J1(alpha_ref),   r'$J_1$ (log ponderada)',                'plasma_r'),
    (J2(alpha_ref),   r'$J_2$ (lineal normalizada)',           'plasma_r'),
    (J3(alpha_ref),   r'$J_3$ (Tchebyshev)',                   'plasma_r'),
    (J5(alpha_ref, 5.0), r'$J_5$ (log + penalización $\theta_B$)', 'plasma_r'),
]

lbl = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
for ax, (data, title, cmap), lb in zip(axes.flat, cost_maps, lbl):
    valid = np.isfinite(data)
    vmin, vmax = np.nanpercentile(data[valid], 5), np.nanpercentile(data[valid], 95)
    im = ax.contourf(XX, YY, data, levels=30, cmap=cmap,
                     vmin=vmin, vmax=vmax, extend='both')
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, format='%.1f')

    # Mínimo de la función (posición óptima del Rx)
    if valid.any():
        idx = np.nanargmin(data)
        iy, ix = np.unravel_index(idx, data.shape)
        ax.plot(XX[iy, ix], YY[iy, ix], 'o', color='white',
                ms=7, mec='black', mew=1.2, zorder=5)

    ax.plot(x_Tx, y_Tx, '*', color='gold',  ms=11, mec='black', mew=0.8, zorder=6)
    ax.plot(0,    0,    '+', color='white',  ms=9,  mew=2.0,    zorder=6)
    ax.set_title(f'{lb}  {title}', pad=4)
    ax.set_aspect('equal')
    ax.grid(True, lw=0.4)

for ax in axes[1, :]:
    ax.set_xlabel('$x_{Rx}$ [m]')
for ax in axes[:, 0]:
    ax.set_ylabel('$y_{Rx}$ [m]')

legend_elems = [
    Line2D([0],[0], marker='*', color='gold', ms=9, mec='black', lw=0, label='Tx'),
    Line2D([0],[0], marker='+', color='gray', ms=9,  mew=2,      lw=0, label='Blanco'),
    Line2D([0],[0], marker='o', color='white',ms=7, mec='black', lw=0, label='Rx óptimo'),
]
fig.legend(handles=legend_elems, loc='lower center', ncol=3, fontsize=9,
           bbox_to_anchor=(0.5, -0.03))
plt.tight_layout()
save('fig_cf01_landscapes.pdf')

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 2 — Fronteras de Pareto en el espacio de objetivos
# ══════════════════════════════════════════════════════════════════════════════
alphas = np.linspace(0.0, 1.0, 60)

def pareto_curve(cost_fn):
    f1_pts, f2_pts = [], []
    for a in alphas:
        cost = cost_fn(a)
        if not np.any(np.isfinite(cost)):
            continue
        idx = np.nanargmin(cost)
        iy, ix = np.unravel_index(idx, cost.shape)
        f1_pts.append(np.nanmean(f1[max(0,iy-1):iy+2, max(0,ix-1):ix+2]))
        f2_pts.append(np.nanmean(f2[max(0,iy-1):iy+2, max(0,ix-1):ix+2]))
    return np.array(f1_pts), np.array(f2_pts)

f1_J1, f2_J1 = pareto_curve(J1)
f1_J2, f2_J2 = pareto_curve(J2)
f1_J3, f2_J3 = pareto_curve(J3)

# J4: punto único (sin α)
cost4 = J4()
idx4  = np.nanargmin(cost4)
iy4, ix4 = np.unravel_index(idx4, cost4.shape)
f1_J4, f2_J4 = float(f1[iy4, ix4]), float(f2[iy4, ix4])

# Punto ideal (mínimo teórico de ambos objetivos simultáneamente)
f1_ideal = np.nanmin(f1)
f2_ideal = np.nanmin(f2)

fig, ax = plt.subplots(figsize=(7.5, 5.5))
ax.plot(f1_J1, f2_J1, '-o', color=PALETTE[0], ms=3.5, lw=1.6,
        label=r'$J_1$ — log ponderada')
ax.plot(f1_J2, f2_J2, '-s', color=PALETTE[1], ms=3.5, lw=1.6,
        label=r'$J_2$ — lineal normalizada')
ax.plot(f1_J3, f2_J3, '-^', color=PALETTE[2], ms=3.5, lw=1.6,
        label=r'$J_3$ — Tchebyshev')
ax.plot(f1_J4, f2_J4, 'D', color=PALETTE[3], ms=9, mec='black', mew=0.8, zorder=5,
        label=r'$J_4$ — distancia eucl. (sin $\alpha$)')
ax.plot(f1_ideal, f2_ideal, 'k*', ms=12, label='Punto ideal', zorder=6)

# Anotar α en J1
for i, a in enumerate(alphas):
    if abs(a - 0.30) < 0.01 or abs(a - 0.60) < 0.01:
        ax.annotate(fr'$\alpha={a:.2f}$', xy=(f1_J1[i], f2_J1[i]),
                    xytext=(6, 4), textcoords='offset points', fontsize=8)

ax.set_xlabel(r'$f_1 = -\log_{10}\,P_r$ (pérdida de potencia)', fontsize=11)
ax.set_ylabel(r'$f_2 = \log_{10}(\delta_{xy}\cdot\delta_z)\;\mathrm{[m^2]}$', fontsize=11)
ax.set_title('Fronteras de Pareto en el espacio de objetivos', fontsize=11)
ax.legend(fontsize=9, loc='upper right')
ax.grid(True, lw=0.4)
plt.tight_layout()
save('fig_cf02_pareto_fronts.pdf')

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 3 — Sensibilidad al parámetro α
# ══════════════════════════════════════════════════════════════════════════════
def opt_position(cost_fn, alphas):
    rho_opt, phi_opt, dxy_opt, dz_opt = [], [], [], []
    for a in alphas:
        cost = cost_fn(a)
        if not np.any(np.isfinite(cost)):
            rho_opt.append(np.nan); phi_opt.append(np.nan)
            dxy_opt.append(np.nan); dz_opt.append(np.nan)
            continue
        idx = np.nanargmin(cost)
        iy, ix = np.unravel_index(idx, cost.shape)
        xo, yo = XX[iy, ix], YY[iy, ix]
        rho_opt.append(np.sqrt(xo**2 + yo**2))
        phi_opt.append(np.degrees(np.arctan2(yo, xo)))
        dxy_opt.append(np.nanmean(dxy_grid[max(0,iy-1):iy+2, max(0,ix-1):ix+2]))
        dz_opt.append(np.nanmean(dz_grid[max(0,iy-1):iy+2, max(0,ix-1):ix+2]))
    return (np.array(rho_opt), np.array(phi_opt),
            np.array(dxy_opt), np.array(dz_opt))

alphas_fine = np.linspace(0.0, 1.0, 80)
rho1, phi1, dxy1, dz1 = opt_position(J1, alphas_fine)
rho2, phi2, dxy2, dz2 = opt_position(J2, alphas_fine)
rho3, phi3, dxy3, dz3 = opt_position(J3, alphas_fine)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

ax = axes[0]
ax.plot(alphas_fine, rho1, color=PALETTE[0], lw=1.8, label=r'$J_1$')
ax.plot(alphas_fine, rho2, color=PALETTE[1], lw=1.8, ls='--', label=r'$J_2$')
ax.plot(alphas_fine, rho3, color=PALETTE[2], lw=1.8, ls=':', label=r'$J_3$')
ax.axvline(0.3, color='gray', lw=0.8, ls='--', alpha=0.7, label=r'$\alpha=0{,}30$ (ref.)')
ax.set_xlabel(r'$\alpha$ (peso de resolución)', fontsize=11)
ax.set_ylabel(r'$\rho_{Rx}^{*}$ [m] (distancia radial al blanco)', fontsize=11)
ax.set_title('(a)  Distancia radial del Rx óptimo vs. $\\alpha$', fontsize=10)
ax.legend(fontsize=9)
ax.grid(True, lw=0.4)
ax.set_xlim(0, 1)

ax = axes[1]
ax.plot(alphas_fine, phi1, color=PALETTE[0], lw=1.8, label=r'$J_1$')
ax.plot(alphas_fine, phi2, color=PALETTE[1], lw=1.8, ls='--', label=r'$J_2$')
ax.plot(alphas_fine, phi3, color=PALETTE[2], lw=1.8, ls=':', label=r'$J_3$')
ax.axhline(np.degrees(np.arctan2(y_Tx, x_Tx)), color='gold',
           lw=1.0, ls='--', label='Azimut Tx')
ax.axvline(0.3, color='gray', lw=0.8, ls='--', alpha=0.7)
ax.set_xlabel(r'$\alpha$ (peso de resolución)', fontsize=11)
ax.set_ylabel(r'$\phi_{Rx}^{*}$ [°] (azimut del Rx óptimo)', fontsize=11)
ax.set_title('(b)  Azimut del Rx óptimo vs. $\\alpha$', fontsize=10)
ax.legend(fontsize=9)
ax.grid(True, lw=0.4)
ax.set_xlim(0, 1)

plt.tight_layout()
save('fig_cf03_alpha_sensitivity.pdf')

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 4 — Norma del gradiente (condicionamiento)
# ══════════════════════════════════════════════════════════════════════════════
# Calcula ‖∇J‖ numéricamente en la grilla
dx = XG[1] - XG[0]

def grad_norm(cost):
    gy, gx = np.gradient(np.where(np.isfinite(cost), cost, np.nan),
                         dx, dx)
    return np.sqrt(gx**2 + gy**2)

gn1 = grad_norm(J1(alpha_ref))
gn2 = grad_norm(J2(alpha_ref))
gn3 = grad_norm(J3(alpha_ref))

# Perfil a lo largo de la línea horizontal y_Rx = 0
mid = N // 2

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Panel izquierdo: norma del gradiente a lo largo de y=0
ax = axes[0]
valid_mask = np.isfinite(gn1[mid,:]) & np.isfinite(gn2[mid,:]) & \
             np.isfinite(gn3[mid,:])
xs = XG[valid_mask]
for gn, col, lbl in zip([gn1, gn2, gn3], PALETTE[:3],
                         [r'$J_1$ (log pond.)', r'$J_2$ (lineal norm.)',
                          r'$J_3$ (Tchebyshev)']):
    g = gn[mid, valid_mask]
    g_norm = g / np.nanmax(g)   # normalizar para comparación
    ax.plot(xs, g_norm, color=col, lw=1.8, label=lbl)

ax.set_xlabel(r'$x_{Rx}$ [m]  (perfil $y_{Rx}=0$)', fontsize=11)
ax.set_ylabel(r'$\|\nabla J\|$ normalizado', fontsize=11)
ax.set_title('(a)  Norma del gradiente a lo largo de $y_{Rx}=0$', fontsize=10)
ax.legend(fontsize=9)
ax.grid(True, lw=0.4)

# Panel derecho: mapa de ‖∇J₁‖ con isolíneas de costo
ax = axes[1]
cf = ax.contourf(XX, YY, gn1, levels=30, cmap='YlOrRd')
plt.colorbar(cf, ax=ax, fraction=0.046, pad=0.04, label=r'$\|\nabla J_1\|$')
cs = ax.contour(XX, YY, J1(alpha_ref), levels=15, colors='white',
                linewidths=0.5, alpha=0.5)
ax.clabel(cs, fmt='%.1f', fontsize=7, inline=True)
ax.plot(x_Tx, y_Tx, '*', color='gold', ms=11, mec='black', mew=0.8, zorder=5)
ax.plot(0, 0, '+', color='black', ms=9, mew=2, zorder=5)
ax.set_xlabel('$x_{Rx}$ [m]', fontsize=11)
ax.set_ylabel('$y_{Rx}$ [m]', fontsize=11)
ax.set_title(r'(b)  Mapa de $\|\nabla J_1\|$ con isolíneas de $J_1$', fontsize=10)
ax.set_aspect('equal')
ax.grid(True, lw=0.4)

plt.tight_layout()
save('fig_cf04_gradient_norm.pdf')

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 5 — Efecto de la penalización de Brewster (J₅)
# ══════════════════════════════════════════════════════════════════════════════
gammas = np.concatenate([np.linspace(0.0, 1.0, 20),
                          np.linspace(1.0, 20.0, 30)])

theta_opt_list  = []
Pr_loss_list    = []
delta_prod_list = []

for gm in gammas:
    cost5 = J5(alpha_ref, gm)
    if not np.any(np.isfinite(cost5)):
        theta_opt_list.append(np.nan)
        Pr_loss_list.append(np.nan)
        delta_prod_list.append(np.nan)
        continue
    idx = np.nanargmin(cost5)
    iy, ix = np.unravel_index(idx, cost5.shape)
    theta_opt_list.append(np.degrees(theta_R_grid[iy, ix]))
    Pr_loss_list.append(float(f1[iy, ix]))
    delta_prod_list.append(float(f2[iy, ix]))

theta_opts  = np.array(theta_opt_list)
Pr_losses   = np.array(Pr_loss_list)
delta_prods = np.array(delta_prod_list)

theta_B_deg = np.degrees(theta_B)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

ax = axes[0]
ax.plot(gammas, theta_opts, color=PALETTE[4], lw=2.0, label=r'$\theta_{R}^{*}(\gamma)$')
ax.axhline(theta_B_deg, color='black', lw=1.0, ls='--',
           label=fr'$\theta_B = {theta_B_deg:.1f}°$')
ax.set_xlabel(r'$\gamma$ [rad$^{-2}$] — intensidad de la penalización', fontsize=11)
ax.set_ylabel(r'$\theta_R^{*}$ [°] — ángulo de incidencia en Rx óptimo', fontsize=11)
ax.set_title(r'(a)  Convergencia de $\theta_R^*$ hacia $\theta_B$ al aumentar $\gamma$',
             fontsize=10)
ax.legend(fontsize=9)
ax.grid(True, lw=0.4)

ax = axes[1]
l1, = ax.plot(gammas, Pr_losses,   color=PALETTE[0], lw=1.8,
              label=r'$-\log_{10}P_r$ (pérdida de potencia)')
ax.set_xlabel(r'$\gamma$ [rad$^{-2}$]', fontsize=11)
ax.set_ylabel(r'$-\log_{10}\,P_r$ [dB equiv.]', color=PALETTE[0], fontsize=11)
ax.tick_params(axis='y', labelcolor=PALETTE[0])

ax2 = ax.twinx()
l2, = ax2.plot(gammas, delta_prods, color=PALETTE[1], lw=1.8, ls='--',
               label=r'$\log_{10}(\delta_{xy}\delta_z)$ (resolución)')
ax2.set_ylabel(r'$\log_{10}(\delta_{xy}\delta_z)\;\mathrm{[m^2]}$',
               color=PALETTE[1], fontsize=11)
ax2.tick_params(axis='y', labelcolor=PALETTE[1])

ax.set_title(r'(b)  Costo de la penalización: pérdida de potencia y resolución vs. $\gamma$',
             fontsize=10)
lines = [l1, l2]
ax.legend(lines, [l.get_label() for l in lines], fontsize=9, loc='center right')
ax.grid(True, lw=0.4)
ax2.spines['right'].set_visible(True)

plt.tight_layout()
save('fig_cf05_brewster_penalty.pdf')

print('\nTodas las figuras generadas correctamente en utils/')
