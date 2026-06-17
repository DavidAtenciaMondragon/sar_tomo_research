"""
generate_figures_radon.py
Genera todas las figuras para models/explicacion_radon.md
Guarda en: models/figures_radon/*.png

Ejecutar desde la raiz del proyecto:
    python utils/generate_figures_radon.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Rectangle
from skimage.transform import radon, iradon, rescale
from skimage.data import shepp_logan_phantom
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(SCRIPT_DIR, '..', 'models', 'figures_radon')
os.makedirs(FIG_DIR, exist_ok=True)

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
    'obj':    '#8c510a',
    'line':   '#d7191c',
    'axis':   '#2c2c2c',
    'rot':    '#2c7bb6',
    'fill':   '#fdae61',
    'green':  '#1a9641',
    'purple': '#762a83',
    'gray':   '#636363',
}

def save(name):
    path = os.path.join(FIG_DIR, name)
    plt.savefig(path, bbox_inches='tight', pad_inches=0.12)
    plt.close()
    print('saved', path)


# ──────────────────────────────────────────────────────────────────────────
# Fig 1b: Definicion de la recta L_{theta,s} a partir de theta y s
# ──────────────────────────────────────────────────────────────────────────
def fig_line_definition():
    fig, ax = plt.subplots(figsize=(6.5, 6.5))

    theta = np.deg2rad(35)
    L = 2.6
    s0 = 1.0

    n_dir = np.array([np.cos(theta), np.sin(theta)])     # direccion normal n-hat
    perp_dir = np.array([-np.sin(theta), np.cos(theta)]) # direccion a lo largo de la recta
    foot = s0 * n_dir                                     # pie de la perpendicular = s*n-hat

    # ejes x, y
    ax.annotate('', xy=(L, 0), xytext=(-0.3, 0),
                 arrowprops=dict(arrowstyle='-|>', color=C['axis'], lw=1.5))
    ax.annotate('', xy=(0, L), xytext=(0, -0.3),
                 arrowprops=dict(arrowstyle='-|>', color=C['axis'], lw=1.5))
    ax.text(L+0.08, -0.18, '$x$', fontsize=14)
    ax.text(-0.30, L+0.05, '$y$', fontsize=14)

    # direccion normal n-hat = (cos theta, sin theta), como linea punteada desde el origen
    ax.plot([0, L*n_dir[0]], [0, L*n_dir[1]], color=C['rot'], lw=1.3, ls='--',
            label=r'Normal $\hat{n}=(\cos\theta,\sin\theta)$')

    # vector s*n-hat desde el origen hasta el pie de la perpendicular
    ax.annotate('', xy=tuple(foot), xytext=(0, 0),
                 arrowprops=dict(arrowstyle='-|>', color=C['rot'], lw=2.2))
    mid = 0.5*foot
    ax.text(mid[0]-0.22, mid[1]+0.10, '$s$', color=C['rot'], fontsize=14)

    # angulo theta entre el eje x y la normal
    arc_t = np.linspace(0, theta, 30)
    ax.plot(0.45*np.cos(arc_t), 0.45*np.sin(arc_t), color=C['rot'], lw=1.2)
    ax.text(0.50, 0.10, r'$\theta$', color=C['rot'], fontsize=13)

    # la recta L_{theta,s}: perpendicular a n-hat, pasando por foot
    line_a = foot - 1.5*perp_dir
    line_b = foot + 1.5*perp_dir
    ax.plot([line_a[0], line_b[0]], [line_a[1], line_b[1]], color=C['line'], lw=2.5,
            label=r'$L_{\theta,s}:\ x\cos\theta+y\sin\theta=s$')

    # dos puntos genericos sobre la recta (cualquier punto de L satisface la ecuacion)
    for t, lbl, off in [(1.0, '$(x_1,y_1)$', (0.10, 0.05)),
                         (-0.8, '$(x_2,y_2)$', (0.10, -0.22))]:
        P = foot + t*perp_dir
        ax.plot(*P, 'o', color=C['gray'], ms=5, zorder=5)
        ax.text(P[0]+off[0], P[1]+off[1], lbl, fontsize=11, color=C['gray'])

    # pie de la perpendicular: punto de L mas cercano al origen, a distancia s sobre n-hat
    ax.plot(*foot, 'o', color=C['rot'], ms=5, zorder=5)

    ax.text(-L+0.05, -L+0.25,
            r'$L_{\theta,s}$ es perpendicular a $\hat{n}$ y lo cruza a'
            '\n'
            r'distancia $s$ del origen. Cualquier punto $(x,y)$'
            '\n'
            r'sobre $L_{\theta,s}$ cumple $x\cos\theta+y\sin\theta=s$',
            fontsize=10, color=C['gray'])

    ax.set_xlim(-L-0.3, L+0.5)
    ax.set_ylim(-L-0.3, L+0.6)
    ax.set_aspect('equal')
    ax.set_title(r'Definicion de la recta $L_{\theta,s}$', pad=18)
    ax.axis('off')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.06), frameon=False,
              fontsize=10, ncol=1)
    save('fig1b_definicion_recta.png')


# ──────────────────────────────────────────────────────────────────────────
# Fig 1: Geometria de la transformada de Radon (sistema rotado s-u, rayo)
# ──────────────────────────────────────────────────────────────────────────
def fig_geometry():
    fig, ax = plt.subplots(figsize=(6, 6))

    theta = np.deg2rad(35)
    L = 2.6

    # ejes originales x, y
    ax.annotate('', xy=(L, 0), xytext=(-L, 0),
                 arrowprops=dict(arrowstyle='-|>', color=C['axis'], lw=1.5))
    ax.annotate('', xy=(0, L), xytext=(0, -L),
                 arrowprops=dict(arrowstyle='-|>', color=C['axis'], lw=1.5))
    ax.text(L+0.05, -0.15, '$x$', fontsize=14)
    ax.text(-0.25, L+0.02, '$y$', fontsize=14)

    # ejes rotados s (a lo largo del rayo) y u (perpendicular, direccion de proyeccion)
    s_dir = np.array([np.cos(theta), np.sin(theta)])
    u_dir = np.array([-np.sin(theta), np.cos(theta)])

    ax.annotate('', xy=tuple(L*s_dir), xytext=tuple(-L*s_dir),
                 arrowprops=dict(arrowstyle='-|>', color=C['rot'], lw=1.8))
    ax.annotate('', xy=tuple(L*u_dir), xytext=tuple(-0.3*u_dir),
                 arrowprops=dict(arrowstyle='-|>', color=C['green'], lw=1.8))
    ax.text(*(L*s_dir + np.array([0.05, 0.05])), '$s$', color=C['rot'], fontsize=14)
    ax.text(*(L*u_dir + np.array([0.05, 0.05])), '$u$', color=C['green'], fontsize=14)

    # objeto f(x,y): un par de manchas elipticas
    obj1 = plt.matplotlib.patches.Ellipse((0.1, 0.05), 1.6, 1.0, angle=15,
                                           fc=C['obj'], alpha=0.35, ec=C['obj'])
    obj2 = plt.matplotlib.patches.Ellipse((0.5, -0.4), 0.6, 0.4, angle=-20,
                                           fc=C['obj'], alpha=0.55, ec=C['obj'])
    ax.add_patch(obj1)
    ax.add_patch(obj2)
    ax.text(0.0, 0.75, '$f(x,y)$', fontsize=13, color=C['obj'])

    # rayo (linea de integracion) a distancia s0 del origen, en direccion u
    s0 = 0.55
    p0 = s0 * u_dir
    ray_a = p0 - 1.6*s_dir
    ray_b = p0 + 1.6*s_dir
    ax.plot([ray_a[0], ray_b[0]], [ray_a[1], ray_b[1]], color=C['line'], lw=2.5,
            label='Rayo de integracion')

    # marcar distancia s0 desde origen perpendicular al rayo
    ax.annotate('', xy=tuple(p0), xytext=(0, 0),
                 arrowprops=dict(arrowstyle='-|>', color=C['gray'], lw=1.3, ls='--'))
    ax.text(0.18, 0.38, '$s_0$', color=C['gray'], fontsize=12)

    # angulo theta entre x y s
    arc_t = np.linspace(0, theta, 30)
    ax.plot(0.45*np.cos(arc_t), 0.45*np.sin(arc_t), color=C['rot'], lw=1.2)
    ax.text(0.55, 0.13, r'$\theta$', color=C['rot'], fontsize=13)

    # punto sobre el rayo
    ax.plot(*p0, 'o', color=C['gray'], ms=4)

    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_aspect('equal')
    ax.set_title('Geometria de un rayo de proyeccion')
    ax.axis('off')
    ax.legend(loc='lower right', frameon=False, fontsize=10)
    save('fig1_geometria.png')


# ──────────────────────────────────────────────────────────────────────────
# Fig 2: Punto fuente y su sinograma (curva sinusoidal)
# ──────────────────────────────────────────────────────────────────────────
def fig_point_sinogram():
    N = 256
    img = np.zeros((N, N))
    cx, cy = N//2 + 40, N//2 + 25   # punto desplazado del centro
    img[cy, cx] = 1.0

    # suavizar un poco para visualizacion (disco pequeno)
    yy, xx = np.mgrid[0:N, 0:N]
    r2 = (xx-cx)**2 + (yy-cy)**2
    img = np.exp(-r2/(2*3.0**2))

    thetas = np.linspace(0., 180., 180, endpoint=False)
    sino = radon(img, theta=thetas, circle=True)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))

    axes[0].imshow(img, cmap='gray', extent=(-N//2, N//2, -N//2, N//2), origin='lower')
    axes[0].set_title(r'Objeto: fuente puntual en $(x_0,y_0)$')
    axes[0].set_xlabel('$x$')
    axes[0].set_ylabel('$y$')

    im = axes[1].imshow(sino, cmap='inferno', aspect='auto',
                         extent=(thetas[0], thetas[-1], -sino.shape[0]//2, sino.shape[0]//2),
                         origin='lower')
    axes[1].set_title('Sinograma $p(\\theta, s)$')
    axes[1].set_xlabel(r'$\theta$ (grados)')
    axes[1].set_ylabel('$s$ (pixeles)')

    # superponer la curva analitica s = x0 cos(theta) + y0 sin(theta)
    x0 = (cx - N//2)
    y0 = (cy - N//2)
    th_rad = np.deg2rad(thetas)
    s_curve = x0*np.cos(th_rad) + y0*np.sin(th_rad)
    axes[1].plot(thetas, s_curve, color='cyan', lw=1.8, ls='--',
                  label=r'$s=x_0\cos\theta+y_0\sin\theta$')
    axes[1].legend(loc='upper right', fontsize=9)

    fig.colorbar(im, ax=axes[1], shrink=0.8, label='Amplitud')
    fig.tight_layout()
    save('fig2_punto_sinograma.png')


# ──────────────────────────────────────────────────────────────────────────
# Fig 3: Phantom de Shepp-Logan y su sinograma
# ──────────────────────────────────────────────────────────────────────────
def fig_phantom_sinogram():
    phantom = shepp_logan_phantom()
    phantom = rescale(phantom, scale=0.6, mode='reflect', channel_axis=None)

    thetas = np.linspace(0., 180., max(phantom.shape), endpoint=False)
    sino = radon(phantom, theta=thetas, circle=True)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))

    axes[0].imshow(phantom, cmap='gray')
    axes[0].set_title('Phantom Shepp-Logan: $f(x,y)$')
    axes[0].axis('off')

    im = axes[1].imshow(sino, cmap='inferno', aspect='auto',
                         extent=(thetas[0], thetas[-1], -sino.shape[0]//2, sino.shape[0]//2),
                         origin='lower')
    axes[1].set_title('Sinograma $p(\\theta,s) = \\mathcal{R}\\{f\\}$')
    axes[1].set_xlabel(r'$\theta$ (grados)')
    axes[1].set_ylabel('$s$ (pixeles)')
    fig.colorbar(im, ax=axes[1], shrink=0.8, label='Proyeccion')

    fig.tight_layout()
    save('fig3_phantom_sinograma.png')
    return phantom, sino, thetas


# ──────────────────────────────────────────────────────────────────────────
# Fig 4: Teorema de la Rebanada de Fourier (Fourier Slice Theorem)
# ──────────────────────────────────────────────────────────────────────────
def fig_fourier_slice(phantom):
    F = np.fft.fftshift(np.fft.fft2(phantom))
    mag = np.log1p(np.abs(F))

    N = phantom.shape[0]
    theta_deg = 35
    theta = np.deg2rad(theta_deg)

    # proyeccion a ese angulo y su FFT 1D
    sino_line = radon(phantom, theta=[theta_deg], circle=True)[:, 0]
    P = np.fft.fftshift(np.fft.fft(sino_line))
    freqs = np.fft.fftshift(np.fft.fftfreq(len(sino_line)))

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4))

    axes[0].plot(np.linspace(-N//2, N//2, len(sino_line)), sino_line, color=C['line'])
    axes[0].set_title(r'Proyeccion $p(\theta_0,s)$, $\theta_0=%d^\circ$' % theta_deg)
    axes[0].set_xlabel('$s$')
    axes[0].set_ylabel('Amplitud')

    axes[1].plot(freqs, np.abs(P), color=C['rot'])
    axes[1].set_title(r'$|P(\theta_0,\omega)| = |\mathrm{FT}_{1D}\{p(\theta_0,\cdot)\}|$')
    axes[1].set_xlabel(r'$\omega$ (frecuencia espacial)')
    axes[1].set_ylabel('Magnitud')

    axes[2].imshow(mag, cmap='gray', extent=(-N//2, N//2, -N//2, N//2), origin='lower')
    L = N//2 - 2
    axes[2].plot([-L*np.cos(theta), L*np.cos(theta)],
                  [-L*np.sin(theta), L*np.sin(theta)],
                  color=C['line'], lw=2.0,
                  label=r'Linea radial a $\theta_0$')
    axes[2].set_title(r'$|F(u,v)| = |\mathrm{FT}_{2D}\{f\}|$ (escala log)')
    axes[2].set_xlabel('$u$')
    axes[2].set_ylabel('$v$')
    axes[2].legend(loc='upper right', fontsize=9)

    fig.suptitle('Teorema de la rebanada de Fourier: '
                  r'$P(\theta,\omega) = F(\omega\cos\theta,\,\omega\sin\theta)$', y=1.03)
    fig.tight_layout()
    save('fig4_fourier_slice.png')


# ──────────────────────────────────────────────────────────────────────────
# Fig 5: Filtro rampa |omega| (dominio de frecuencia y ejemplo de ventana)
# ──────────────────────────────────────────────────────────────────────────
def fig_ramp_filter():
    N = 257
    omega = np.fft.fftshift(np.fft.fftfreq(N))
    ramp = np.abs(omega)

    # ventanas comunes aplicadas al filtro rampa
    hann = ramp * (0.5 + 0.5*np.cos(np.pi*omega/omega.max()))
    hamming = ramp * (0.54 + 0.46*np.cos(np.pi*omega/omega.max()))

    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    ax.plot(omega, ramp, color=C['line'], label=r'Ideal: $|\omega|$ (Ram-Lak)')
    ax.plot(omega, hann, color=C['rot'], ls='--', label='Hann')
    ax.plot(omega, hamming, color=C['green'], ls=':', label='Hamming')
    ax.set_xlabel(r'Frecuencia espacial $\omega$')
    ax.set_ylabel(r'$|H(\omega)|$')
    ax.set_title('Filtro de retroproyeccion filtrada (FBP)')
    ax.legend(frameon=False)
    ax.grid(True)
    fig.tight_layout()
    save('fig5_filtro_rampa.png')


# ──────────────────────────────────────────────────────────────────────────
# Fig 6: Reconstruccion - retroproyeccion simple vs filtrada
# ──────────────────────────────────────────────────────────────────────────
def fig_reconstruction(phantom, sino, thetas):
    recon_unfiltered = iradon(sino, theta=thetas, filter_name=None, circle=True)
    recon_filtered = iradon(sino, theta=thetas, filter_name='ramp', circle=True)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.4))

    axes[0].imshow(phantom, cmap='gray')
    axes[0].set_title('Original $f(x,y)$')
    axes[0].axis('off')

    axes[1].imshow(recon_unfiltered, cmap='gray')
    axes[1].set_title('Retroproyeccion simple\n(sin filtrar, borrosa)')
    axes[1].axis('off')

    axes[2].imshow(recon_filtered, cmap='gray')
    axes[2].set_title('Retroproyeccion filtrada\n(filtro rampa)')
    axes[2].axis('off')

    fig.tight_layout()
    save('fig6_reconstruccion.png')


# ──────────────────────────────────────────────────────────────────────────
# Fig 7: Retroproyeccion acumulativa - "estrella" de pocas proyecciones
# ──────────────────────────────────────────────────────────────────────────
def fig_backprojection_buildup():
    N = 200
    img = np.zeros((N, N))
    yy, xx = np.mgrid[0:N, 0:N]
    cx, cy = N//2, N//2
    r2 = (xx-cx)**2 + (yy-cy)**2
    img = np.exp(-r2/(2*4.0**2))

    n_list = [1, 4, 16, 180]
    fig, axes = plt.subplots(1, len(n_list), figsize=(13, 3.6))

    for ax, n in zip(axes, n_list):
        thetas = np.linspace(0., 180., n, endpoint=False)
        sino = radon(img, theta=thetas, circle=True)
        recon = iradon(sino, theta=thetas, filter_name=None, circle=True)
        ax.imshow(recon, cmap='gray')
        ax.set_title(f'{n} proyeccion(es)')
        ax.axis('off')

    fig.suptitle('Retroproyeccion simple de un punto: superposicion de "estrellas"', y=1.02)
    fig.tight_layout()
    save('fig7_retroproyeccion_acumulada.png')


if __name__ == '__main__':
    fig_line_definition()
    fig_geometry()
    fig_point_sinogram()
    phantom, sino, thetas = fig_phantom_sinogram()
    fig_fourier_slice(phantom)
    fig_ramp_filter()
    fig_reconstruction(phantom, sino, thetas)
    fig_backprojection_buildup()
    print('Listo. Figuras en', FIG_DIR)
