"""Estilo comun de figuras (paleta categorica fija: mono=azul, bi=naranja, fusion=aqua)."""
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from . import config

C = {"mono": "#2a78d6", "bi": "#eb6834", "fus": "#1baf7a", "fusN": "#1baf7a",
     "s4": "#eda100", "s5": "#e87ba4", "s6": "#008300", "s7": "#6250d6", "s8": "#e34948",
     "ink": "#0b0b0b", "ink2": "#52514e", "grid": "#d9d8d4"}
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
LAB = {"mono": "Mono (Ref2)", "bi": "Bi", "fus": "Mono+Bi (2N)", "fusN": "Mono+Bi (N, mismo presupuesto)"}
SEQ = LinearSegmentedColormap.from_list("seq_azul", ["#fcfcfb", "#9cc3ee", "#2a78d6", "#0d3770", "#06172f"])
DR = 30.0   # rango dinamico comun de todas las imagenes [dB]

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9, "legend.fontsize": 8,
    "axes.edgecolor": C["ink2"], "axes.labelcolor": C["ink"], "xtick.color": C["ink2"],
    "ytick.color": C["ink2"], "axes.grid": True, "grid.color": C["grid"], "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0,
    "lines.markersize": 5, "figure.dpi": 110, "savefig.dpi": 200, "savefig.bbox": "tight",
    "legend.frameon": False, "image.cmap": "Blues",
})


def save(fig, name, outdir=None):
    """Guarda PDF (vectorial) y PNG en doc/figures/relatorio_1610 y en la carpeta del experimento."""
    for d in filter(None, [config.FIGDOC, outdir]):
        os.makedirs(d, exist_ok=True)
        fig.savefig(os.path.join(d, name + ".pdf"))
        fig.savefig(os.path.join(d, name + ".png"))
    plt.close(fig)


def db_img(A, ref):
    return 20 * np.log10(np.maximum(np.abs(A) / ref, 1e-12))


def show_slice(ax, X, Y, A, ref, title="", xlabel="x [m]", ylabel="y [m]", dr=DR):
    im = ax.pcolormesh(X, Y, db_img(A, ref).T, cmap=SEQ, vmin=-dr, vmax=0, shading="auto", rasterized=True)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return im
