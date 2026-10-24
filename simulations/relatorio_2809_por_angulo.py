"""Desglose por ángulo de la comparación monoestático vs biestático especular.

Post-proceso de las salidas de verifica_ganancia_brewster.m (una corrida con los
cinco objetos) para el relatorio doc/relatorio_2809.tex:
  - figura fig_resumen_por_angulo.png: diferencia B-M por posición y objeto en
    orientación nominal y en el Monte Carlo (IC 95 %), y Pr{B>M} por posición;
  - diferencia energética por banda de elevación con IC 95 % bootstrap pareado
    (mismo remuestreo que el agregado global de verifica_ganancia_brewster.m);
  - filas LaTeX de las tablas por objeto y por posición del apéndice.

Uso: python relatorio_2809_por_angulo.py [carpeta_corrida]
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
CORRIDA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    AQUI, "io", "verifica_ganancia_brewster", "20260928_140704_objetos")
STAMP = os.path.basename(CORRIDA).split("_objetos")[0]
FIGDIR = os.path.join(AQUI, "..", "doc", "figures", "relatorio_2809")

OBJETOS = ["cubo", "caja", "cilindro", "esfera", "cono"]
ETIQ = {"cubo": "Cubo", "caja": "Caja", "cilindro": "Cilindro", "esfera": "Esfera", "cono": "Cono"}
COLOR = {"cubo": "#0072BD", "caja": "#D95319", "cilindro": "#EDB120", "esfera": "#7E2F8E", "cono": "#77AC30"}
MARCA = {"cubo": "s", "caja": "D", "cilindro": "o", "esfera": "^", "cono": "v"}
BANDAS = [(20, 35), (40, 60), (65, 85)]   # elevaciones [°], ambos extremos incluidos
EPS_B = 26.31                              # elevación de Brewster en el arco [°]
NBOOT = 2000
rng = np.random.default_rng(20260928)


def leer(obj):
    nom = pd.read_csv(os.path.join(CORRIDA, obj, f"verifica_ganancia_brewster_{STAMP}.csv"))
    mc = pd.read_csv(os.path.join(CORRIDA, obj, f"montecarlo_resumen_{STAMP}.csv"))
    raw = sio.loadmat(os.path.join(CORRIDA, obj, f"montecarlo_{STAMP}.mat"))
    PrM = raw["K_mono"] * raw["sigma_mono_m2"]      # [Nmc x Ne] W
    PrB = raw["K_esp"] * raw["sigma_esp_m2"]
    return nom, mc, PrM, PrB


def dif_banda(PrM, PrB, cols):
    """Diferencia energética B-M [dB] sumando las columnas dadas, con IC 95 % pareado."""
    d = 10 * np.log10(PrB[:, cols].sum() / PrM[:, cols].sum())
    sM = PrM[:, cols].sum(axis=1)
    sB = PrB[:, cols].sum(axis=1)
    ib = rng.integers(0, len(sM), size=(NBOOT, len(sM)))
    r = 10 * np.log10(sB[ib].sum(axis=1) / sM[ib].sum(axis=1))
    return d, np.percentile(r, 2.5), np.percentile(r, 97.5)


def f(x, dec=2, signo=True):
    s = f"{x:+.{dec}f}" if signo else f"{x:.{dec}f}"
    return f"${s}$" if s.startswith("-") or s.startswith("+") else s


datos = {o: leer(o) for o in OBJETOS}

# ─── figura ───────────────────────────────────────────────────────────────────
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})
fig, axs = plt.subplots(3, 1, figsize=(7.2, 9.0), sharex=True, constrained_layout=True)
dx = {o: (i - 2) * 0.45 for i, o in enumerate(OBJETOS)}
for o in OBJETOS:
    nom, mc, _, _ = datos[o]
    k = nom.eps_deg < 90
    e = nom.eps_deg[k].to_numpy()
    kw = dict(color=COLOR[o], marker=MARCA[o], ms=4.5, lw=1.2, label=ETIQ[o])
    axs[0].plot(e, (nom.Pr_bi_dBm - nom.Pr_mono_dBm)[k], **kw)
    d = mc.dif_media_dB[k].to_numpy()
    lo, hi = mc.IC95_dif_inf[k].to_numpy(), mc.IC95_dif_sup[k].to_numpy()
    axs[1].errorbar(e + dx[o], d, yerr=[d - lo, hi - d], capsize=2, elinewidth=0.9, **kw)
    axs[2].plot(e, mc.P_esp_mayor[k], **kw)
for ax in axs:
    ax.axvline(EPS_B, color="0.4", ls="--", lw=0.9)
axs[0].axhline(0, color="k", lw=0.8)
axs[1].axhline(0, color="k", lw=0.8)
axs[2].axhline(0.5, color="k", lw=0.8)
axs[0].text(EPS_B + 0.6, 0.97, r"$\epsilon_B$", transform=axs[0].get_xaxis_transform(), va="top")
axs[0].set_ylabel(r"$\Delta_k = P_r^B - P_r^M$ [dB]")
axs[1].set_ylabel(r"$\Delta_k$ (media, IC 95 %) [dB]")
axs[2].set_ylabel(r"$\Pr\{P_r^B > P_r^M\}$")
axs[2].set_xlabel(r"Elevación $\epsilon$ [°]")
axs[2].set_ylim(0, 1)
axs[0].set_title("(a) Orientación nominal", loc="left")
axs[1].set_title(r"(b) Monte Carlo: diferencia de potencia media", loc="left")
axs[2].set_title("(c) Monte Carlo: probabilidad de que el especular reciba más", loc="left")
axs[2].set_xticks(range(20, 90, 5))
for ax in axs[1:]:
    for a, b in BANDAS[1:]:
        ax.axvline(a - 2.5, color="0.75", lw=0.6, ls=":")
axs[0].set_ylim(-13, 45)
axs[0].legend(ncol=5, fontsize=8, loc="upper right")
os.makedirs(FIGDIR, exist_ok=True)
fig.savefig(os.path.join(FIGDIR, "fig_resumen_por_angulo.png"), dpi=200)

# ─── diferencias por banda ────────────────────────────────────────────────────
filas = []
for o in OBJETOS:
    nom, mc, PrM, PrB = datos[o]
    eps = mc.eps_deg.to_numpy()
    for a, b in BANDAS + [(20, 85)]:
        cols = np.where((eps >= a) & (eps <= b))[0]
        dn = 10 * np.log10(nom.Pr_bi_dBm.pipe(lambda s: 10 ** (s / 10))[cols].sum()
                           / nom.Pr_mono_dBm.pipe(lambda s: 10 ** (s / 10))[cols].sum())
        d, lo, hi = dif_banda(PrM, PrB, cols)
        pg = (PrB[:, cols] > PrM[:, cols]).mean()
        filas.append(dict(objeto=o, eps_min=a, eps_max=b, dif_nominal_dB=dn, dif_mc_dB=d,
                          ic95_inf=lo, ic95_sup=hi, P_esp_mayor=pg))
bandas = pd.DataFrame(filas)
bandas.to_csv(os.path.join(CORRIDA, f"resumen_por_banda_{STAMP}.csv"), index=False)
print(bandas.round(3).to_string(index=False))

# ─── filas LaTeX ──────────────────────────────────────────────────────────────
print("\n%% Tabla: Delta_k nominal por posición y objeto")
nomAll = {o: datos[o][0] for o in OBJETOS}
for i, e in enumerate(nomAll["cubo"].eps_deg):
    if e == 90:
        continue
    vals = [nomAll[o].Pr_bi_dBm[i] - nomAll[o].Pr_mono_dBm[i] for o in OBJETOS]
    print(f"{e:.0f} & " + " & ".join(f(v) for v in vals) + r" \\")

print("\n%% Tabla: bandas (nominal | MC con IC | Pr)")
for o in OBJETOS:
    s = bandas[bandas.objeto == o]
    celdas = []
    for _, r in s.iterrows():
        sig = r"\sig" if (r.ic95_inf > 0 or r.ic95_sup < 0) else ""
        celdas.append(f"{f(r.dif_mc_dB)}{sig}")
    print(f"{ETIQ[o]} & " + " & ".join(celdas) + r" \\")

for o in OBJETOS:
    nom, mc, _, _ = datos[o]
    print(f"\n%% Apéndice {o}: nominal")
    for i, r in nom.iterrows():
        print(f"{r.eps_deg:.0f} & {f(r.sigma_mono_dBsm)} & {f(r.sigma_bi_dBsm)} & "
              f"{f(r.Pr_mono_dBm)} & {f(r.Pr_bi_dBm)} & {f(r.Pr_bi_dBm - r.Pr_mono_dBm)} & "
              f"{f(r.Pr_bi_dBm - r.Pr_superficie_dBm)} \\\\")
    print(f"%% Apéndice {o}: Monte Carlo")
    for i, r in mc.iterrows():
        if r.eps_deg == 90:
            continue
        sig = r"\sig" if (r.IC95_dif_inf > 0 or r.IC95_dif_sup < 0) else ""
        print(f"{r.eps_deg:.0f} & {f(r.Pr_mono_media_dBm)} & {f(r.Pr_esp_media_dBm)} & "
              f"{f(r.dif_media_dB)}{sig} & $[{r.IC95_dif_inf:+.2f},\\,{r.IC95_dif_sup:+.2f}]$ & "
              f"{f(r.Pr_mono_mediana_dBm)} & {f(r.Pr_esp_mediana_dBm)} & "
              f"{r.P_esp_mayor:.2f} & {r.aciertos_mono:.0f} / {r.aciertos_esp:.0f} \\\\")
