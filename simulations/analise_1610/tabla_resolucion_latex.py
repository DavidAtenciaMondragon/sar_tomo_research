"""Tabla LaTeX de verificacion de las formas cerradas (Apendice A) a partir de
io/analise_1610/expB/expB_resolucion_analitica.csv -> doc/relatorio_1610_tablaResol.tex
Uso: python tabla_resolucion_latex.py   (despues de t_resolucion_analitica.py)
"""
import os
import pandas as pd

from bisar import config

d = pd.read_csv(os.path.join(config.outdir("expB"), "expB_resolucion_analitica.csv"))
SEL = {"centro": ["B dphi=0", "B dphi=45", "B dphi=90", "B dphi=135", "B dphi=165",
                  "radio 1.3x, dphi=0", "Rx 60 m mas bajo, dphi=0", "Rx casi nadir (0.2x), dphi=0"],
       "descentrado": ["B dphi=0", "B dphi=45", "B dphi=90", "B dphi=135", "B dphi=165",
                       "radio 1.3x, dphi=0", "Rx 60 m mas bajo, dphi=0", "Rx casi nadir (0.2x), dphi=0"]}
NOM = {"B dphi=0": r"$\Delta\phi=0^{\circ}$ (mono)", "B dphi=45": r"$\Delta\phi=45^{\circ}$",
       "B dphi=90": r"$\Delta\phi=90^{\circ}$", "B dphi=135": r"$\Delta\phi=135^{\circ}$",
       "B dphi=165": r"$\Delta\phi=165^{\circ}$", "radio 1.3x, dphi=0": r"$\Delta\phi=0$, $\kappa=1.3$",
       "Rx 60 m mas bajo, dphi=0": r"$\Delta\phi=0$, $z_{Rx}=z_{Tx}-60$\,m",
       "Rx casi nadir (0.2x), dphi=0": r"$\Delta\phi=0$, $\kappa=0.2$"}


def pct(a, b):
    return f"{100 * (a / b - 1):+.1f}"


L = [r"\begin{table}[H]", r"\centering\footnotesize\setlength{\tabcolsep}{3.6pt}",
     r"\caption{Verificación de las formas cerradas (\ref{eq:A_dxy}) y (\ref{eq:A_dz}) frente a las resoluciones "
     r"medidas en la imagen (back-projection con refracción exacta, cortes finos de 2 y 5\,mm por el pico). "
     r"$\varepsilon$: error relativo de la forma cerrada. CSV: \texttt{expB/expB\_resolucion\_analitica.csv}.}",
     r"\label{tab:A_verif}",
     r"\begin{tabular}{@{}llccccccccc@{}}", r"\toprule",
     r"Blanco & Rx & $\bar\psi_{Rx}$ & $\langle S\rangle$ & $\delta_x$ & $\delta_y$ & $\delta_{xy}^{\,\mathrm{cerr}}$ "
     r"& $\varepsilon_x/\varepsilon_y$ [\%] & $\delta_z$ & $\delta_z^{\,\mathrm{cerr}}$ & $\varepsilon_z$ [\%] \\",
     r" & & [°] & & [m] & [m] & [m] & & [m] & [m] & \\", r"\midrule"]
for b, casos in SEL.items():
    for i, c in enumerate(casos):
        r = d[(d.blanco == b) & (d.caso == c)].iloc[0]
        lab = ({"centro": r"$(0,0,-2)$", "descentrado": r"$(6,4,-2)$"}[b]) if i == 0 else ""
        L.append(f"{lab} & {NOM[c]} & {r.psi_R_media_deg:.1f} & {r.S_media:.3f} & {r.delta_x_medido:.3f} & "
                 f"{r.delta_y_medido:.3f} & {r.delta_xy_cerrada:.3f} & "
                 f"{pct(r.delta_xy_cerrada, r.delta_x_medido)}/{pct(r.delta_xy_cerrada, r.delta_y_medido)} & "
                 f"{r.delta_z_medido:.3f} & {r.delta_z_cerrada:.3f} & {pct(r.delta_z_cerrada, r.delta_z_medido)} \\\\")
    if b == "centro":
        L.append(r"\midrule")
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
dst = os.path.join(os.path.dirname(os.path.dirname(config.FIGDOC)), "relatorio_1610_tablaResol.tex")
open(dst, "w", encoding="utf-8").write("\n".join(L) + "\n")
print("->", dst)
