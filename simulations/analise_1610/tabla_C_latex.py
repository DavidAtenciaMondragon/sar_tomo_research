"""Tablas LaTeX coloreadas del Exp. C a partir de expC/expC_scr_imagen.csv:
  doc/relatorio_1610_tablaC_scr.tex   SCR en imagen, S1, s = 2 cm: TM/TE x F1/F2 x incidencia x canal
                                      (verde = mejor canal de la fila; rojo = SCR < 0 dB)
  doc/relatorio_1610_tablaC_err.tex   SCR frente al error de incidencia en theta_pB (TM, S1 y S0)
  doc/relatorio_1610_tablaC_f1f2.tex  F1 frente a F2: diferencia media de SCR por canal y victorias
Uso: python tabla_C_latex.py
Colores definidos en relatorio_1610.tex: mejor (verde), peor (rojo).
"""
import os

import numpy as np
import pandas as pd

from bisar import config

DOC = os.path.dirname(os.path.dirname(config.FIGDOC))
d = pd.read_csv(os.path.join(config.outdir("expC"), "expC_scr_imagen.csv"))
d["ch"] = np.where(d.canal == "mono", "mono", d.canal + d.dphi.fillna(0).astype(int).astype(str))
CH = ["mono", "bi90", "bi135", "bi180", "fus90", "fus135"]
HEAD = ["Mono", "Bi 90°", "Bi 135°", "Ref1 180°", "Mono+bi 90°", "Mono+bi 135°"]
THB = d[(d.suelo == "S1") & (d.tipo == "pB") & (d.err_deg == 0)].theta_c.iloc[0]


def cell(v, best, fmt="{:.1f}"):
    s = fmt.format(v).replace("-", "$-$")
    if best:
        return r"\cellcolor{mejor}\textbf{" + s + "}"
    if v < 0:
        return r"\cellcolor{peor}" + s
    return s


def theta_lab(th):
    return r"$\theta_{pB}$=63.45°" if abs(th - THB) < 0.01 else (f"{th:.1f}°" if th % 1 else f"{th:.0f}°")


# ------------------------------------------------------------- tabla principal
base = d[(d.suelo == "S1") & (d.sigma_h == 0.02) & (d.tipo.isin(["fuera", "pB"])) &
         ((d.err_deg == 0) | d.err_deg.isna())]
T = base.pivot_table(index=["pol", "familia", "theta_c"], columns="ch", values="SCR_img_dB")[CH]
L = [r"\begin{table}[H]", r"\centering\small",
     r"\caption{SCR en imagen [dB] en el vóxel del blanco: S1, $s=2$\,cm, superficie congelada (3 realizaciones). "
     r"\colorbox{mejor}{\textbf{verde}}: mejor canal de la fila; \colorbox{peor}{rojo}: SCR $<0$ (la superficie "
     r"domina). $\leftarrow$: fila en el pseudo-Brewster. CSV: \texttt{expC/expC\_scr\_imagen.csv}.}",
     r"\label{tab:C_scr}", r"\begin{tabular}{@{}llrrrrrr@{}}", r"\toprule",
     r"Pol. & Familia, $\theta_c$ & " + " & ".join(HEAD) + r" \\", r"\midrule"]
for pol in ["TM", "TE"]:
    for fam in ["F1", "F2"]:
        sub = T.loc[(pol, fam)]
        for i, (th, r) in enumerate(sub.iterrows()):
            ib = int(np.argmax(r.values))
            vals = " & ".join(cell(v, j == ib) for j, v in enumerate(r.values))
            lab = (pol if (fam == "F1" and i == 0) else "")
            fl = f"{fam}, {theta_lab(th)}"
            if abs(th - THB) < 0.01:
                fl = r"\textbf{" + fl + r"} $\leftarrow$"
            L.append(f"{lab} & {fl} & {vals} \\\\")
        if not (pol == "TE" and fam == "F2"):
            L.append(r"\cmidrule(l){2-8}" if fam == "F1" else r"\midrule")
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
open(os.path.join(DOC, "relatorio_1610_tablaC_scr.tex"), "w", encoding="utf-8").write("\n".join(L) + "\n")

# ------------------------------------------------------------- F1 frente a F2
rows = []
for ch, hd in zip(CH, HEAD):
    a = T.xs("F1", level="familia")[ch]
    b = T.xs("F2", level="familia")[ch]
    diff = (a - b).values
    rows.append((hd, diff.mean(), diff.min(), diff.max(), int((diff > 0).sum()), len(diff)))
L = [r"\begin{table}[H]", r"\centering\small",
     r"\caption{F1 (hélice) frente a F2 (cono): $\mathrm{SCR}_{F1}-\mathrm{SCR}_{F2}$ en la misma incidencia "
     r"central, polarización y canal (10 pares por canal: 5 incidencias $\times$ TM/TE, Tabla~\ref{tab:C_scr}). "
     r"\colorbox{mejor}{verde}: F1 gana.}",
     r"\label{tab:C_f1f2}", r"\begin{tabular}{@{}lrrrc@{}}", r"\toprule",
     r"Canal & media [dB] & mín. [dB] & máx. [dB] & F1 gana \\", r"\midrule"]
for hd, m, mn, mx, w, n in rows:
    c = r"\cellcolor{mejor}" if m > 0 else r"\cellcolor{peor}"
    L.append(f"{hd} & {c}$+${m:.1f} & {mn:.1f} & {mx:.1f} & {w}/{n} \\\\".replace("$+$-", "$-$").replace(" -", " $-$"))
tot_w = sum(r[4] for r in rows); tot_n = sum(r[5] for r in rows)
L += [r"\midrule", f"Todos & $+${np.mean([r[1] for r in rows]):.1f} & & & {tot_w}/{tot_n} \\\\",
      r"\bottomrule", r"\end{tabular}", r"\end{table}"]
open(os.path.join(DOC, "relatorio_1610_tablaC_f1f2.tex"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("F1 vs F2:", [(r[0], round(r[1], 1), f"{r[4]}/{r[5]}") for r in rows], f"total {tot_w}/{tot_n}")

# ------------------------------------------------------------- robustez frente al error angular
E = d[(d.tipo == "pB") & (d.pol == "TM") & (d.sigma_h == 0.02)]
CHE = ["mono", "bi135", "bi180", "fus135"]
HE = ["Mono", "Bi 135°", "Ref1 180°", "Mono+bi 135°"]
L = [r"\begin{table}[H]", r"\centering\small",
     r"\caption{Robustez del SCR en imagen [dB] frente al error de incidencia $\delta$ alrededor de $\theta_{pB}$ "
     r"(TM, $s=2$\,cm). Pares $a$\,/\,$b$: $\delta<0$\,/\,$\delta>0$. \colorbox{mejor}{verde}: mejor canal de la "
     r"fila; \colorbox{peor}{rojo}: SCR $<0$.}",
     r"\label{tab:C_err_scr}", r"\begin{tabular}{@{}llrrrr@{}}", r"\toprule",
     r"Suelo, familia & $|\delta|$ & " + " & ".join(HE) + r" \\", r"\midrule"]
for sk, fam in [("S1", "F1"), ("S1", "F2"), ("S0", "F2")]:
    S = E[(E.suelo == sk) & (E.familia == fam)]
    errs = sorted(set(np.abs(S.err_deg.dropna())))
    for i, e in enumerate(errs):
        vals, means = [], []
        for ch in CHE:
            neg = S[(S.err_deg == -e) & (S.ch == ch)].SCR_img_dB
            pos = S[(S.err_deg == e) & (S.ch == ch)].SCR_img_dB
            vv = [neg.iloc[0], pos.iloc[0]] if e > 0 else [pos.iloc[0]]
            vals.append(vv); means.append(np.mean(vv))
        ib = int(np.argmax(means))
        cells = []
        for j, vv in enumerate(vals):
            s = " / ".join(f"{v:.1f}" for v in vv).replace("-", "$-$")
            if j == ib:
                s = r"\cellcolor{mejor}\textbf{" + s + "}"
            elif max(vv) < 0:
                s = r"\cellcolor{peor}" + s
            cells.append(s)
        lab = f"{sk}, {fam}" if i == 0 else ""
        L.append(f"{lab} & {e:g}° & " + " & ".join(cells) + r" \\")
    if not (sk == "S0"):
        L.append(r"\midrule")
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
open(os.path.join(DOC, "relatorio_1610_tablaC_err.tex"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("->", DOC)
