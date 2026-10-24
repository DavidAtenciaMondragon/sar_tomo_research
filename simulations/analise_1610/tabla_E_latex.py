"""Genera doc/relatorio_1610_tablaE.tex a partir de expE/expE_resumen.csv y expE2/expE_resumen.csv."""
import os
import pandas as pd

from bisar import config

IO = config.IO
DOC = os.path.dirname(config.FIGDOC.rstrip("/\\").rsplit(os.sep, 1)[0])
OUT = os.path.join(os.path.dirname(os.path.dirname(config.FIGDOC)), "relatorio_1610_tablaE.tex")
MET = [("SNR_dB", "$\\Delta$SNR", "{:+.1f}"), ("contrast_dB", "$\\Delta$Contr.", "{:+.1f}"),
       ("PSLR_dB", "$\\Delta$PSLR", "{:+.2f}"), ("vol3dB_m3", "$\\Delta V$", "{:+.1f}"),
       ("loc_err_h_m", "$\\Delta$loc", "{:+.3f}")]


def block(path, caption, label):
    d = pd.read_csv(path)
    pw = d[d.metrica == "potencia_media_bi_vs_mono_dB"]
    d = d[d.canal == "fusN"]
    lines = ["\\begin{table}[H]", "\\centering\\small", f"\\caption{{{caption}}}", f"\\label{{{label}}}",
             "\\resizebox{\\textwidth}{!}{\\begin{tabular}{@{}lrr|" + "r" * (2 * len(MET)) + "@{}}", "\\toprule",
             "Objeto & $\\Delta\\phi$ & $\\bar P_{bi}/\\bar P_{mono}$ & " +
             " & ".join(f"\\multicolumn{{2}}{{c}}{{{lab}}}" for _, lab, _ in MET) + " \\\\",
             " & & [dB] & " + " & ".join("med. & P$_{mej}$" for _ in MET) + " \\\\", "\\midrule"]
    for obj in d.objeto.unique():
        for dphi in sorted(d.dphi.unique()):
            row = [obj.replace("_", " "), f"{dphi}°"]
            p = pw[(pw.objeto == obj) & (pw.dphi == dphi)]
            row.append(f"{p.media_dif.values[0]:+.1f}" if len(p) else "---")
            for m, _, fmt in MET:
                x = d[(d.objeto == obj) & (d.dphi == dphi) & (d.metrica == m)]
                if len(x):
                    row += [fmt.format(x.mediana_dif.values[0]), f"{100 * x.P_mejora.values[0]:.0f}\\,\\%"]
                else:
                    row += ["---", "---"]
            lines.append(" & ".join(row) + " \\\\")
        lines.append("\\midrule")
    lines[-1] = "\\bottomrule"
    lines += ["\\end{tabular}}", "\\end{table}", ""]
    return "\n".join(lines)


txt = "% generado por simulations/analise_1610/tabla_E_latex.py\n"
cap = ("Monte Carlo SO(3) por adquisici\\'on: mediana de (fusi\\'on N $-$ mono) y probabilidad de mejora "
       "P$_{{mej}}$ (Rx {rx}). $\\Delta V$ en dB; loc.\\ en m. N = {n}.")
txt += block(os.path.join(IO, "expE", "expE_resumen.csv"), cap.format(rx="TM co-polar", n="1000 (cubo), 200"), "tab:E_TM")
p2 = os.path.join(IO, "expE2", "expE_resumen.csv")
if os.path.exists(p2):
    txt += block(p2, cap.format(rx="alineado, doble polarizaci\\'on", n="200"), "tab:E_AL")
open(OUT, "w", encoding="utf-8").write(txt)
print("escrito", OUT)
