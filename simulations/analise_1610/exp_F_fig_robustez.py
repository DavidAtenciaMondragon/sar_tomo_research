"""Re-dibuja la figura de robustez del Exp. F desde expF_robustez.csv (eje logaritmico)."""
import pandas as pd
import matplotlib.pyplot as plt

from bisar import config, plotting as P

out = config.outdir("expF")
dr = pd.read_csv(f"{out}/expF_robustez.csv")
nom = dr[(dr.perturbacion == "nominal") & (dr.canal == "fus")].iloc[0]
mono = dr[(dr.perturbacion == "nominal") & (dr.canal == "mono")].iloc[0]
groups = [("pos_ruido_m", "ruido de posicion Rx [cm]", 100), ("pos_deriva_m", "deriva de posicion Rx [cm]", 100),
          ("fase_deg", "ruido de fase Rx [°]", 1), ("sinc_offset_ns", "offset de sincronismo [ns]", 1),
          ("sinc_deriva_ns", "deriva de sincronismo [ns]", 1)]
fig, axs = plt.subplots(1, 3, figsize=(13, 3.5))
fig.subplots_adjust(wspace=0.3)
for (gname, lab, fct), col in zip(groups, P.SERIES):
    d = dr[(dr.perturbacion == gname) & (dr.canal == "fus")].sort_values("valor")
    for a, met in zip(axs, ["PSLR_dB", "ISLR_dB", "loc_err_m"]):
        a.plot(d.valor * fct, d[met], marker="o", color=col, label=lab)
for a, met, yl in zip(axs, ["PSLR_dB", "ISLR_dB", "loc_err_m"], ["PSLR fusion [dB]", "ISLR fusion [dB]", "error de localizacion [m]"]):
    a.axhline(nom[met], color=P.C["fus"], ls="--", lw=1, label="fusion nominal")
    a.axhline(mono[met], color=P.C["mono"], lw=1.4, label="mono")
    a.set_xscale("log"); a.set_xlabel("magnitud de la perturbacion (unidades en la leyenda)"); a.set_ylabel(yl)
axs[0].legend(fontsize=6.5)
fig.suptitle("Exp. F: robustez de la fusion con la trayectoria propuesta (errores solo en el canal bi, sin calibrar)", y=1.03)
P.save(fig, "F_robustez", out)
