"""Experimento F2 - Calibracion de fase entre canales y fusion coherente vs incoherente.

Hallazgo del Exp. F: un offset de sincronismo constante dt equivale (a la portadora) a una fase
constante phi0 = 2 pi fc dt entre el canal bi y el mono (1 ns -> 153 deg) y puede casi anular la
fusion coherente.  Aqui, con la trayectoria propuesta del Exp. F:
  (a) barrido de phi0 = 0..180 deg en la fusion coherente (I_m + e^{j phi0} I_b)/2;
  (b) calibracion: phi0_est = arg( sum_{lobulo ppal.} I_m^* I_b ) sobre el lobulo principal del
      dispersor mas brillante (equivale a alinear la fase en un reflector fuerte o con la senal
      directa Tx->Rx), y se re-fusiona.  Nota: estimarla sobre TODA la imagen falla (~180 deg),
      porque dominan los anillos de lobulos secundarios en contrafase;
  (c) fusion incoherente (|I_m|^2 + |I_b|^2)/2 ("sumar potencias"), insensible a la fase.
Salidas: expF/expF2_fase.csv, figura F2_fase.  Uso: python exp_F2_calibracion_fase.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bisar import config, geometry, imaging, metrics, plotting as P, runner, scene

p = config.load("expF")
out = config.outdir("expF")
sysm = scene.System(p)
tgt = np.array(p["blanco"])
t, Ptx = geometry.tx_trajectory(p)
T = t[-1]
best = pd.read_csv(f"{out}/expF_candidatos.csv").sort_values("J", ascending=False).iloc[0]
Prx = geometry.rx_from_tx(Ptx, best.dphi0 + best.A * (2 * t / T - 1), best.kappa)
rg = p["rejilla_psf"]
g = imaging.Grid.around(tgt, rg["hx"], rg["dx"], rg["z_lo"], rg["z_hi"], rg["dz"])
s_m, _ = scene.target_echo(sysm, Ptx, Ptx, tgt)
s_b, _ = scene.target_echo(sysm, Ptx, Prx, tgt)
I_m = runner.bp(sysm, Ptx, Ptx, s_m, g.vox, gate=True)
I_b = runner.bp(sysm, Ptx, Prx, s_b, g.vox, gate=True)
rows = [{"caso": "mono", "phi0_deg": np.nan, **metrics.psf_metrics(g, I_m, tgt)}]
ML = np.abs(I_m) >= 0.5 * np.abs(I_m).max()          # lobulo principal del dispersor mas brillante


def phase_est(Ib):
    """Fase residual estimada sobre el lobulo principal del blanco mas brillante."""
    return np.angle(np.vdot(I_m[ML], Ib[ML]))
for ph in [0, 30, 60, 90, 120, 150, 180]:
    Ib = I_b * np.exp(1j * np.radians(ph))
    rows.append({"caso": "coherente_sin_calibrar", "phi0_deg": ph, **metrics.psf_metrics(g, 0.5 * (I_m + Ib), tgt)})
    est = phase_est(Ib)
    rows.append({"caso": "coherente_calibrada", "phi0_deg": ph, "phi0_est_deg": np.degrees(est),
                 **metrics.psf_metrics(g, 0.5 * (I_m + Ib * np.exp(-1j * est)), tgt)})
# incoherente: la metrica de PSF se calcula sobre la raiz de la intensidad media
inc = np.sqrt(0.5 * (np.abs(I_m) ** 2 + np.abs(I_b) ** 2))
rows.append({"caso": "incoherente", "phi0_deg": np.nan, **metrics.psf_metrics(g, inc, tgt)})
# ademas: offset de sincronismo de 1 ns calibrado (fase + desplazamiento en rango de 0.3 m)
s_b1 = s_b * np.exp(-2j * np.pi * np.outer(np.full(len(Ptx), 1e-9), sysm.freqs))
I_b1 = runner.bp(sysm, Ptx, Prx, s_b1, g.vox, gate=True)
est = phase_est(I_b1)
rows.append({"caso": "sinc_1ns_sin_calibrar", "phi0_deg": np.nan, **metrics.psf_metrics(g, 0.5 * (I_m + I_b1), tgt)})
rows.append({"caso": "sinc_1ns_calibrada_fase", "phi0_deg": np.nan, "phi0_est_deg": np.degrees(est),
             **metrics.psf_metrics(g, 0.5 * (I_m + I_b1 * np.exp(-1j * est)), tgt)})
d = pd.DataFrame(rows)
d.to_csv(f"{out}/expF2_fase.csv", index=False)
print(d[["caso", "phi0_deg", "phi0_est_deg", "PSLR_dB", "ISLR_dB", "vol3dB_m3", "fwhm_x", "contrast_dB",
         "loc_err_m", "peak_amp"]].round(4).to_string())

fig, axs = plt.subplots(1, 3, figsize=(12, 3.3))
for caso, col, lab in [("coherente_sin_calibrar", P.C["bi"], "coherente sin calibrar"),
                       ("coherente_calibrada", P.C["fus"], "coherente calibrada")]:
    x = d[d.caso == caso]
    for a, met in zip(axs, ["PSLR_dB", "ISLR_dB", "contrast_dB"]):
        a.plot(x.phi0_deg, x[met], marker="o", color=col, label=lab)
for a, met, yl in zip(axs, ["PSLR_dB", "ISLR_dB", "contrast_dB"], ["PSLR [dB]", "ISLR [dB]", "contraste [dB]"]):
    a.axhline(d[d.caso == "mono"][met].values[0], color=P.C["mono"], label="mono")
    a.axhline(d[d.caso == "incoherente"][met].values[0], color=P.SERIES[6], ls="--", label="incoherente (suma de potencias)")
    a.set_xlabel("fase residual entre canales $\\phi_0$ [°]"); a.set_ylabel(yl); a.set_xticks(range(0, 181, 30))
axs[0].legend(fontsize=7)
fig.suptitle("Exp. F2: la fusion coherente exige calibrar UNA fase entre canales; sumar potencias no basta", y=1.03)
P.save(fig, "F2_fase", out)
