# analise_1610: battery of experiments for the 03/10/2026 guidance

Answers the guidance in `Orientacoes_David_SAR_biestatico_03-10-2026_final.docx.pdf`.
The resulting report is `doc/relatorio_1610.tex`.

## Requirements
Python 3.11+ with `numpy`, `scipy`, `pandas`, `matplotlib`, `numba` (the Anaconda environment already has them).
The validation scripts against the MATLAB pipeline need MATLAB R2021a (only `tests/valida_bp_matlab.m`).

## Reproduction
```
cd simulations/analise_1610
python tests/t_rcs.py            # PO vs Table 2 of the 28/09 report and the specular plate
python tests/t_bp_matlab.py      # BP vs MATLAB pipeline GS/PROC (10 GHz, mono0/offset0/offset180)
python tests/t_superficie.py     # Kirchhoff facet sum vs image theory; diffuse SPM vs sigma0 integral
python run_all.py                # A..F plus complements, in order (about 4-5 h on 16 cores)
```
Each experiment can also be run on its own (`python exp_B_barrido_dphi.py`, ...).
Complements added after reviewing the first results (also run by `run_all.py`):

| Script | Why it was added |
|---|---|
| `exp_C_figuras_enlace.py` | Fig. 7 with S0/S1 robustness to angular and attitude error (analytic) |
| `exp_D2_rejilla_densa.py` | Horizontal separability with spacing near the resolution (the 1 m grid does not discriminate) |
| `exp_D3_rx_alineado.py` | Nominal objects with an aligned-polarization Rx (dual-pol Rx) |
| `BISAR_EXP=expE2 python exp_E_montecarlo.py` | SO(3) Monte Carlo with an aligned-polarization Rx |
| `exp_F2_calibracion_fase.py` | Inter-channel phase calibration and coherent vs incoherent fusion |
| `exp_F_fig_robustez.py` | Robustness figure (log axis) from `expF_robustez.csv` |

Processing note: C and F use Hann spectral weighting (`"ventana_espectral": "hann"`).
Without it, the range sidelobes of the nadir reflection accumulate coherently on the axis of a
constant-altitude orbit. The rectangular-weighting C results are kept in `io/analise_1610/expC_rect`.

## Structure
| Path | Contents |
|---|---|
| `params/common.json` | Parameters shared by ALL experiments (Tx trajectory, 425 MHz, 50 MHz, soils, roughness) |
| `params/expX.json` | Single parameter file per experiment (inherits from `common.json`) |
| `bisar/` | Engine: geometry, EM (complex Fresnel, pseudo-Brewster, SPM), bistatic PO, echoes, BP, metrics |
| `exp_A_referencias.py` | Ref2, Ref1 and Ref1+Ref2 fusion |
| `exp_B_barrido_dphi.py` | Sweep dphi = 0:15:180 (mono, bi, 2N fusion and equal-budget N fusion) |
| `exp_B_fig_geometria.py` | Exp. B geometry figure: Tx and Rx fly the same helix, offset by dphi in azimuth |
| `t_resolucion_analitica.py` | Closed-form delta_x, delta_y, delta_z (`bisar/resolucion.py`, report Appendix A) vs. the image; writes `expB/expB_resolucion_analitica.csv` |
| `tabla_resolucion_latex.py` | Builds `doc/relatorio_1610_tablaResol.tex` from that CSV |
| `exp_C_brewster.py` | Pseudo-Brewster: reflectance, SPM sigma0, SCR (link and image), angular/attitude error, image cost |
| `exp_C_fig_geometria.py` | Exp. C geometry figure: F1/F2 families, incidence per pulse, bistatic specular point and its incidence vs dphi |
| `exp_C_fig_trayectorias3d.py` | Exp. C full 3D trajectories: base helix, F1 and F2 with the constant-incidence cones |
| `tabla_C_latex.py` | Builds the colour-coded Exp. C tables (`doc/relatorio_1610_tablaC_*.tex`) from `expC/expC_scr_imagen.csv` |
| `exp_D_biblioteca.py` | P0/P1/P2 target library, 64-point grid |
| `exp_E_montecarlo.py` | Monte Carlo in SO(3), one orientation per acquisition |
| `exp_F_trayectoria.py` | Pareto front, weighted J, final trajectory and robustness |
| `../io/analise_1610/expX/` | CSV with ALL metrics + `params_resueltos_expX.json` (parameters actually used) |
| `../../doc/figures/relatorio_1610/` | Figures (PDF vector + PNG) |

## Conventions (guidance §5.1 and §11)
* Complex signals in all channels; mono and bi of the **same pulse** (same Tx trajectory, target and physics).
* Each image is the coherent **mean** over its samples. The fusion is `fus = (I_mono + I_bi)/2` (2N samples).
  The `fusN` fusion uses mono on odd pulses and bi on even pulses, so it has the same N as mono.
  The SNR is computed by propagating the thermal noise, so the sample-count gain stays explicit (+3 dB) and separate.
* All images share the same dynamic range (30 dB) and the same reference (monostatic peak).
* Monte Carlo draws one orientation per acquisition, the same in mono and bi.
* The target received power is **not** used as a criterion or a variable anywhere.
* The Rx flies at the same height as the Tx, `z_Rx(t) = z_Tx(t)`. Its radius is `kappa * rho_Tx(t)` and its azimuth is `phi_Tx(t) + dphi(t)`.
