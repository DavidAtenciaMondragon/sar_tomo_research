# SAR Tomography Research

Research project on 3D image formation via **bistatic helical SAR with two-medium propagation** (air and soil with Snell refraction). The goal is to analytically derive and numerically validate 3D spatial resolution formulas for subsurface target detection using UAV-borne radar systems, and to optimize the receiver's flight path to jointly maximize received power (via the Brewster angle) and 3D imaging resolution.

## Overview

The system models two radar platforms (transmitter and receiver) flying helical cone trajectories around a buried target volume. The electromagnetic signal crosses a flat air-soil interface (or a stack of soil layers) governed by Snell's law. A backprojection algorithm reconstructs the 3D image from the bistatic phase history, and a separate optimization pipeline computes, for a fixed transmitter trajectory, the receiver position that best balances received power and spatial resolution at every instant of the flight.

**Key results:** analytical resolution formulas validated against simulation with < 12% error across all tested configurations; receiver flight-path optimization exploiting the Brewster angle, generalized to layered (multislab) soils.

## Repository Structure

```
sar_tomo_research/
├── simulations/
│   ├── run_snell_pipeline.m              # Bistatic 2-medium image formation (GS→PROC)
│   ├── run_snell_multislab_pipeline.m    # Same, generalized to N soil layers
│   ├── run_snell_volumetrico_pipeline.m  # Volumetric target variant
│   ├── run_plano_de_voo.m                # Receiver flight-path optimization (2-medium)
│   ├── run_plano_de_voo_multislab.m      # Receiver flight-path optimization (N layers)
│   ├── run_barrido_alpha_grid_multislab.m# Sensitivity sweep (power/resolution weight, grid size)
│   ├── gs/                               # Signal Generator (GS)
│   ├── proc/                             # Image Processor (backprojection)
│   ├── flightpath/                       # Flight-path optimizer (cost function, search, plots)
│   ├── common/                           # Shared functions (Snell solvers, multislab TM)
│   ├── tools/                            # I/O utilities (JSON, binary)
│   ├── parametros/                       # System configuration (JSON)
│   └── io/                               # Simulation outputs (CSV, figures)
├── doc/                                  # LaTeX/Markdown technical documentation and PhD dissertation
├── models/                               # Theoretical derivations (Markdown)
├── hypotheses/                           # Evolution of resolution hypotheses
├── disertaciones/                        # Dissertation structure guidelines
├── papers/                               # Reference literature (PDF)
├── utils/                                # Figure-generation scripts (Python)
└── notes/                                # Paper summaries and literature review
```

## Physical Model

- **Geometry:** bistatic SAR, TX and RX on separate helical cone trajectories
- **Interface:** flat air-soil boundary at z = 0, refractive index n₂ = √εᵣ
- **Optical path:** R_OP = d₁^TX + n₂·d₂^TX + n₂·d₂^RX + d₁^RX
- **Phase history:** S(f, t) = A · exp(−j · 2πf/c · R_OP)
- **Image formation:** coherent backprojection over a 3D processing grid

## Simulation Pipeline

```
JSON parameters → [GS] phase history (.mat) → [PROC] 3D image → −3dB resolution (CSV)
```

1. Edit parameters in `simulations/parametros/` (frequency, trajectory geometry, target position)
2. Run the pipeline from MATLAB:
   ```matlab
   run_snell_pipeline.m
   ```
3. Results are written to `simulations/io/snell/`

Two parameter sets are available:
- `_espiral` — X-band (10 GHz), reference configuration
- `_espiral_EMIRADOS` — P-band (425 MHz), realistic penetration depth for buried targets

## Receiver Flight-Path Optimization

For a transmitter flying a fixed conical-helix trajectory, the receiver's horizontal position is
a free degree of freedom at every instant. `run_plano_de_voo.m` (two-medium soil) and
`run_plano_de_voo_multislab.m` (N-layer soil) search, for each transmitter position, the receiver
position that minimizes a log-normalized cost combining received power (maximized near the
Brewster angle) and expected 3D resolution, using multi-start constrained optimization:

```
J = -(1 - α)·log10(Pr) + α·log10(δ_xy · δ_z)
```

1. Edit parameters in `simulations/parametros/*_plano_voo.json` (helix geometry, soil layers, α weight)
2. Run from MATLAB: `run_plano_de_voo.m` or `run_plano_de_voo_multislab.m`
3. Results (optimized receiver trajectory, power, resolution, cost) are written to
   `simulations/io/plan_vuelo/` or `simulations/io/plan_vuelo_multislab/`, in georeferenced
   tabular form ready for import into a GIS (e.g. QGIS) for flight planning and field validation.

Full theoretical background, reasoning, configuration parameters and benefits in
[`doc/proposta_voo_bistatico.tex`](doc/proposta_voo_bistatico.tex).

## Validated Resolution Formulas

| Resolution | Formula | Error vs. simulation |
|---|---|---|
| Vertical δ_z | c / (2 W_z), where W_z = n₂B cosθ_t + [tomographic term] | < 1% |
| Horizontal δ_xy | 0.60λ / (π sinψ₀ \|cos(ΔΦ/2)\|) | < 12% |

Full derivation in [`models/derivacion_modelo_resolucion.md`](models/derivacion_modelo_resolucion.md).

## Requirements

- MATLAB R2019b or later (`fmincon`/Optimization Toolbox required for the flight-path optimizer)
- Python 3 (optional, for figure generation via `utils/generate_figures*.py`)
- A LaTeX distribution (MiKTeX or TeX Live) to compile the documents under `doc/` and `hypotheses/`

## Documentation

| File | Contents |
|---|---|
| `models/modelo_fisico_base.md` | Complete physical model specification |
| `models/derivacion_modelo_resolucion.md` | Step-by-step resolution derivation (88 equations) |
| `models/analisis_simulacion.md` | Simulation validation and measured results |
| `models/multiple_slab_tm.md` | ABCD transfer-matrix formalism for N-layer TM transmittance |
| `notes/modelo_conceptual.md` | Literature synthesis and comparison with prior work |
| `hypotheses/hipotesis_v5.md` | Latest validated resolution hypothesis, including off-axis targets |
| `doc/proposta_voo_bistatico.tex` | Receiver flight-path optimization: theory, reasoning, parameters and benefits |
| `doc/explicacion_plan_de_vuelo.tex` | Full derivation of the two-medium flight-path optimizer |
| `doc/plan_de_vuelo_multislab.md` | Multislab extension of the flight-path optimizer |
| `doc/dissertacao_tomo_sar.tex` | PhD dissertation draft (Portuguese) |
