# SAR Tomography Research

Research project on 3D image formation via **bistatic helical SAR with two-medium propagation** (air and soil with Snell refraction). The goal is to analytically derive and numerically validate 3D spatial resolution formulas for subsurface target detection using UAV-borne radar systems.

## Overview

The system models two radar platforms (transmitter and receiver) flying synchronized helical cone trajectories around a buried target. The electromagnetic signal crosses a flat air-soil interface governed by Snell's law. A backprojection algorithm reconstructs the 3D image from the bistatic phase history.

**Key results:** analytical resolution formulas validated against simulation with < 12% error across all tested configurations.

## Repository Structure

```
sar_tomo_research/
├── simulations/
│   ├── run_snell_pipeline.m              # Main pipeline entry point
│   ├── gs/                               # Signal Generator (GS)
│   ├── proc/                             # Image Processor (backprojection)
│   ├── common/                           # Shared functions (Snell solvers)
│   ├── tools/                            # I/O utilities (JSON, binary)
│   ├── parametros/                       # System configuration (JSON)
│   └── io/                               # Simulation outputs (CSV, figures)
├── models/                               # Theoretical derivations (Markdown)
├── hypotheses/                           # Evolution of resolution hypotheses
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

## Validated Resolution Formulas

| Resolution | Formula | Error vs. simulation |
|---|---|---|
| Vertical δ_z | c / (2 W_z), where W_z = n₂B cosθ_t + [tomographic term] | < 1% |
| Horizontal δ_xy | 0.60λ / (π sinψ₀ \|cos(ΔΦ/2)\|) | < 12% |

Full derivation in [`models/derivacion_modelo_resolucion.md`](models/derivacion_modelo_resolucion.md).

## Requirements

- MATLAB R2019b or later (no additional toolboxes required)
- Python 3 (optional, for figure generation via `utils/generate_figures.py`)

## Documentation

| File | Contents |
|---|---|
| `models/modelo_fisico_base.md` | Complete physical model specification |
| `models/derivacion_modelo_resolucion.md` | Step-by-step resolution derivation (88 equations) |
| `models/analisis_simulacion.md` | Simulation validation and measured results |
| `notes/modelo_conceptual.md` | Literature synthesis and comparison with prior work |
| `hypotheses/hipotesis_v5.md` | Latest validated hypothesis including off-axis targets |
