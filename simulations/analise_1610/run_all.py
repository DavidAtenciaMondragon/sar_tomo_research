"""Ejecuta todos los experimentos de analise_1610 en orden (A -> F y complementos)."""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
STEPS = [
    ("exp_A_referencias.py", {}),
    ("exp_B_barrido_dphi.py", {}),
    ("exp_B_fig_geometria.py", {}),
    ("t_resolucion_analitica.py", {}),
    ("tabla_resolucion_latex.py", {}),
    ("exp_C_brewster.py", {}),
    ("exp_C_figuras_enlace.py", {}),
    ("exp_C_fig_geometria.py", {}),
    ("exp_C_fig_trayectorias3d.py", {}),
    ("tabla_C_latex.py", {}),
    ("exp_D_biblioteca.py", {}),
    ("exp_D2_rejilla_densa.py", {}),
    ("exp_D3_rx_alineado.py", {}),
    ("exp_E_montecarlo.py", {"BISAR_EXP": "expE"}),
    ("exp_E_montecarlo.py", {"BISAR_EXP": "expE2"}),
    ("exp_F_trayectoria.py", {}),
    ("exp_F_fig_robustez.py", {}),
    ("exp_F2_calibracion_fase.py", {}),
]
for s, env in STEPS:
    t0 = time.time()
    print(f"=== {s} {env}", flush=True)
    subprocess.run([sys.executable, os.path.join(HERE, s)], check=True, cwd=HERE, env={**os.environ, **env})
    print(f"=== {s} OK ({time.time()-t0:.0f} s)", flush=True)
