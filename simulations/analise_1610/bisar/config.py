"""Carga de parámetros: cada experimento tiene un JSON que hereda de common.json."""
import copy
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                       # simulations/analise_1610
PARAMS = os.path.join(ROOT, "params")
SIM = os.path.dirname(ROOT)                        # simulations/
IO = os.path.join(SIM, "io", "analise_1610")
FIGDOC = os.path.join(os.path.dirname(SIM), "doc", "figures", "relatorio_1610")


def _merge(a, b):
    out = copy.deepcopy(a)
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load(name):
    """Carga params/<name>.json fusionado sobre params/common.json."""
    with open(os.path.join(PARAMS, "common.json"), encoding="utf-8") as f:
        base = json.load(f)
    if name == "common":
        return base
    with open(os.path.join(PARAMS, f"{name}.json"), encoding="utf-8") as f:
        exp = json.load(f)
    return _merge(base, exp)


def outdir(name):
    d = os.path.join(IO, name)
    os.makedirs(d, exist_ok=True)
    os.makedirs(FIGDOC, exist_ok=True)
    return d


def save_resolved(p, name):
    """Guarda el conjunto de parámetros efectivamente usado (reproducibilidad)."""
    d = outdir(name)
    with open(os.path.join(d, f"params_resueltos_{name}.json"), "w", encoding="utf-8") as f:
        json.dump(p, f, indent=2, ensure_ascii=False)
