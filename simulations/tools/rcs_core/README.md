# estudio_rcs — RCS monoestática de objetos PEC simples

Instalación (Windows, Python 3.13):

```
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
```

## Script principal: `simular.py`

```
.venv\Scripts\python simular.py cubo --dim 1.0 --freq 400e6
```

Calcula la RCS TM-TM en todos los ángulos (az 0–360°, el −90…90°) y abre en el navegador
un HTML interactivo con el escenario 3D (radar, k, campo E), la malla, el diagrama 3D de RCS,
el mapa az×el y los cortes. Opciones: `--pol TE`, `--metodo po|mom|auto`, `--paso 2`,
`--vista AZ EL`, `--dim` (cilindro/cono: `radio largo`), `--archivo x.stl`.
Ayuda completa: `.venv\Scripts\python simular.py -h`.

En MATLAB: `r = load('resultados/<nombre>.mat')`. Las matrices `r.rcs_dbsm` y `r.rcs_m2` son
Nel × Naz (filas = `r.el_deg`, columnas = `r.az_deg`), así que `imagesc(r.az_deg, r.el_deg, r.rcs_dbsm)`
funciona directamente. Ver `leer_rcs.m`.

| Archivo | Qué hace |
|---|---|
| `simular.py` | **Main**: simulación y visualización 3D interactiva (HTML); guarda también `.mat` y `.npz`. |
| `leer_rcs.m` | Ejemplo de MATLAB: carga el `.mat`, muestra un resumen, interpola en un ángulo y dibuja el mapa, los cortes, la malla y el diagrama 3D. |
| `rcs/po.py` | Solver de Óptica Física (PO): integral exacta por triángulo y núcleo numba. Calcula un mapa completo de az×el en menos de 1 s. |
| `rcs/mom.py` | Solver de onda completa (MoM, EFIE) basado en **bempp-cl**: una factorización LU y después cualquier número de ángulos. Da VV, HH, VH y HV. |
| `rcs/geometry.py` | Mallas de cubo, caja, cilindro, esfera, cono y placa con gmsh, y carga de STL/OBJ con trimesh. |
| `rcs/analytic.py` | Serie de Mie de la esfera PEC (referencia exacta). |
| `comparar_po_mom.py` | Cortes de PO frente a MoM para cubo y cilindro → `resultados/po_vs_mom.png`. |
| `barrido_frecuencia.py` | RCS frente a frecuencia con PO, MoM y Mie → `resultados/barrido_frecuencia.png`. |


Convención: el azimut se mide en el plano xy desde +x y la elevación desde el plano xy.
V = TM (E en el plano que contiene z y la dirección de propagación) y H = TE (E horizontal); TM-TM = VV. Las unidades son SI; σ en m² (dBsm con `dbsm()`).

## Cuándo usar cada método
- **PO**: el objeto debe medir ≳ 2–3 λ. Da < 0.5 dB en los lóbulos especulares (caras
  normales, generatriz del cilindro). No modela bordes, ondas de superficie ni rebotes
  múltiples, y en monoestática HH = VV. Fuera de especular puede equivocarse en 5–20 dB.
- **MoM**: es exacto (validado frente a Mie con un error de 0.1–0.25 dB con malla λ/10). La memoria crece como 16·N²
  bytes (N ≈ 1.5 × triángulos), así que en un PC sirve para objetos de hasta ~2–3 λ.
