# Escenario biestático con múltiples medios (multicapa de suelo)

**Fecha:** 2026-06-15
**Estado:** Implementado y validado numéricamente (pipeline `run_snell_multislab_pipeline.m`)
**Referencias:** `simulations/run_snell_pipeline.m` (caso base de 1 capa), `models/multiple_slab_tm.md`

---

## 1. Motivación

`run_snell_pipeline.m` modela la propagación a través de **un único medio de
suelo** (aire → suelo homogéneo, índice de refracción `n` constante). Este
documento describe la generalización a **N capas de suelo planas y
paralelas**, manteniendo el mismo criterio de arquitectura (generador de
señales GS + procesador PROC) y la misma coherencia física en el cálculo de
atraso, fase y energía transmitida.

El nuevo pipeline es `run_snell_multislab_pipeline.m`, que ejecuta:

- `simulations/gs/GS_snell_multislab_test_script.m` (generador de señales)
- `simulations/proc/PROC_snell_multislab_test_script.m` (procesador, back-projection)

con el set de parámetros `*_multislab.json` en `simulations/parametros/`.

---

## 2. Escenario de la primera simulación

```
                 Tx, Rx (trayectoria espiral)
                        \\        //
   ============================================  z = 0   (superficie)
        Medio 1 (suelo):  n1 = 2,  d1 = 1 m
   --------------------------------------------  z = -1 m
        Medio 2 (suelo):  n2 = 3,  (semi-infinito)
                         *  <- target, z0 = -1.5 m
   ============================================
```

- **Aire**: `n0 = 1` (semi-infinito).
- **Capa de suelo 1**: índice de refracción `n1 = 2`, espesor `d1 = 1 m`.
- **Capa de suelo 2**: índice de refracción `n2 = 3`, se extiende desde
  `z = -1 m` hacia abajo (semi-infinita a efectos del modelo).
- **Target**: enterrado a `1.5 m` de profundidad desde la superficie, es
  decir a `0.5 m` dentro de la capa 2 (`n = 3`).

Esto se configura en `simulations/parametros/system_multislab.json` mediante
el arreglo `CapasSuelo`:

```json
"CapasSuelo": [
    {"n": 2, "d": 1},
    {"n": 3, "d": 1}
]
```

y en `simulations/parametros/target_multislab.json`:

```json
"pos": [0, 0, -1.5]
```

El espesor declarado para la **capa que contiene al target** (`d` de la
última entrada de `CapasSuelo` necesaria) no se usa para la geometría —
la profundidad real del target (`z0`) determina cuánto de esa capa atraviesa
el rayo. Solo deben listarse las capas hasta llegar a la del target.

---

## 3. Geometría del rayo: ley de Snell generalizada (N capas)

Para interfaces planas, horizontales e infinitas, la componente horizontal
del vector de onda se conserva en cada interfaz (ley de Snell, ver
`models/multiple_slab_tm.md`, Sección 2):

$$n_0 \sin\theta_0 = n_1 \sin\theta_1 = \dots = n_K \sin\theta_K = p$$

donde `p` es el **parámetro de rayo**, constante a lo largo de toda la
trayectoria (aire + todas las capas atravesadas).

Si el rayo atraviesa segmentos de altura/espesor vertical `h_i` (el primero
en aire, los siguientes en cada capa de suelo, y el último es la fracción de
la capa del target hasta la profundidad real del target), la **proyección
horizontal total** debe igualar la distancia horizontal `u` entre el
radar y el target:

$$u = \sum_{i=0}^{K} h_i \tan\theta_i, \qquad \theta_i = \arcsin\!\left(\frac{p}{n_i}\right)$$

`calculateSlantRangeMultislab.m` (en `simulations/common/`) resuelve esta
ecuación para `p` mediante **bisección vectorizada** (la función es
monótonamente creciente en `p`, por lo que la raíz es única). Una vez
encontrado `p`:

- Ángulo de incidencia en aire: `theta_inc = asin(p / n0)`.
- Rango oblicuo de cada segmento: `R_i = h_i / cos(theta_i)`.

Esta es la generalización directa del método de punto fijo de
`calculateSlantRange.m` (caso de 1 capa), pero formulada como una búsqueda en
1 parámetro (`p`) en vez de iterar sobre posiciones — válido para cualquier
número de capas porque el invariante de Snell es global (Sección 2 de
`multiple_slab_tm.md`).

### Caso degenerado de reflexión total interna

Si en algún segmento `p > n_i`, no existe ángulo real → el rayo no puede
alcanzar esa capa con esa geometría. La función limita la búsqueda de `p` al
rango `[0, min(n_i))`, consistente con la condición de TIR descrita en
`multiple_slab_tm.md` Sección 2.

---

## 4. Atraso (delay) y fase

Para cada enlace (Tx→target y target→Rx) se obtiene un vector de rangos
oblicuos `R_i` y sus índices de refracción `n_i` (con `n_0 = 1` para el tramo
en aire). El **camino óptico total** (suma de `n_i · R_i`, coherente con
`v_i = c/n_i`) es:

$$\text{path} = \sum_{i=0}^{K} n_i R_i$$

**Atraso total ida+vuelta** (Tx → target → Rx):

$$t = \frac{1}{c}\left(\text{path}_{Tx} + \text{path}_{Rx}\right)$$

**Fase acumulada**:

$$\phi = -\frac{2\pi}{\lambda_0}\left(\text{path}_{Tx} + \text{path}_{Rx}\right)$$

con `λ0 = c / f_portadora` (longitud de onda en vacío). Estas expresiones son
la generalización directa de las usadas en `GS_snell_test_script.m`
(`t = (1/c)*(r1Tx + n*r2Tx + r1Rx + n*r2Rx)`), reemplazando la suma de dos
términos (`r1 + n·r2`) por la suma sobre todos los segmentos/capas
atravesados.

---

## 5. Energía de transmisión (coeficientes TM multicapa)

Para el aporte de amplitud del target se calcula la **transmitancia de
potencia** de la onda al penetrar desde el aire hasta la capa que contiene al
target, usando el formalismo ABCD de `models/multiple_slab_tm.md`.

A diferencia del modelo original (que asume entrada y salida en vacío,
`Z_in = Z_out`, Sección 6), aquí:

- **Medio de entrada**: aire, `n_in = n0 = 1`.
- **Slabs intermedios**: todas las capas de suelo *por encima* de la capa del
  target (en este escenario, solo la capa 1: `n = 2`, `d = 1 m`).
- **Medio de salida**: la capa que contiene al target (`n_out = 3`),
  tratada como semi-infinita.

`multislabTMcoef.m` (en `simulations/common/`) implementa la versión
generalizada de las Secciones 7–8 **sin** la simplificación `Z_out = Z_in`
(explícitamente prevista como extensión en la Sección 11 del documento de
referencia):

$$\Gamma_{TM} = \frac{A Z_{out} + B - Z_{in}(C Z_{out} + D)}{A Z_{out} + B + Z_{in}(C Z_{out} + D)}$$

$$\tau_{TM} = \frac{2 Z_{out}}{A Z_{in} + B + C Z_{in} Z_{out} + D Z_{out}}$$

$$T_{TM} = |\tau_{TM}|^2 \cdot \frac{\mathrm{Re}\{Z_{in}\}}{\mathrm{Re}\{Z_{out}\}}, \qquad R_{TM} = |\Gamma_{TM}|^2, \qquad R_{TM}+T_{TM}=1$$

El ángulo de incidencia `theta_in` usado en `multislabTMcoef` es el mismo
`theta_inc` (en aire) obtenido de `calculateSlantRangeMultislab` para cada
posición del radar — por lo tanto varía punto a punto a lo largo de la
trayectoria espiral.

### Aplicación al dato crudo

El coeficiente complejo `Tau_TM` (transmisión de campo) se aplica como factor
de amplitud/fase adicional en cada enlace:

```matlab
auxData(IND) = strTarget.rcs .* Tau_Tx .* Tau_Rx .* exp(1i*phi);
```

donde `Tau_Tx` y `Tau_Rx` son los coeficientes de transmisión evaluados en el
ángulo de incidencia de cada trayectoria (Tx y Rx, respectivamente). Esto
captura, de forma coherente con `t` y `phi`, la atenuación y el desfase que
introduce la interfaz aire-suelo (y cualquier capa intermedia) sobre la señal
que efectivamente alcanza al target y regresa.

**Simplificación adoptada (documentada explícitamente):** se usa el mismo
`Tau_TM` (calculado "de entrada", aire → capa del target) para ambos sentidos
de cada enlace (ida y vuelta dentro del enlace Tx↔target y target↔Rx). Una
generalización futura podría calcular por separado el coeficiente de
"salida" (capa del target → aire) resolviendo el sistema con `Z_in` y `Z_out`
intercambiados.

---

## 6. Verificación numérica

Pruebas de consistencia ejecutadas sobre `multislabTMcoef.m` y
`calculateSlantRangeMultislab.m` para el escenario `n=[2,3]`, `d=[1,1]`,
`z0=-1.5`:

- **Cierre geométrico**: `sum(h_i * tan(theta_i)) == u` (distancia horizontal)
  se verifica exactamente.
- **Invariante de Snell**: `R_i = h_i / cos(theta_i)` coincide con el valor
  devuelto por la función para cada segmento.
- **Conservación de energía**: `R_TM + T_TM = 1` (error de máquina, `~1e-16`).
- **Caso límite Fresnel**: con 0 slabs intermedios e incidencia normal,
  `Gamma_TM = (n_in - n_out)/(n_in + n_out)`, reproduciendo Fresnel clásico
  (`n_in=1, n_out=3` → `Gamma = -0.5`, ✓).

Ejecución completa del pipeline (`run_snell_multislab_pipeline.m`) con la
grilla `gridSnell` (`dxy=0.04, lxy=0.12, dz=0.08, lz=0.24`, centrada en el
target):

```
Target ubicado en la capa 2 (n = 3.000).
Pico detectado en: X = 0.0000 m, Y = 0.0400 m, Z = -1.5800 m
Resolucion X (-3 dB): 0.0311 m
Resolucion Y (-3 dB): 0.0247 m
Resolucion Z (-3 dB): 0.3152 m
```

El pico de back-projection queda correctamente localizado en torno a la
posición real del target `(0, 0, -1.5)`, con desviaciones del orden del
tamaño de celda de la grilla (`dxy=0.04`, `dz=0.08`).

---

## 7. Estructura de archivos nuevos

| Archivo | Rol |
|---|---|
| `simulations/run_snell_multislab_pipeline.m` | Orquestador GS→PROC |
| `simulations/gs/GS_snell_multislab_test_script.m` | Generador de señales (multicapa) |
| `simulations/proc/PROC_snell_multislab_test_script.m` | Procesador back-projection (multicapa) |
| `simulations/common/calculateSlantRangeMultislab.m` | Geometría de rayo N-capas (Snell + Fermat) |
| `simulations/common/multislabTMcoef.m` | Coeficientes TM (ABCD) con `Z_in != Z_out` |
| `simulations/parametros/system_multislab.json` | Define `CapasSuelo` (lista de `{n, d}`) y grilla |
| `simulations/parametros/target_multislab.json` | Posición y RCS del target |
| `simulations/parametros/radarTx_multislab.json`, `radarRx_multislab.json` | Trayectorias Tx/Rx (idénticas a `_espiral`) |

Salidas (figuras, raw data, datos procesados) se guardan en
`simulations/io/snell_multislab/`. En particular:

- `trayectorias_snell_multislab.png`: trayectorias Tx/Rx y, en la leyenda, las
  interfaces del medio multicapa (superficie, interfaz capa 1|2 y límite
  declarado de la última capa).
- `perfil_capas_multislab.png`: esquema 2D (perfil) del medio multicapa con
  cada capa coloreada, su `n` y espesor `d`, y la posición del target —
  pensado para visualizar de forma clara la estructura de capas, ya que en la
  vista 3D de la trayectoria el espesor de las capas (~1-2 m) es
  imperceptible frente a la escala de la espiral (~100-170 m).

---

## 8. Limitaciones y extensiones futuras

- Interfaces **planas, horizontales e infinitas** (sin DEM/topografía), igual
  que el caso de 1 capa.
- Medios **sin pérdidas** (`eps_r` real) y **no magnéticos**, igual que
  `multiple_slab_tm.md`.
- Polarización **TM** únicamente.
- El target debe estar contenido en alguna de las capas declaradas en
  `CapasSuelo`; si su profundidad excede el espesor total declarado, la
  función `calculateSlantRangeMultislab` lanza un error explícito.
- El factor de transmisión se aplica de forma simplificada (mismo `Tau_TM`
  para ida y vuelta de cada enlace, ver Sección 5). Una extensión natural es
  calcular también la transmitancia de "salida" (capa del target → aire) para
  cada enlace por separado.
- El pipeline soporta de forma genérica `M >= 1` capas (no solo 2): basta
  extender `CapasSuelo` en el JSON. La geometría (`calculateSlantRangeMultislab`)
  y los coeficientes TM (`multislabTMcoef`) ya están escritos para `N` capas
  arbitrarias.
