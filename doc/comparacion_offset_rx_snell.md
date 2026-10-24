# Interpretación de Resultados: Offset Azimutal del Rx en `run_snell_pipeline.m`

**Propósito**: Interpretar y explicar los resultados de la comparación entre dos configuraciones de la hélice del receptor (`NorthOffset = 180°`, geometría "especular", vs. `NorthOffset = 0°`) en el pipeline `run_snell_pipeline.m`, evaluando tanto la resolución de imagen (back-projection) como la potencia biestática recibida (ecuación de radar completa).

---

## 1. Configuración evaluada

### 1.1 Geometría del Tx (fija en ambos escenarios)

La hélice del Tx se derivó del óptimo `psi0* = 53.39°` hallado por `run_otimiza_tx.m` (escenario `R0fixed`, γ=0.30, corrida `otimiza_tx_20260912_092942`), aplicando la condición de hélice óptima β=ψ₀ (`otimiza_tx.tex`, §2.3) sobre la geometría nominal (R₀=188.68 m, B_helix=47.17 m) de `radarTx_espiral.json`:

| Parámetro | Valor |
|---|---|
| RaioMenorEsp / RaioMaiorEsp | 137.40 m / 165.53 m |
| AltMenorEsp / AltMaiorEsp | 93.58 m / 131.45 m |
| Vt | 113.50 m/s (ajustado para igualar la duración T del pase Rx, ver §4.4) |
| NorthOffset | 0° |

### 1.2 Geometría del Rx (dos escenarios)

Radio/altura sin modificar (172.5/147.5 m, 120/80 m), solo cambia el ángulo de inicio azimutal:

| Escenario | NorthOffset Rx | Suffix de parámetros |
|---|---|---|
| Especular | 180° | `_espiral_offset180` |
| Alineado | 0° | `_espiral_offset0` |

### 1.3 Target

Punto único en `[0, 0, −5]` m, RCS σ=1 m² — **nótese que el target está exactamente sobre el eje vertical de rotación de ambas hélices** (x=y=0). Esta propiedad es la clave de la interpretación en §3.

---

## 2. Resultados numéricos

### 2.1 Resolución (−3 dB, back-projection numérico)

| Eje | offset=180° | offset=0° | Razón |
|---|---|---|---|
| X | 0.1912 m | 0.01518 m | 12.6× |
| Y | 0.1948 m | 0.01519 m | 12.8× |
| Z | 0.2056 m | 0.2056 m | 1.0× |

### 2.2 Potencia biestática recibida (ecuación de radar completa)

| Escenario | P_r máx | P_r media | P_r mín |
|---|---|---|---|
| offset=180° | −74.99 dBm | −75.09 dBm | −75.36 dBm |
| offset=0° | −74.99 dBm | −75.09 dBm | −75.36 dBm |

Diferencia: **0.00 dB** en las tres estadísticas — las series temporales P_r(t) son bit-a-bit idénticas entre ambos escenarios.

---

## 3. Por qué la resolución cambia radicalmente y la energía no cambia nada

Este es el resultado central y, a primera vista, contraintuitivo: **mover el Rx 180° en azimut no cambia la potencia recibida en absoluto, pero cambia la resolución horizontal en un orden de magnitud**. Ambos efectos tienen la misma causa raíz — la posición del target sobre el eje de rotación — pero actúan sobre cantidades físicas distintas.

### 3.1 Por qué la energía es invariante: simetría axial del presupuesto de potencia

La ecuación de radar biestática implementada (idéntica a `flightpath/objective_function.m`, ver `doc/verifica_plan_de_vuelo.md` §2.4) es:

$$P_{r} = \frac{P_t \cdot G_t \cdot G_r \cdot \sigma \cdot T_1 \cdot T_2 \cdot \lambda^2}{(4\pi)^3 \cdot R_T^2 \cdot R_R^2}$$

Cada factor de esta ecuación depende **únicamente de la distancia radial ρ(t) y la altura z(t)** de la plataforma respecto al eje vertical que pasa por el target — nunca del ángulo azimutal absoluto θ(t):

- **R_T, R_R** (rango total refractado): dependen de la distancia horizontal `D = ‖posición_xy − target_xy‖`. Como el target está en (0,0), `D = ρ(t) = R(t)` (el radio instantáneo de la espiral), sin importar en qué azimut se encuentre la plataforma.
- **T₁, T₂** (transmitancia de Fresnel/TM): dependen del ángulo de incidencia en el punto de refracción de Fermat, que a su vez depende solo de (D, z) — la misma invariancia.
- **G_t** (ganancia Tx): el boresight apunta siempre al origen (`[0,0,0]`, prácticamente coincidente con el target a 5 m de profundidad), así que el ángulo fuera de eje depende de la geometría (Tx, target) — invariante ante rotación del Rx, y de hecho invariante entre ambos escenarios porque la trayectoria del Tx es **la misma** en los dos casos.
- **G_r** (ganancia Rx): en `objective_function.m` se asume que la antena del Rx apunta siempre exactamente al target (antena idealmente orientada), por lo que G_r es constante en el pico del patrón independientemente de la posición del Rx.

`funcao_espiral.m` construye `R(t)` y `Pz(t)` en función solo del tiempo (vía `Vr`, `Vd`); el `NorthOffset` solo desplaza `θ(t) = ... + NorthOffset`, es decir, **rota la trayectoria alrededor del eje z sin alterar ρ(t) ni z(t)**. Como el target está exactamente sobre ese eje, dos hélices Rx idénticas salvo el `NorthOffset` recorren *el mismo conjunto de pares (ρ, z) en el mismo instante t* — sencillamente en un plano azimutal distinto. El resultado es que **P_r(t) es matemáticamente idéntica en ambos escenarios**, no una coincidencia numérica.

**Consecuencia práctica**: en esta geometría (target sobre el eje), el offset azimutal del Rx es un grado de libertad "gratis" desde el punto de vista energético — no hay ningún costo de potencia asociado a elegirlo.

### 3.2 Por qué la resolución sí cambia: diversidad angular del par Tx–Rx

La resolución horizontal (derivada en `doc/explicacion_resolucion.tex`, §7.2, y usada en `flightpath/calculateBistaticResolution.m`) es:

$$\delta_{xy}(\Delta\phi) = \frac{0.60\,\lambda}{\pi\,\sin\psi_0\cdot|\cos(\Delta\phi/2)|}$$

donde **Δφ es la separación azimutal entre Tx y Rx vista desde el target** — justamente la cantidad que el `NorthOffset` relativo controla directamente. A diferencia de la potencia, la resolución de un sistema SAR/tomográfico depende de la **diversidad angular de observación** (cuánto "rodean" las plataformas al blanco), no de la distancia radial pura. Es el equivalente angular de la apertura sintética: dos aperturas casi colineales (Δφ→180°, Tx y Rx en lados opuestos vistos desde el target) degeneran, porque `cos(Δφ/2)=cos(90°)→0` hace `δ_xy→∞`.

Con el Tx en `NorthOffset=0°`:
- **offset=180°** → Δφ≈180° → término `|cos(90°)|≈0` → resolución horizontal fuertemente degradada (tabla de `doc/verifica_plan_de_vuelo.md` §5.4 predice exactamente este comportamiento para Δφ→180°).
- **offset=0°** → Δφ≈0° → término `|cos(0°)|=1` → resolución horizontal óptima (mínima δ_xy alcanzable para ese ψ₀).

Esto coincide con el resultado empírico del back-projection (0.015 m vs 0.19 m, factor ~12.7×), que es consistente con la razón teórica `|cos(0°)|/|cos(90°+ε)|` para el Δφ efectivo entre ambas trayectorias.

La resolución vertical δ_z, en cambio, **no depende de Δφ** (solo de ψ₀ del Rx y B_helix, ver `doc/verifica_plan_de_vuelo.md` §2.6.1), lo que explica por qué es idéntica (0.2056 m) en ambos escenarios.

### 3.3 Síntesis: dos cantidades físicas, dos simetrías distintas

| Cantidad | Depende de | Invariante ante NorthOffset porque... |
|---|---|---|
| P_r (energía) | ρ(t), z(t) (posición **radial**) | target sobre el eje de rotación ⇒ ρ(t) no depende del azimuto |
| δ_xy (resolución horiz.) | Δφ (separación **angular** Tx–Rx) | NorthOffset controla Δφ directamente |
| δ_z (resolución vert.) | ψ₀, B_helix | no depende de azimuto en absoluto |

El offset azimutal del Rx es, en esta geometría, un parámetro que **solo afecta resolución, nunca energía**. Esto convierte la elección en una decisión sin trade-off: no hay ninguna razón energética para preferir `offset=180°` sobre `offset=0°`, y sí una razón de resolución (factor ~13×) para preferir `offset=0°`.

---

## 4. Hallazgos metodológicos adicionales

### 4.1 `run_snell_pipeline.m` no es sensible a la energía real

`GS_snell_test_script.m` inyecta amplitud constante (`strTarget.rcs`) al generar los datos crudos, sin modelar `T₁`, `T₂`, ganancia de antena ni *spreading loss* (`common/...`, línea `auxData(IND) = strTarget.rcs.*exp(1i*phi)`). Por diseño, el pico de la imagen back-projected (`report3dBResolution.m`) es **insensible a la geometría Tx/Rx** — confirmado empíricamente: el pico de `|outputData|` fue exactamente 411562 en ambos escenarios (0.000 dB de diferencia), a pesar de la enorme diferencia de resolución. Este pipeline sirve para validar foco/geometría de imagen, no presupuesto de enlace. La comparación de energía físicamente objetiva requirió replicar la ecuación de radar de `flightpath/objective_function.m` en un script separado (`verifica_energia_snell_offsets.m`).

### 4.2 Bug corregido en `flightpath/getAntennaGain.m`

Durante el cálculo de energía apareció un outlier: `Gr ≈ 1.7×10⁻⁵⁰` en la muestra `t=0` del escenario `offset=180°`, mientras que `T₁`, `T₂`, `Gt`, `R_T`, `R_R` eran idénticos a los de `offset=0°` en ese mismo instante. Causa raíz: el *wrap* de ángulo azimutal en `getAntennaGain.m` solo corregía el caso `Δaz > π` (restando 2π) pero no el caso simétrico `Δaz < −π` (sumando 2π):

```matlab
% Antes (incompleto):
az_delta_pattern(az_delta_pattern > pi) = az_delta_pattern(az_delta_pattern > pi) - 2*pi;

% Después (wrap simétrico correcto):
az_delta_pattern(az_delta_pattern > pi)  = az_delta_pattern(az_delta_pattern > pi)  - 2*pi;
az_delta_pattern(az_delta_pattern < -pi) = az_delta_pattern(az_delta_pattern < -pi) + 2*pi;
```

En el caso de borde (`az_Tx_ref≈0`, `az_Tx_tg≈2π−ε`), `Δaz≈−2π+ε` no se envolvía correctamente, cayendo en el extremo del patrón de radiación gaussiano (ganancia casi nula). Esta función es compartida por `objective_function.m`, `run_plano_de_voo.m` y `run_otimiza_tx.m`; cualquier resultado previo con una caída de potencia puntual e inexplicable en una geometría de borde específica podría estar afectado por este mismo artefacto.

### 4.3 Sincronización de longitudes Tx/Rx en `funcao_espiral.m`

Al aplicar la geometría óptima del Tx (radios más pequeños que el Rx original), la duración del pase `T = (R_top−R_bottom)·N·2π/(Vt·ln(R_top/R_bottom))` dejó de coincidir entre Tx y Rx, causando un error de tamaños incompatibles en `GS_snell_test_script.m` (línea `t = (1/c)*(r1Tx + n*r2Tx + r1Rx + n*r2Rx)`, que asume igual número de muestras). Se corrigió ajustando `Vt` del Tx (113.50 m/s) para igualar `T` entre ambas hélices, preservando `NumVoltasEsp=2` en ambas.

### 4.4 Bug de nomenclatura de figuras

`GS_snell_test_script.m` y `PROC_snell_test_script.m` guardaban las figuras (`trayectorias_snell.png`, `slice_abs_output_snell.png`, `isosurface_output_snell.png`) con nombre fijo, sin sufijo de escenario — a diferencia de los `.mat`/`.csv` de salida, que sí incluyen `paramSuffix`. Corridas sucesivas con distinto `paramSuffix` sobrescribían las figuras de la corrida anterior. Se corrigió incluyendo `paramSuffix` en los cuatro nombres de archivo de figura.

---

## 5. Conclusión y recomendación

Para la geometría evaluada (target sobre el eje de rotación de las hélices Tx/Rx), **`NorthOffset=0°` domina estrictamente a `NorthOffset=180°`**: misma potencia recibida (0.00 dB de diferencia, por simetría axial del presupuesto de enlace) y resolución horizontal ~13× mejor (por diversidad angular Δφ≈0° vs Δφ≈180°). No existe trade-off físico que favorezca la configuración especular en este caso — su motivación original (maximizar Δφ pensando en diversidad de aspecto) no se traduce en ninguna ventaja de potencia porque el target está centrado en el eje, y sí introduce una penalización severa de resolución vía el término `1/|cos(Δφ/2)|`.

**Generalización**: esta conclusión es específica de targets sobre el eje de simetría. Para un target descentrado (x,y ≠ 0), ρ(t) — y por tanto P_r(t) — sí dependería del azimuto relativo entre Tx y Rx, reintroduciendo un trade-off potencia-resolución genuino, análogo al descrito en `doc/verifica_plan_de_vuelo.md` §2.6 para el problema de optimización continua del Rx.

---

## 6. Trazabilidad

| Artefacto | Ruta |
|---|---|
| Parámetros Tx/Rx/target/sistema por escenario | `simulations/parametros/{radarTx,radarRx,target,system}_espiral_offset{180,0}.json` |
| Script de energía física | `simulations/verifica_energia_snell_offsets.m` |
| Series temporales P_r, T1, T2, Gt, Gr, R_T, R_R | `simulations/io/snell/processed/energia_biestatica_espiral_offset{180,0}.csv` |
| Resumen energía | `simulations/io/snell/processed/comparacion_energia_offsets_resumen.csv` |
| Resolución −3dB | `simulations/io/snell/processed/bistatic_sar_snell_resolution_3dB_espiral_offset{180,0}.csv` |
| Figuras (trayectoria, slice, isosuperficie, comparación Pr(t)) | `simulations/io/snell/figure/*_espiral_offset{180,0}.png`, `comparacion_energia_offsets.png` |
| Fix aplicado | `simulations/flightpath/getAntennaGain.m` (wrap azimutal simétrico) |
