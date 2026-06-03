# Análisis de Representatividad de la Simulación MATLAB

**Fecha:** 2026-06-01 (actualizado tras correcciones del usuario)  
**Script analizado:** `simulations/run_snell_pipeline.m` → `GS_snell_test_script.m` + `PROC_snell_test_script.m`  
**Configuración activa:** `paramSuffix = '_espiral'`  
**Escenario objetivo:** SAR biestático · Vuelo helicoidal · Dos medios (aire/suelo) · Interfaz plana · Objeto enterrado

---

## Estado actual — Resumen rápido

| Elemento verificado | Estado | Referencia |
|---------------------|:------:|-----------|
| Estructura PROC sin `end` prematuro | ✅ | Línea 188 cierra correctamente |
| Biestático real (NorthOffset TX=0°, RX=180°) | ✅ | `funcao_espiral` aplica offset |
| Índice de refracción n=2 | ✅ | GS línea 82 |
| `report3dBResolution` llamada | ✅ | PROC línea 186, genera CSV |
| Target puntual único | ✅ | `target_espiral.json`: `[20,0,-5]` |
| Interfaz plana en z=0 | ✅ | Ambos solvers asumen z_interfaz=0 |
| Ley de Snell iterativa | ✅ | `calculateSlantRange` (GS) y `slantRange` (PROC) |
| Resultado de resolución disponible | ✅ | `bistatic_sar_snell_resolution_3dB_espiral.csv` |

---

## 1. Verificación de cada corrección

### 1.1 PROC_snell_test_script.m — sin bug estructural

El script cierra la función en la línea 188 (último `end`), con todo el código de procesamiento correctamente dentro:

```
línea   1: function PROC_snell_test_script(paramSuffix)
...
línea 186: report3dBResolution(X, Y, Z, outputData, resolutionCsvFile);
línea 187: 
línea 188: end   ← cierra la función correctamente
```

La función es invocada correctamente por `run_snell_pipeline.m`:
```matlab
PROC_snell_test_script(paramSuffix);
```

### 1.2 Biestático real — NorthOffset implementado

`funcao_espiral` ahora acepta el parámetro `NorthOffset` y lo aplica al ángulo inicial:

```matlab
theta = (Vt/Vr)*(log(R_top)-log(R_top-Vr*t)) + deg2rad(NorthOffset);
```

Con TX (`NorthOffset=0`) y RX (`NorthOffset=180`), las dos trayectorias quedan desplazadas 180° en azimut → en cada instante $t_n$, TX y RX están en **lados opuestos** de la hélice.

Los parámetros físicos son idénticos en TX y RX (mismo radio, altura, velocidad), por lo que la separación biestática es estrictamente azimutal: el ángulo biestático visto desde el objetivo varía con la posición del helicoide.

### 1.3 n=2 verificado

```matlab
% GS_snell_test_script.m, línea 82:
n = 2;
```

El PROC carga `n` del archivo `.mat` (guardado por el GS) y lo usa como `n2 = n`. Consistencia GS ↔ PROC confirmada.

### 1.4 Medición de resolución

`report3dBResolution` se llama con el grid de salida y el CSV de resultados:

```matlab
% PROC_snell_test_script.m, línea 185-186:
resolutionCsvFile = fullfile(processedDir, sprintf('bistatic_sar_snell_resolution_3dB%s.csv', paramSuffix));
report3dBResolution(X, Y, Z, outputData, resolutionCsvFile);
```

Resultado actual en `io/snell/processed/bistatic_sar_snell_resolution_3dB_espiral.csv`:

| Eje | Resolución −3dB | Límites [m] | Pico detectado |
|-----|----------------|------------|----------------|
| X | **0.138 m** | [19.931, 20.069] | X=20, Y=0, Z=−5 |
| Y | **0.140 m** | [−0.070, +0.070] | X=20, Y=0, Z=−5 |
| Z | **0.301 m** | [−5.168, −4.867] | X=20, Y=0, Z=−5 |

El pico se detecta exactamente en la posición del target `[20, 0, −5]`. ✅

### 1.5 Target puntual

```json
{"target": [{"pos": [20, 0, -5], "rcs": 1}]}
```

Un único dispersor puntual a 5 m de profundidad, desplazado 20 m en X respecto al eje de la hélice. El GS lo procesa como un solo punto (no hay `lineSegment`).

---

## 2. Consideraciones físicas restantes

### 2.1 ⚠️ Frecuencia 10 GHz no penetra suelos reales a 5 m

Los parámetros `_espiral` usan `FreqPortadora = 10e9` Hz (banda X, λ=3 cm).

**Profundidad de penetración en suelo típico (ε_r=4, σ≈0.01 S/m):**

$$\delta_{skin} = \frac{1}{\alpha} = \frac{c}{\omega\sqrt{\varepsilon_r/2\,(\sqrt{1+(\sigma/\omega\varepsilon)^2}-1)}}$$

Para 10 GHz y suelo húmedo: $\delta_{skin} \approx$ 2–5 cm. Para 425 MHz (P-band): $\delta_{skin} \approx$ 0.5–2 m.

| Frecuencia | Banda | λ | Penetración típica | Target a 5 m |
|------------|-------|---|--------------------|--------------|
| 10 GHz | X | 3 cm | ~2–5 cm | ❌ Físicamente irreal |
| 425 MHz | P | 70 cm | ~0.5–2 m | ✅ Adecuado |

**La simulación no modela atenuación dieléctrica** (el modelo asume solo cambio de velocidad, sin absorción), por lo que es matemáticamente válida pero **no representa un escenario físico realista** a 10 GHz para un objeto a 5 m de profundidad. Si el propósito es estudiar la geometría de la PSF y la resolución, esto puede aceptarse como simplificación. Si el propósito es validar la detectabilidad del objeto, la frecuencia debe cambiarse a P-band (configuración `_espiral_EMIRADOS`).

### 2.2 ⚠️ La geometría bistática NorthOffset=180° tiene un k-espacio diferente al monoestático

Con TX y RX siempre en lados opuestos del helicoide, los vectores $\hat{\mathbf{e}}_2^{TX}$ y $\hat{\mathbf{e}}_2^{RX}$ tienen componentes horizontales que se **cancelan parcialmente**, mientras que las componentes verticales se **refuerzan**. Esto produce:

- **Cobertura horizontal** ($k_x, k_y$): reducida respecto a monoestático → peor resolución horizontal
- **Cobertura vertical** ($k_z$): amplificada (~doble que monoestático) → mejor resolución vertical

**Comparación con predicciones monoestáticas de `resolucion_v2.md`:**

Para target en el eje (x_P=0), la cancelación sería total → sin resolución horizontal. Para target off-axis (x_P=20 m), la cancelación es parcial y la resolución medida es:

$$\delta_{xy}^{med} \approx 0.14\text{ m} \gg \delta_{xy}^{mono} \approx 0.003\text{ m} \quad\text{(fórmula monoestática, 10 GHz)}$$

Este comportamiento es **físicamente coherente** con la geometría bistática NorthOffset=180°: las fórmulas analíticas de `resolucion_v2.md` aplican al caso monoestático o biestático genérico, pero no predicen directamente este caso de TX/RX diametralmente opuestos.

**Para el caso biestático con NorthOffset=180°**, el vector de número de onda instantáneo es:

$$\mathbf{k}(f,t) = \frac{2\pi f\,n_2}{c}\!\left(\hat{\mathbf{e}}_2^{TX}(\alpha) + \hat{\mathbf{e}}_2^{RX}(\alpha+\pi)\right)$$

La componente horizontal de $\hat{\mathbf{e}}_2^{TX}$ apunta desde $\mathbf{Q}_{TX}$ hacia P, y la de $\hat{\mathbf{e}}_2^{RX}$ apunta en sentido casi opuesto → cancelación parcial de $k_x$ y $k_y$.

### 2.3 ✅ Coherencia de la resolución vertical medida

La resolución vertical $\delta_z = 0.301$ m proviene de dos contribuciones:

**Contribución de ancho de banda** (B=50 MHz, n₂=2, ψ₀≈57°, cos θ_t,0≈0.91):

$$\Delta k_z^{(B)} \approx \frac{4\pi n_2 B}{c}\cos\theta_{t,0} = \frac{4\pi \times 2 \times 50\times10^6}{3\times10^8}\times 0.91 \approx 3.8\text{ m}^{-1}$$

$$\delta_z^{(B)} = \frac{2\pi}{3.8} \approx 1.65\text{ m} \quad\text{(si solo hubiera BW)}$$

**Contribución geométrica del helicoide cónico** (variación de look angle de ψ=49.7° a ψ=63.7°):

$$\Delta\psi \approx 14°, \quad \Delta\cos\theta_t \approx 0.030$$

$$\Delta k_z^{(geom)} \approx \frac{4\pi f_0 n_2}{c}\cdot\Delta\cos\theta_t = \frac{4\pi \times 10^{10}}{3\times10^8}\times 2 \times 0.030 \approx 25\text{ m}^{-1}$$

$$\delta_z^{(geom)} = \frac{2\pi}{25} \approx 0.25\text{ m}$$

La **contribución geométrica domina** sobre la de BW ($\Delta k_z^{(geom)} \gg \Delta k_z^{(B)}$). La suma da:

$$\Delta k_z \approx 3.8 + 25 = 28.8\text{ m}^{-1} \Rightarrow \delta_z \approx 0.22\text{ m}$$

El valor medido $\delta_z = 0.301$ m es del mismo orden (diferencia por el carácter off-axis del target, la geometría bistática y el criterio de resolución). ✅ **Coherente físicamente.**

### 2.4 ℹ️ Target off-axis — las fórmulas de `resolucion_v2.md` no aplican directamente

El target está en `[20, 0, -5]` (x_P=20 m ≠ 0). Las expresiones analíticas de `resolucion_v2.md` fueron derivadas bajo la hipótesis **"target en el eje de la hélice"** (x_P=y_P=0).

Para target off-axis, el look angle efectivo es diferente para cada azimut del sensor, la PSF se inclina ligeramente (efecto near-range demostrado por Góes 2022), y la simetría azimutal se rompe. Los resultados numéricos de la simulación son válidos; lo que no aplica directamente son las predicciones analíticas cerradas.

### 2.5 ℹ️ `ratio=1` hardcodeado, no usa `strSystem.ratioUp=2`

```matlab
% PROC_snell_test_script.m, línea 66:
ratio = 1;   % ← sistema tiene ratioUp=2 definido pero no se usa
```

Sin upsampling, la resolución se mide directamente sobre la grilla procesada. Con `ratio=2`, la grilla se interpola en rango antes del backprojection, mejorando la precisión del muestreo de la PSF. No es un error crítico pero podría mejorar ligeramente la precisión de la medición de δz.

---

## 3. Parámetros de la Configuración `_espiral`

| Parámetro | TX | RX | Unidad |
|-----------|----|----|--------|
| Frecuencia portadora | 10,000 | 10,000 | MHz |
| Ancho de banda | 50 | — | MHz |
| PRF | 400 | 400 | Hz |
| Vueltas | 2 | 2 | — |
| Radio menor (cima) | 147.5 | 147.5 | m |
| Radio mayor (base) | 172.5 | 172.5 | m |
| Altura mayor | 120 | 120 | m |
| Altura menor | 80 | 80 | m |
| NorthOffset | **0** | **180** | ° |
| n (índice de refracción) | 2 | 2 | — |

**Target:** `[20, 0, -5]` m · RCS = 1 · Profundidad real en suelo: 5 m

**Grid de procesado (`gridSnell`):**
- dxy = 0.04 m, span ±0.12 m → 7 puntos en X e Y (centrado en x=20, y=0)
- dz = 0.08 m, span ±0.24 m → 7 puntos en Z (centrado en z=−5)

---

## 4. Conclusión de Representatividad

| Elemento del escenario | Representado | Observación |
|------------------------|:---:|------------|
| Trayectoria helicoidal cónica | ✅ | TX y RX con parámetros idénticos, offset azimutal 180° |
| Configuración biestática | ✅ | NorthOffset diferente → trayectorias distintas |
| Dos medios homogéneos (n₂=2) | ✅ | n=2 en GS; cargado en PROC |
| Interfaz plana en z=0 | ✅ | Ambos solvers iterativos asumen z_interfaz=0 |
| Objeto enterrado (z<0) | ✅ | Target en z=−5 m |
| Refracción Snell en interfaz | ✅ | Solver iterativo convergente en GS y PROC |
| Formación de imagen 3D (BP) | ✅ | Loop sobre grid 7×7×7 en PROC |
| Corrección de fase biestática | ✅ | `R1t + n₂·R2t + n₂·R2r + R1r` en fase GS y PROC |
| Medición de resolución −3dB | ✅ | CSV generado: δx=14cm, δy=14cm, δz=30cm |
| Frecuencia adecuada para 5m profundidad | ⚠️ | 10 GHz (X-band): penetración ~cm; P-band necesaria para escenario realista |
| Atenuación dieléctrica en suelo | ❌ | No modelada — la simulación ignora la absorción del suelo |
| Patrón de antena | ❌ | No aplicado — amplitud uniforme |
| Pérdida de propagación (1/R²) | ❌ | No aplicada — RCS independiente de distancia |

**Veredicto:** La simulación es **representativa de la geometría** del escenario (helicoide cónico biestático, dos medios, refracción, objeto enterrado) y produce resultados de resolución coherentes con el modelo físico. La principal limitación es la **frecuencia de 10 GHz**, que no es realista para detección a 5 m de profundidad en suelo. Para validar el escenario completo con penetración real, usar la configuración `_espiral_EMIRADOS` (425 MHz, P-band).
