# Verificación y Crítica del Script `run_plano_de_voo.m`

**Propósito del script**: Optimizar la trayectoria del receptor en un sistema SAR bistático, buscando la posición de Rx que maximiza la potencia recibida proveniente de targets subterráneos, explotando el ángulo de Brewster para maximizar la transmitancia en la interfaz aire–suelo (polarización TM).

---

## 1. Descripción general del sistema

### 1.1 Parámetros clave

| Parámetro | Valor | Notas |
|---|---|---|
| Frecuencia portadora | 400 MHz (banda P) | Corregido de 10 GHz; ver §4.2 |
| Longitud de onda | λ = 0.75 m | c/400 MHz |
| Potencia Tx | 100 kW | |
| Velocidad del vehículo | 120 m/s | |
| PRF | 400 Hz | |
| Altitud de la espiral | 80–120 m | Radio 147.5–172.5 m |
| Índice de refracción suelo | n₂ = 2 | n₁ = 1 (aire); εᵣ = 4, lossless |
| Ángulo de Brewster teórico | θ_B ≈ 63.43° | arctan(n₂/n₁) = arctan(2) |
| Profundidad de los targets | 2–10 m | 4×4×4 = 64 puntos |
| RCS | σ = 1 m² | Isótropo |
| Radio de búsqueda Rx | ±200 m | Espacio de búsqueda 2D |

### 1.2 Flujo del script

```
1. Cargar parámetros → trayectoria espiral Tx
2. Decimar posiciones Tx (factor 200) → ~50 posiciones de evaluación
3. Para cada posición Tx decimada:
   a. Calcular ganancias Tx hacia todos los targets
   b. Construir geometría de búsqueda Brewster
   c. Lanzar fmincon desde 5 puntos candidatos
   d. Seleccionar Rx_opt = argmax(Pr_total)
4. Interpolar (PCHIP) trayectoria completa del Rx
5. Guardar resultados y visualizar
```

---

## 2. Ecuacionamiento del sistema

### 2.1 Ángulo de Brewster (polarización TM)

Para una interfaz plana entre dos medios dieléctricos lossless con índices n₁ y n₂, el ángulo de incidencia para el cual el coeficiente de reflexión TM (p-polarización) es nulo es:

$$\theta_B = \arctan\!\left(\frac{n_2}{n_1}\right)$$

Para n₁ = 1, n₂ = 2:

$$\theta_B = \arctan(2) \approx 63.43°$$

A este ángulo, la onda incidente TM se transmite **íntegramente** al segundo medio (T = 1, R = 0), lo que constituye la condición ideal para penetración subsuperficial. Nótese que θ_B es independiente de la frecuencia para medios lossless.

### 2.2 Coeficientes de Fresnel (TM)

La ley de Snell rige el ángulo de transmisión:

$$n_1 \sin\theta_1 = n_2 \sin\theta_2 \quad \Rightarrow \quad \theta_2 = \arcsin\!\left(\frac{n_1}{n_2}\sin\theta_1\right)$$

Los coeficientes de campo eléctrico para polarización TM son:

$$r_{TM} = \frac{n_2\cos\theta_1 - n_1\cos\theta_2}{n_2\cos\theta_1 + n_1\cos\theta_2}, \qquad t_{TM} = \frac{2n_1\cos\theta_1}{n_2\cos\theta_1 + n_1\cos\theta_2}$$

Las transmitancia y reflectancia de **potencia**, que conservan energía (T + R = 1), son:

$$T = \frac{n_2\cos\theta_2}{n_1\cos\theta_1}|t_{TM}|^2, \qquad R = |r_{TM}|^2$$

Verificación numérica en θ₁ = θ_B = arctan(2):

- θ₂ = arcsin(sin 63.43°/2) ≈ arcsin(0.447) ≈ 26.57°
- r_TM = (2·cos 63.43° − 1·cos 26.57°) / (2·cos 63.43° + 1·cos 26.57°) = (0.894 − 0.894)/1.788 = **0** ✓
- t_TM = 2·1·cos(63.43°) / 1.788 = 0.894/1.788 ≈ **0.5**
- T = (2·cos 26.57°)/(1·cos 63.43°)·(0.5)² = (1.788/0.447)·0.25 = 4·0.25 = **1** ✓

La implementación en `calculateTMcoef.m` es **correcta**.

### 2.3 Principio de Fermat y punto de refracción

El punto de cruce en la interfaz z = 0 se determina minimizando el tiempo óptico total:

$$\mathcal{T}(x, y) = n_1\,\|\mathbf{P}_1 - \mathbf{r}\| + n_2\,\|\mathbf{r} - \mathbf{P}_2\|, \quad \mathbf{r} = (x, y, 0)$$

donde P₁ y P₂ son las posiciones del emisor y receptor respectivamente. La condición de estacionariedad ∂T/∂x = ∂T/∂y = 0 es equivalente a la **ley de Snell** para una superficie plana. La implementación en `calculateRefractionPointFermat.m` mediante `fminunc` con quasi-Newton es **correcta**.

### 2.4 Ecuación de radar bistática (forma corregida)

La función objetivo `objective_function.m` evalúa para cada target i la potencia bistática con propagación refractada en dos tramos:

$$\boxed{P_{r,i} = \frac{P_t \cdot G_t^{(i)} \cdot G_r^{(i)} \cdot \sigma \cdot T_1^{(i)} \cdot T_2^{(i)} \cdot \lambda^2}{(4\pi)^3 \cdot R_T^2 \cdot R_R^2}}$$

donde las distancias efectivas de spreading integran la longitud total de cada camino refractado:

$$R_T = R_1 + R_2 \quad \text{(Tx → interfaz → target)}, \qquad R_R = R_3 + R_4 \quad \text{(target → interfaz → Rx)}$$

con:
- R₁ = |Tx − P_int1|: tramo aéreo de ida
- R₂ = |P_int1 − target|: tramo subterráneo de ida
- R₃ = |target − P_int2|: tramo subterráneo de vuelta
- R₄ = |P_int2 − Rx|: tramo aéreo de vuelta
- T₁, T₂: transmitancias de potencia TM en cada cruce de interfaz

La potencia total acumulada sobre el volumen de targets es:

$$P_r^{total} = \sum_{i=1}^{N_{tg}} P_{r,i}$$

Esta formulación asegura que el factor de pérdida de espacio libre $1/(4\pi R)^2$ se aplica a la distancia total recorrida por el rayo en cada tramo, coherente con la física de propagación esférica en medios en los que la interfaz es un refractor —no una fuente secundaria.

### 2.5 Justificación del factor (4π)³

La ecuación de radar bistática estándar (derivada de la ecuación de Friis) descompone la potencia recibida como:

$$P_r = \underbrace{\frac{P_t G_t}{4\pi R_T^2}}_{\text{densidad de potencia en el target}} \cdot \sigma \cdot \underbrace{\frac{G_r \lambda^2}{(4\pi)^2 R_R^2}}_{\text{apertura efectiva del Rx}}$$

El producto de los tres factores geométricos $4\pi \cdot (4\pi)^2 = (4\pi)^3$ surge naturalmente. Las transmitancias T₁ y T₂ se agregan como modificadores de la densidad de potencia incidente en la interfaz, y son dimensionalmente consistentes dado que son coeficientes de potencia (adimensionales, ∈ [0,1]).

---

## 3. Verificación de la coherencia del modelo

### 3.1 Coeficientes de Fresnel — CORRECTO ✓

Las ecuaciones en `calculateTMcoef.m` corresponden exactamente a la forma estándar de las ecuaciones de Fresnel para TM. La transmitancia de potencia conserva energía (T + R = 1) y el cálculo del ángulo de Brewster es correcto.

### 3.2 Fermat y ángulo de incidencia — CORRECTO ✓

El punto de refracción se determina correctamente minimizando el tiempo óptico. El ángulo de incidencia se calcula respecto a la normal superficial (ẑ = [0,0,1]) usando el valor absoluto del producto escalar, lo que es correcto para superficies horizontales.

### 3.3 Ganancia de antena — APROXIMACIÓN ACEPTABLE ✓

El patrón de radiación es gaussiano separable con semiancho a −3 dB definido por las aperturas configuradas (50° en elevación, 14° en azimut). El boresight se orienta siempre hacia el origen [0,0,0], lo cual es una aproximación válida si el área de interés está centrada en ese punto.

### 3.4 Estrategia multi-start — RAZONABLE ✓

Los 5 puntos iniciales cubren: (1) opuesto a los targets, (2–3) perpendiculares, (4–5) región guiada por el ángulo de Brewster con jitter aleatorio. Esta estrategia reduce el riesgo de converger a mínimos locales, aunque no lo elimina completamente dado que `fmincon` es un optimizador local.

### 3.5 Ecuación de radar — CORREGIDO ✓

El denominador corregido usa las distancias totales de cada camino refractado: (R₁+R₂)² · (R₃+R₄)². Ver §4.1 para el análisis del error original.

### 3.6 Frecuencia portadora — CORREGIDO ✓

Cambiada a 400 MHz (P-band), compatible con penetración subsuperficial a 2–10 m en suelo seco. Ver §4.2 para el análisis detallado.

---

## 4. Críticas identificadas y estado de corrección

### 4.1 [CORREGIDO] Inconsistencia dimensional en la ecuación de radar

**Error original**: La implementación usaba cuatro distancias independientes al cuadrado en el denominador:

$$\text{original (incorrecto):}\quad \frac{1}{(4\pi)^3 \cdot R_1^2\, R_2^2\, R_3^2\, R_4^2}$$

Esto equivale a modelar la interfaz aire–suelo como una **fuente puntual secundaria** que reemite omnidireccionalmente en cada tramo, lo que no corresponde a la física de la refracción dieléctrica. El error introduce una penalización adicional de spreading que crece como el producto cruzado de las distancias parciales, sobreestimando la atenuación por un factor:

$$\text{Error multiplicativo} = \frac{R_T^2\,R_R^2}{R_1^2\,R_2^2\,R_3^2\,R_4^2} \cdot R_T^2\,R_R^2 = \left(\frac{R_1+R_2}{R_1\,R_2}\right)^2 \left(\frac{R_3+R_4}{R_3\,R_4}\right)^2$$

**Corrección aplicada** en `objective_function.m`:

```matlab
% Antes (incorrecto):
((4*pi)^3 * R_Tx_int1^2 * R_int1_Tg^2 * R_Tg_int2^2 * R_int2_Rx^2)

% Después (correcto):
R_T = R_Tx_int1 + R_int1_Tg;
R_R = R_Tg_int2 + R_int2_Rx;
((4*pi)^3 * R_T^2 * R_R^2)
```

**Nota**: El error era sistemático en todas las evaluaciones de la función objetivo, de modo que el argmax (posición óptima del Rx) seguía siendo representativo de la geometría óptima. Sin embargo, los **valores de potencia en dBm** reportados eran incorrectos en escala absoluta y la corrección es indispensable para cualquier presupuesto de enlace o análisis de SNR.

### 4.2 [CORREGIDO] Incompatibilidad frecuencia–penetración en suelo

**Error original**: La frecuencia portadora era 10 GHz (banda X), λ ≈ 3 cm.

La profundidad de penetración efectiva (distancia a la que la potencia cae a 1/e²) en un medio dieléctrico con pérdidas es:

$$\delta_p = \frac{1}{\alpha}, \qquad \alpha = \omega\sqrt{\frac{\mu\varepsilon'}{2}\left(\sqrt{1 + \left(\frac{\varepsilon''}{\varepsilon'}\right)^2} - 1\right)}$$

Para suelo seco (εᵣ ≈ 4, tanδ ≈ 0.01):

| Frecuencia | λ | Atenuación (suelo seco) | δ_p |
|---|---|---|---|
| 10 GHz (X) | 3 cm | ~25 dB/m | ~0.17 m |
| 400 MHz (P) | 75 cm | ~0.5 dB/m | ~8 m |

A 10 GHz, los targets a 2–10 m de profundidad sufrirían pérdidas de **50–500 dB** solo por absorción dieléctrica en el suelo, haciendo la detección físicamente inviable.

**Corrección aplicada**: Frecuencia portadora cambiada a **400 MHz** (banda P) en `radarTx_espiral_plano_voo.json` y `radarRx_espiral_plano_voo.json`.

**Justificación de 400 MHz**:
- λ = 0.75 m → resolución en rango con B = 50 MHz: ΔR = c/(2B) ≈ 3 m
- Penetración en suelo seco: δ_p ~ 6–10 m, compatible con targets a 2–10 m
- El ángulo de Brewster θ_B = arctan(2) ≈ 63.43° no cambia (independiente de f para lossless)
- Sistemas reales de P-band GPR-SAR operan en este rango (e.g., BIOMASS ESA ~435 MHz)

**Impacto en la ecuación de Friis**: el factor λ² aumenta de (0.03)² = 9×10⁻⁴ m² a (0.75)² = 0.5625 m², lo que representa una ganancia de ~28 dB en la potencia recibida, compensa parcialmente las mayores pérdidas de spreading a las distancias del sistema.

### 4.3 [CRÍTICA ACTIVA] Modelo de suelo lossless

El índice de refracción n₂ = 2 es puramente real, asumiendo un dieléctrico sin pérdidas. En suelo real:

$$n_2 = n' - jn'' = \sqrt{\varepsilon_r\left(1 - j\tan\delta\right)}$$

La parte imaginaria introduce:
1. Un **pseudo-ángulo de Brewster** donde |r_TM|_min > 0 (el cero exacto desaparece).
2. Atenuación exponencial de amplitud exp(−αR₂) en el subsuelo que debe multiplicar las distancias subterráneas.

Para una primera aproximación en el estudio de geometría bistática, el modelo lossless es aceptable. Sin embargo, para resultados cuantitativos de SNR o presupuesto de enlace debe incorporarse tanδ y la atenuación por tramo.

### 4.4 [CRÍTICA ACTIVA] Restricción de altitud del receptor

El receptor se fuerza a la misma altitud que el transmisor (`Rx_z = PzT(tx_idx)`), reduciendo la optimización a un problema 2D (x,y). Consecuencias:
- No se explora el espacio de altitudes óptimas del receptor.
- No hay restricciones de velocidad ni aceleración del vehículo receptor.
- La interpolación PCHIP garantiza suavidad matemática pero no trayectoria dinámicamente factible.

### 4.5 [OBSERVACIÓN] Reproducibilidad

Los puntos candidatos 4 y 5 en `chooseCandidateRxPos.m` usan `randn` y `rand` sin semilla fija. Para reproducibilidad científica, agregar `rng(semilla_fija)` al inicio de `run_plano_de_voo.m`.

### 4.6 [OBSERVACIÓN] Polarización implícita

El script asume polarización TM (p-polarización) pero los parámetros JSON no lo declaran explícitamente. Si la antena emite en TE, el ángulo de Brewster no existe y la premisa del sistema colapsa. Debe documentarse como hipótesis de diseño.

---

## 5. Estimación numérica de la geometría

### 5.1 Verificación del ángulo de incidencia en la espiral

Para un Tx a altitud h = 100 m y radio espiral r ≈ 160 m, el ángulo de incidencia sobre la superficie (plano z = 0) desde la proyección vertical es:

$$\theta_{inc} = \arctan\!\left(\frac{r}{h}\right) = \arctan\!\left(\frac{160}{100}\right) \approx 58°$$

Diferencia con el ángulo de Brewster: |63.43° − 58°| ≈ **5.4°**. La espiral no alcanza exactamente θ_B. Para alcanzar θ_B a h = 100 m se requeriría:

$$r_{Brewster} = h \cdot \tan(\theta_B) = 100 \cdot 2 = 200 \text{ m}$$

El radio actual (147.5–172.5 m) está ~15% por debajo de lo necesario en la trayectoria del Tx. La optimización del Rx compensa este déficit geométrico buscando posiciones que acerquen el ángulo efectivo del tramo de retorno a θ_B.

### 5.2 Distancia de Brewster para el receptor

`buildRxSearchGeometry` calcula la distancia en planta desde el centro de los targets hasta el Rx óptimo bajo la hipótesis de θ_B:

$$d_{Brewster} = \frac{h_{Rx} + |z_{tg}^{medio}|}{\tan\theta_B} = \frac{h_{Rx} + 6}{2}$$

Para h_Rx = 100 m y profundidad media de 6 m: **d_Brewster = 53 m**. Esta es la distancia horizontal esperada del receptor al centro del área de targets si la geometría fuera puramente de Brewster.

### 5.3 Presupuesto de enlace estimado (post-corrección)

Con los parámetros corregidos y una geometría aproximada (R_T = R_R = 160 m, T₁ = T₂ = 1 en θ_B):

$$P_r \approx \frac{10^5 \cdot 4 \cdot 4 \cdot 1 \cdot (0.75)^2}{(4\pi)^3 \cdot (160)^2 \cdot (160)^2} \approx \frac{10^5 \cdot 16 \cdot 0.5625}{(4\pi)^3 \cdot 6.55 \times 10^9}$$

$$P_r \approx \frac{9 \times 10^5}{(1984) \cdot 6.55 \times 10^9} \approx 6.9 \times 10^{-11} \text{ W} \approx -71.6 \text{ dBm}$$

Este valor es una cota superior (suelo lossless, T = 1 en ambas interfaces, σ = 1 m²) y sirve como referencia para evaluar el margen de SNR del sistema.

---

## 6. Síntesis y veredicto

| Aspecto | Estado | Severidad |
|---|---|---|
| Coeficientes de Fresnel TM | Correcto ✓ | — |
| Transmitancia de potencia (T + R = 1) | Correcto ✓ | — |
| Principio de Fermat (punto de refracción) | Correcto ✓ | — |
| Ángulo de incidencia respecto a normal | Correcto ✓ | — |
| Ecuación de radar — denominador R₁²R₂²R₃²R₄² | **Corregido** ✓ | Alta |
| Frecuencia 10 GHz para penetración subsuperficial | **Corregido** ✓ | Alta |
| Modelo de suelo lossless | Aproximación activa | Media |
| Restricción altitud Rx = altitud Tx | Simplificación activa | Media |
| Polarización TM no declarada en parámetros | Hipótesis implícita | Baja |
| Reproducibilidad (semilla aleatoria) | Pendiente | Baja |

**Conclusión**: Tras las correcciones aplicadas, la simulación es **física y matemáticamente coherente** en su estructura principal. El principio de explotar el ángulo de Brewster para maximizar la transmitancia TM en la interfaz aire–suelo está correctamente formulado, la ecuación de radar ahora usa las distancias totales de propagación correctas, y la frecuencia de 400 MHz es compatible con la penetración subsuperficial a las profundidades modeladas. Los valores de potencia reportados son cuantitativamente válidos bajo la hipótesis de suelo lossless, que constituye el límite superior físicamente realizable del sistema.
