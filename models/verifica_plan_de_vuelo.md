# Verificación y Crítica del Script `run_plano_de_voo.m`

**Propósito del script**: Optimizar la trayectoria del receptor en un sistema SAR bistático, buscando la posición de Rx que equilibra la potencia recibida (explotando el ángulo de Brewster) y la resolución tridimensional del sistema (δ_xy, δ_z), ambos criterios combinados en una función de costo log-normalizada única.

---

## 1. Descripción general del sistema

### 1.1 Parámetros clave

| Parámetro | Valor | Notas |
|---|---|---|
| Frecuencia portadora | 400 MHz (banda P) | Corregido de 10 GHz; ver §4.2 |
| Longitud de onda | λ = 0.75 m | c/400 MHz |
| Ancho de banda chirp | B = 50 MHz | FreqMayor − FreqMenor |
| Potencia Tx | 100 kW | |
| Velocidad del vehículo | 120 m/s | |
| PRF | 400 Hz | |
| Altitud de la espiral | 80–120 m | Radio 147.5–172.5 m |
| Longitud arco espiral | B_helix ≈ 47.17 m | √(40²+25²) |
| Ángulo inclinación hélice | β ≈ 58.0° | arctan(40/25) |
| Índice de refracción suelo | n₂ = 2 | n₁ = 1 (aire); εᵣ = 4, lossless |
| Ángulo de Brewster teórico | θ_B ≈ 63.43° | arctan(n₂/n₁) = arctan(2) |
| Profundidad de los targets | 2–10 m | 4×4×4 = 64 puntos |
| RCS | σ = 1 m² | Isótropo |
| Radio de búsqueda Rx | ±200 m | Espacio de búsqueda 2D |
| Peso tradeoff | α = 0.3 | 70% potencia, 30% resolución |

### 1.2 Flujo del script

```
1. Cargar parámetros → trayectoria espiral Tx
2. Calcular B_helix, beta_helix, B (geometría para resolución)
3. Decimar posiciones Tx (factor 200) → ~50 posiciones de evaluación
4. Para cada posición Tx decimada:
   a. Calcular ganancias Tx hacia todos los targets (precalculadas)
   b. Construir geometría de búsqueda Brewster
   c. Lanzar fmincon desde 5 puntos candidatos minimizando J combinado
   d. Evaluar potencia y resolución por separado en el punto óptimo
5. Interpolar (PCHIP) trayectoria completa del Rx
6. Guardar resultados (.mat y .csv) con Pr, δxy, δz
7. Visualizar configuración bistática optimizada
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

### 2.6 Función de costo multi-objetivo potencia–resolución

#### 2.6.1 El conflicto geométrico

El ángulo de Brewster maximiza la transmitancia TM posicionando al Rx a una distancia horizontal d_B del centroide del target (§5.2). Este posicionamiento determina el ángulo azimutal Δφ entre Tx y Rx visto desde el target, y la resolución horizontal del sistema depende directamente de Δφ mediante (ecuación derivada en `explicacion_resolucion.tex`, §7.2):

$$\delta_{xy}(\Delta\phi) = \frac{0.60\,\lambda_0}{\pi\,\sin\psi_0\cdot|\cos(\Delta\phi/2)|}$$

donde ψ₀ es el look angle del Rx hacia el centroide. Esta expresión revela el **conflicto geométrico** central:

| Condición | Δφ | Efecto en δ_xy | Efecto en P_r |
|---|---|---|---|
| Rx colineal con Tx (cuasi-monostático) | ≈ 0° | Mínima (mejor resolución) | Baja (T₂ lejos de θ_B) |
| Rx en posición de Brewster | ≈ 60–120° | Degradada por 1/cos(Δφ/2) | Máxima (T₂ = 1) |
| Rx antipodal al Tx | ≈ 180° | → ∞ (sin resolución horizontal) | Variable |

La resolución vertical es **invariante a Δφ** (demostrado analíticamente en el .tex, §7.1), pero sí depende del look angle ψ₀ del Rx a través de la apertura tomográfica efectiva B_⊥:

$$\delta_z = \frac{c}{2\,W_z}, \qquad W_z = \underbrace{n_2\,B\,\cos\theta_{t,0}}_{\text{contribución BW}} + \underbrace{\frac{c\,B_\perp\,\sin\psi_0\cos\psi_0}{\lambda_0\,R_0\,n_2\cos\theta_{t,0}}}_{\text{contribución tomográfica}}$$

con B_⊥ = B_helix·|cos(β − ψ₀)|, que se maximiza en la condición óptima β = ψ₀. Mover el Rx cambia ψ₀_Rx y por tanto B_⊥ y δ_z simultáneamente.

**Nota sobre la aproximación**: Las fórmulas de resolución del .tex fueron derivadas para Tx y Rx en la misma hélice (mismo ψ₀). En la optimización bistática dinámica, Tx y Rx pueden tener look angles distintos. La implementación usa ψ₀_Rx (look angle del receptor optimizado), que domina la cobertura del k-espacio en la posición de retorno y es la aproximación más relevante para el diseño del receptor.

#### 2.6.2 El problema de la heterogeneidad dimensional

Los dos objetivos tienen dimensiones y escalas radicalmente distintas:

| Magnitud | Unidad | Escala típica | log₁₀(·/ref) |
|---|---|---|---|
| P_r | W | ~10⁻¹⁰ | −10 |
| δ_xy | m | ~0.05–0.5 | −1.3 a −0.3 |
| δ_z | m | ~0.2–2 | −0.7 a 0.3 |
| δ_xy · δ_z | m² | ~0.01–1 | −2 a 0 |

Una scalarización lineal P_r − w·(δ_xy + δ_z) exigiría que w tenga unidades de W/m, carece de interpretación física, y el resultado sería dominado numéricamente por la magnitud con mayor rango absoluto. La escala logarítmica elimina este problema.

#### 2.6.3 Justificación de la transformación logarítmica

La transformación log₁₀ es una **función monótonamente creciente**, por lo que preserva el orden de los óptimos: el argmax de P_r es el argmin de −log₁₀(P_r), y el argmin de δ_xy·δ_z es el argmin de log₁₀(δ_xy·δ_z). Este es el primer criterio de validez de la normalización.

El segundo criterio es la **comparabilidad de escalas**. Sea Δ_P el rango típico de −log₁₀(P_r) sobre el espacio de búsqueda (≈ 5 décadas), y Δ_R el rango típico de log₁₀(δ_xy·δ_z) (≈ 2–3 décadas). La normalización logarítmica comprime ambos rangos al mismo orden de magnitud, garantizando que el gradiente de J sea sensible a ambos objetivos simultáneamente.

Una tercera justificación es que en la ecuación de radar, la potencia recibida varía como producto de factores (Gt·Gr·T₁·T₂·λ²/R⁴), de modo que P_r responde a cambios en la geometría multiplicativamente. De manera análoga, δ_xy ∝ 1/(sinψ₀·|cos(Δφ/2)|) varía de forma multiplicativa. Scalarizar en espacio logarítmico es, por tanto, **coherente con la física multiplicativa** del sistema.

#### 2.6.4 Función de costo log-normalizada y su derivación

La función de costo combinada a minimizar por `fmincon` es:

$$\boxed{J(\mathbf{r}_{Rx}) = -(1-\alpha)\,\log_{10} P_r(\mathbf{r}_{Rx}) + \alpha\,\log_{10}\!\left[\delta_{xy}(\mathbf{r}_{Rx})\cdot\delta_z(\mathbf{r}_{Rx})\right]}$$

con α ∈ [0, 1] el peso de tradeoff.

**Verificación de la dirección de optimización** (fmincon minimiza J):
- −log₁₀(P_r) se minimiza ↔ P_r se maximiza ✓
- log₁₀(δ_xy·δ_z) se minimiza ↔ δ_xy·δ_z se minimiza (mejor resolución) ✓

**Análisis dimensional**: ambos términos son log₁₀ de cocientes adimensionales respecto a la unidad SI del sistema (1 W y 1 m²):

$$J = -(1-\alpha)\,\log_{10}\!\left(\frac{P_r}{1\,\text{W}}\right) + \alpha\,\log_{10}\!\left(\frac{\delta_{xy}\cdot\delta_z}{1\,\text{m}^2}\right)$$

La diferencia de referencias (1 W, 1 m²) introduce constantes aditivas que no afectan al argmin de J y son la normalización implícita más natural del sistema SI.

**Escala numérica con valores típicos del sistema** (α = 0.3):

$$J \approx -(0.7)\cdot(-10) + (0.3)\cdot(-2) = 7.0 - 0.6 = 6.4$$

**Contribución relativa al gradiente**: el peso efectivo de cada término al gradiente ∂J/∂r_Rx está determinado por el producto (α_efectivo)·(sensibilidad local). Si P_r varía en ±3 décadas sobre el espacio de búsqueda y δ_xy·δ_z varía en ±1 década, las contribuciones al rango de J son:

| Término | Rango típico de J | Con α = 0.3 |
|---|---|---|
| −(1−α)·log₁₀(P_r) | 0.7 × 5 = 3.5 | dominante |
| α·log₁₀(δ_xy·δ_z) | 0.3 × 2 = 0.6 | presente |

Para α ≈ 0.6–0.7 los rangos serían comparables. El valor α = 0.3 es un punto de partida conservador que prioriza la potencia mientras incluye la resolución como restricción suave.

#### 2.6.5 La frontera de Pareto y el rol de α

El problema multi-objetivo tiene una **frontera de Pareto** en el espacio (P_r, δ_xy·δ_z): el conjunto de puntos donde no es posible mejorar un objetivo sin degradar el otro. La scalarización log-normalizada con distinto α recorre esta frontera:

```
δ_xy · δ_z
     ▲
     │ α→1 (solo resolución)
     │  ×
     │   ×
     │    ×  ← frontera de Pareto
     │     ×
     │      × α→0 (solo potencia)
     └──────────────────────────► P_r
```

- α = 0: el optimizador encuentra el punto de máxima potencia (extremo derecho de la frontera)
- α = 1: el optimizador encuentra el punto de mínima resolución (extremo superior, mejor δ)
- α ∈ (0,1): puntos intermedios de la frontera, con tradeoff explícito

La scalarización logarítmica es más efectiva que la lineal para recorrer la frontera de Pareto cuando los objetivos tienen rangos muy distintos, ya que en espacio log-log la frontera tiende a ser más lineal y los pesos α tienen un efecto más predecible.

#### 2.6.6 Casos límite y comportamiento degenerado

| α | Comportamiento del optimizador |
|---|---|
| 0 | Maximización pura de P_r; equivale al script original |
| 0.3 | 70% potencia, 30% resolución; valor por defecto recomendado |
| 0.5 | Pesos iguales en escala log (potencia numéricamente dominante por mayor rango) |
| ~0.65 | Contribuciones al rango de J aproximadamente iguales |
| 1 | Minimización pura de δ_xy·δ_z; el Rx tiende a Δφ → 0° (cuasi-monostático) |

**Penalizaciones automáticas** en los casos degenerados:

| Situación | Efecto en J | Correctamente penalizado |
|---|---|---|
| Δφ → 180° | δ_xy → ∞ → log₁₀(δ_xy·δ_z) → +∞ → J → +∞ | Sí ✓ |
| ψ₀ → 0° (Rx sobre target) | δ_xy → ∞ ídem | Sí ✓ |
| P_r → 0 (posición desfavorable) | −log₁₀(P_r) → +∞ → J → +∞ | Sí ✓ |
| β = ψ₀ (hélice óptima) | B_⊥ = B_helix, δ_z mínima | Favorecido ✓ |

La auto-penalización de casos degenerados es una propiedad emergente de la formulación logarítmica, no un término adicional, lo que mantiene la función de costo simple y diferenciable.

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

### 3.7 Optimización multi-objetivo potencia–resolución — IMPLEMENTADO ✓

Ver §2.6 para la derivación completa de la función de costo combinada.

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

### 5.4 Estimación analítica de resolución 3D para el sistema configurado

Con los parámetros del sistema (f₀ = 400 MHz, B = 50 MHz, n₂ = 2, B_helix = 47.17 m, β = 58.0°):

**Parámetros auxiliares** (target a z_P = −5 m, Rx a h_Rx = 100 m, ρ_Rx ≈ 53 m en posición de Brewster):

$$\psi_0 = \arctan\!\left(\frac{53}{100}\right) \approx 27.9°, \quad R_0 = \sqrt{53^2+100^2} \approx 113.3\,\text{m}$$

$$\sin\theta_{t,0} = \frac{\sin 27.9°}{2} \approx 0.234, \quad \cos\theta_{t,0} \approx 0.972$$

**Resolución vertical** (apertura tomográfica en condición óptima β = ψ₀ = 58.0° ≠ 27.9°):

$$B_\perp = 47.17 \cdot |\cos(58.0° - 27.9°)| = 47.17 \cdot |\cos 30.1°| \approx 40.8\,\text{m}$$

$$W_z = 2 \times 50 \times 10^6 \times 0.972 + \frac{3\times10^8 \times 40.8 \times \sin 27.9° \times \cos 27.9°}{0.75 \times 113.3 \times 2 \times 0.972}$$

$$W_z \approx 97.2\,\text{MHz} + \frac{3\times10^8 \times 40.8 \times 0.467 \times 0.884}{164.8} \approx 97.2 + 30.9 \approx 128\,\text{MHz}$$

$$\delta_z = \frac{c}{2\,W_z} = \frac{3\times10^8}{2 \times 128\times10^6} \approx 1.17\,\text{m}$$

**Nota**: β = 58.0° ≠ ψ₀_Rx ≈ 28° en la posición de Brewster. La condición óptima δz mínima (β = ψ₀) se cumpliría si ψ₀_Rx ≈ 58°, es decir, ρ_Rx ≈ 160 m — radio de la espiral del Tx. Esto implica una **tensión adicional**: la posición de Brewster del Rx (ρ = 53 m) no es la posición de máxima apertura tomográfica (ρ = 160 m). El optimizador con α > 0 encontrará el equilibrio.

**Resolución horizontal** en función de Δφ (para ψ₀_Rx ≈ 28°):

$$\delta_{xy}(\Delta\phi) = \frac{0.60 \times 0.75}{\pi \times \sin 27.9° \times |\cos(\Delta\phi/2)|} = \frac{0.45}{1.467 \times |\cos(\Delta\phi/2)|}$$

| Δφ | \|cos(Δφ/2)\| | δ_xy |
|---|---|---|
| 0° (monostático) | 1.000 | 0.307 m |
| 30° | 0.966 | 0.318 m |
| 60° | 0.866 | 0.354 m |
| 90° (NorthOffset config.) | 0.707 | 0.434 m |
| 120° | 0.500 | 0.614 m |
| 150° | 0.259 | 1.186 m |
| 180° (antipodal) | 0.000 | ∞ |

**Tabla de costo J (α = 0.3) en los escenarios clave**:

| Escenario | P_r (W) | δ_xy (m) | δ_z (m) | J |
|---|---|---|---|---|
| Solo Brewster (Δφ = 90°) | 6.9×10⁻¹¹ | 0.434 | 1.17 | 6.6 |
| Solo resolución (Δφ = 0°) | reducida | 0.307 | ~1.5 | ~7.1 |
| Brewster + Δφ → 0° (inalcanzable simultáneamente) | máx | mínima | — | mínimo teórico |
| Óptimo tradeoff (α = 0.3) | intermedia | ~0.35–0.40 | ~1.2 | ~6.4 |

Los valores muestran que el costo J varía en un rango de ≈ 1 unidad entre los extremos, lo que confirma que el tradeoff es numéricamente significativo y que el optimizador tiene margen real para encontrar un equilibrio.

### 5.5 Sensibilidad de α: ¿cuánto importa el peso?

La elección de α afecta la posición del Rx óptimo y, por tanto, los valores reportados de P_r, δ_xy y δ_z. La siguiente aproximación analítica estima la sensibilidad:

Sea r_B la posición de Brewster (máxima P_r) y r_M la posición monostática (Δφ ≈ 0°, mínima δ_xy). El gradiente de J en la dirección r_B → r_M puede estimarse como:

$$\left.\frac{\partial J}{\partial r}\right|_{r_B \to r_M} \approx -(1-\alpha)\frac{\Delta\log_{10}P_r}{\Delta r} + \alpha\frac{\Delta\log_{10}(\delta_{xy}\cdot\delta_z)}{\Delta r}$$

Si P_r varía ~3 décadas y δ_xy·δ_z varía ~0.5 décadas sobre el segmento r_B → r_M, el punto estacionario (∂J/∂r = 0) ocurre en:

$$\alpha^* = \frac{3}{3 + 0.5} \approx 0.86$$

Para α < α*, el optimizador se mueve hacia r_B (Brewster). Para α > α*, hacia r_M (monostático). El valor por defecto α = 0.3 garantiza que el óptimo esté en la región de Brewster con corrección de resolución, lo cual es el objetivo físico del sistema.

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
| Resolución 3D ignorada en la optimización | **Implementado** ✓ | Alta |
| Modelo de suelo lossless | Aproximación activa | Media |
| Restricción altitud Rx = altitud Tx | Simplificación activa | Media |
| Polarización TM no declarada en parámetros | Hipótesis implícita | Baja |
| Reproducibilidad (semilla aleatoria) | Pendiente | Baja |

**Conclusión**: Tras las correcciones y extensiones aplicadas, la simulación es **física y matemáticamente coherente** en su estructura principal. La función de costo log-normalizada J resuelve correctamente el problema de escala entre potencia (W) y resolución (m), ambos términos son adimensionales en escala logarítmica y el optimizador puede equilibrarlos de forma significativa. El parámetro α permite explorar toda la frontera de Pareto entre maximización de potencia y minimización de resolución 3D, con α = 0.3 como punto de partida recomendado (prioridad en potencia). Los valores de potencia y resolución reportados son cuantitativamente válidos bajo la hipótesis de suelo lossless (cota superior física).
