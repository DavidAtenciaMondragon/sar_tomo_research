# Modelo Conceptual Unificado — SAR/GPR: Resolución, Geometría y Fase

**Generado a partir del análisis comparativo de 22 papers/notas de la carpeta `notes/`**  
**Fecha:** 2026-05-31

---

## Índice

1. [Ecuaciones de resolución más comunes](#1-ecuaciones-de-resolución-más-comunes)
2. [Diferencias monoestático vs. biestático](#2-diferencias-monoestático-vs-biestático)
3. [Modelado de la geometría por tipo de sistema](#3-modelado-de-la-geometría-por-tipo-de-sistema)
4. [Elementos fundamentales comunes](#4-elementos-fundamentales-comunes)
5. [Inconsistencias y diferencias importantes entre papers](#5-inconsistencias-y-diferencias-importantes-entre-papers)
6. [Modelo conceptual unificado](#6-modelo-conceptual-unificado)

---

## 1. Ecuaciones de Resolución Más Comunes

### 1.1 Resolución en Rango — Universal

Aparece en todos los papers SAR y GPR sin excepción:

$$\boxed{\delta_r = \frac{c}{2B}}$$

| Paper | Notación | Observación |
|-------|----------|-------------|
| Moreira et al. 2013 | $\delta_r = c_0/(2B_r)$ | $B_r$ = ancho de banda chirp |
| Gorham & Moore 2010 | $\delta_r = c/(2B)$ con $B=(K{-}1)\Delta f$ | Implementación MATLAB |
| Drone-DInSAR, García-F. 2019, 2024 | $\delta_r = c/(2W)$ | W = BW efectivo |
| Banda et al. 2015 | $\delta_r = c/(2W)$ → 1.8 m a 85 MHz | TomoSAR P-band |
| GPR buried objects 2019 | $\Delta R = v_{prop}/(2\Delta BW)$ | $v_{prop} = c/\sqrt{\varepsilon_r}$ |
| UAV GPR 2024 | $\delta_z \approx c/(2\sqrt{\varepsilon_r}\,BW)$ | Resolución en profundidad |

**Unificación:** $\delta_r = v/(2B)$, donde $v = c/\sqrt{\varepsilon_r}$ es la velocidad en el medio.

---

### 1.2 Resolución en Azimut / Rango Cruzado — SAR Lineal

$$\delta_a = \frac{d_a}{2}$$

Aparece en Moreira et al. 2013 y Drone-DInSAR. Es independiente del rango, válida solo en modo stripmap. La resolución tomográfica (Munson 1983, formulación en el dominio k) es:

$$\delta_x \approx \frac{\pi c}{2\alpha T}, \qquad \delta_y \approx \frac{\pi c}{2\omega_0 \sin\theta_M}$$

Gorham & Moore 2010 da la forma equivalente para spotlight:

$$\delta_x = \frac{\lambda_c}{2\theta_a}$$

donde $\theta_a = (N_p-1)\Delta\theta$ es el ángulo de apertura total.

---

### 1.3 Resolución en Plano — SAR Circular/Espiral

$$\delta_{xy} = \frac{1.12\lambda}{2\pi\sin\psi}$$

De Ishimaru et al. (1998), citado por Góes 2022 y aplicado también al SAR espiral. Válida para blancos en el eje del círculo ($x_t = 0$). Para ángulo $\psi = 45°$: $\delta_{xy} \approx \lambda/4$.

---

### 1.4 Resolución Vertical — SAR Tomografía y Espiral

**Caso tomografía clásica (apertura vertical):**

$$\delta_{LOS_\perp} = \frac{\lambda R_0}{2B}$$

**Caso Circular SAR (una sola vuelta):**

$$\delta_z = \sqrt{\frac{\ln(2)}{\pi}\frac{c}{W\cos\psi}}$$

**Caso Spiral/Multi-Circular SAR (apertura tomográfica $B_\perp$):**

$$\boxed{\delta_z = \sqrt{\frac{\ln(2)}{\pi}\frac{c}{W_z}}, \quad W_z = W\cos\psi_0 + \frac{cB_\perp}{\lambda R_0}\sin\psi_0}$$

Esta es la ecuación central de Góes 2022 (revisada respecto a Ponce et al. 2016).

---

### 1.5 Resolución Bistática (Horne & Yates, Ding & Munson)

$$\rho_u = \frac{c}{f_0 K}, \qquad \rho_v = \frac{c}{2\Delta f \cos(\varphi/2)}$$

donde $K^2 = (\Delta\alpha)^2 + (\Delta\beta)^2 + 2\Delta\alpha\Delta\beta\cos\varphi$ es el factor de apertura biestático. La resolución en rango biestático (dirección $v$) depende del ángulo biestático $\varphi$ a través del factor $\cos(\varphi/2) = |\cos\beta|$ (usando la notación de semi-ángulo $\beta$ de Arikan & Munson).

---

### 1.6 Resolución SAR Universal (dominio número de onda)

$$\Delta A_{SAR} = \frac{\lambda_c}{2(\vartheta_2-\vartheta_1)} \cdot \frac{c}{2B}$$

Expresión de Ulander et al. 2003, que unifica rango y acimut: la resolución es el producto de las resoluciones en frecuencia y ángulo, vista en el dominio de número de onda. Es la base teórica para todos los algoritmos backprojection.

---

### 1.7 Resolución GPR-SAR (Down-Looking, García-Fernández et al.)

$$\delta_u \approx \lambda_c \frac{\sqrt{L^2/4+h^2}}{2L}, \qquad \delta_z \approx \frac{c}{2\sqrt{\varepsilon_r}\,BW}$$

Con $L$ = longitud de la máscara SAR y $h$ = altura de vuelo.

---

## 2. Diferencias Monoestático vs. Biestático

### 2.1 Definición de Distancia

| Sistema | Distancia a píxel $\mathbf{p}$ | Factor de escala k-espacio |
|---------|-------------------------------|---------------------------|
| **Monoestático** | $R = \|\mathbf{r}_{ant} - \mathbf{p}\|$ | $2/c$ |
| **Biestático** | $R_{total} = \|\mathbf{r}_t - \mathbf{p}\| + \|\mathbf{r}_r - \mathbf{p}\|$ | $w_{tr}/c = 2\cos\beta/c$ |
| **GPR subsuperficial** | $R_{eff} = R_{aire} + \sqrt{\varepsilon_r}R_{sub}$ | $2/c$ (sobre camino óptico) |

### 2.2 Fase de la Señal

| Sistema | Fase | Papers |
|---------|------|--------|
| Monoestático | $\Phi = \frac{4\pi}{\lambda}R = \frac{4\pi f}{c}R$ | SAR Tutorial, Basics BP, Toolbox MATLAB |
| Biestático | $\Phi = \frac{2\pi}{\lambda}(R_t+R_r)$ | Arikan, Ding, Bistatic SAR |
| Spiral SAR (interferometría) | $\Delta\Phi = \frac{4\pi}{\lambda}\Delta R$ | Góes 2022, SAR Tutorial |
| DInSAR deformación | $\phi_{def} = \frac{4\pi}{\lambda}d_{LOS}$ | Drone DInSAR |
| GPR con refracción | $\Phi = \frac{4\pi}{\lambda}(R_{aire} + \sqrt{\varepsilon_r}R_{sub})$ | Imoc BP, García-F. 2019 |

### 2.3 Geometría de Iso-Superficies

| Sistema | Iso-rango | Iso-Doppler | Dominio Fourier |
|---------|-----------|-------------|-----------------|
| Monoestático lineal | Semicírculos | Líneas paralelas | Sector anular |
| Monoestático circular | Círculos concéntricos | — | Anillo en $k_x, k_y$ |
| Biestático | **Elipses** (focos en Tx y Rx) | Hipérbolas | Arcos elípticos |
| GPR downward-looking | Hemisferios | — | Semiesfera en $k_x,k_y,k_z$ |

### 2.4 Factor Biestático $w_{tr}$ (Arikan & Munson 1988)

$$w_{tr} = 2|\cos\beta|, \quad \beta = \text{semi-ángulo biestático}$$

- Monoestático: $\beta = 0°$ → $w_{tr} = 2$
- Biestático genérico: $0 < \beta < 90°$ → $0 < w_{tr} < 2$
- Consecuencia directa: el biestático adquiere **menos** soporte en el k-espacio → menor resolución

### 2.5 Problemas Exclusivos del Biestático

1. **Sincronización TX/RX:** requiere <100 ns temporal y osciladores atómicos en frecuencia (Horne & Yates).
2. **Ruido de fase acumulado:** los osciladores de TX y RX son independientes → ruido efectivo = suma de los dos (≈3 dB mayor que monoestático).
3. **Ejes de imagen no ortogonales:** los ejes $u$ y $v$ de resolución biestática no son perpendiculares en general.
4. **Algoritmos adaptados:** el PFA monoestático requiere reproyección al ángulo bisectriz (Horne & Yates); el Fast BP requiere interpolación en bandas de paso en lugar de paso bajo (Ding & Munson).

---

## 3. Modelado de la Geometría por Tipo de Sistema

### 3.1 SAR Lineal (Stripmap/Spotlight)

**Marco de referencia:** origen en el Scene Reference Point (SRP). Ejes $(x,y,z)$ con $x$ = azimut, $y$ = rango cruzado, $z$ = altura.

**Historia de rango** (Moreira et al., Doerry et al.):
$$r(t) = \sqrt{r_0^2 + (vt)^2} \approx r_0 + \frac{(vt)^2}{2r_0}$$

**Señal de acimut** como chirp FM cuadrático, con frecuencia Doppler:
$$f_D = -\frac{2v^2 t}{\lambda r_0}$$

**Suposición clave:** stop-and-go + trayectoria rectilínea uniforme.

---

### 3.2 SAR Spotlight — Formulación Tomográfica

**Marco de referencia:** origen en el centro del parche. El radar observa desde ángulo $\theta$ que varía con el tiempo.

**El teorema de proyección-rebanada** (Munson et al. 1983):
$$C_\theta(t) \approx \frac{A}{2} P_\theta\!\left[\frac{2}{c}(\omega_0 + 2\alpha(t-\tau_0))\right]$$

La señal demodulada = TF 1D de la proyección de $g(x,y)$ → el SAR adquiere rebanadas del k-espacio 2D.

**Datos en el k-espacio:** sector anular $[X_1, X_2] \times [-\theta_M, \theta_M]$.

**Extensión al biestático** (Arikan & Munson 1988): misma formulación con factor $w_{tr}/c$ en lugar de $2/c$.

---

### 3.3 SAR Circular / Espiral (Góes 2022)

**Marco de referencia:** cilíndrico $(\rho, \alpha, z)$ con origen en el centro de la escena.

**Parámetros esenciales:**

| Parámetro | Símbolo | Relación |
|-----------|---------|----------|
| Ángulo de incidencia medio | $\psi_0$ | $\psi_0 = \tan^{-1}(\rho_0/z_0)$ |
| Distancia media | $R_0$ | $R_0 = \sqrt{z_0^2+\rho_0^2}$ |
| Apertura tomográfica efectiva | $B_\perp$ | $B_\perp = B|\cos(\beta-\psi_0)|$ |
| Ángulo de inclinación | $\beta$ | $\beta=90°$ cilíndrica, $\beta=\psi_0$ cónica óptima |

**Trayectoria espiral:** $\rho(t)$, $\alpha(t)$, $z(t)$ parametrizadas por tiempo (Ecs. 4.34–4.44 de Góes).

**Suposición clave:** velocidad tangencial constante $V_0$ >> velocidades radial y vertical.

---

### 3.4 SAR/GPR Down-Looking (García-Fernández, Johansson & Mast)

**Marco de referencia:** Cartesiano ENU o $(x,y,z)$ con $z$ positivo hacia abajo (profundidad).

**Sin refracción (aire):**
$$R_{p,m} = \|\mathbf{r}_t - \mathbf{p}\| + \|\mathbf{r}_r - \mathbf{p}\|$$

**Con refracción (Snell, interfaz plana):**
$$R_{eff} = R_{aire} + \sqrt{\varepsilon_r}\,R_{sub}$$

**Punto de refracción** en la interfaz (principio de Fermat):
$$\mathbf{p}^* = \arg\min_{\mathbf{p}\in\text{interfaz}} \left[\frac{|\mathbf{r}_k - \mathbf{p}|}{c} + \frac{|\mathbf{p} - \mathbf{h}|}{v}\right]$$

**Profundidad aparente vs. real:**
$$d_{aparente} = \sqrt{\varepsilon_r}\cdot d_{real}$$

---

### 3.5 SAR Tomografía Subsuperficial (Banda et al. 2015, Imoc 2023)

Combina la geometría SAR multibaseline (apertura en elevación) con propagación en dos medios. El número de onda vertical es el parámetro clave:

$$k_z = \frac{4\pi B_\perp}{\lambda R\sin\theta}$$

**Corrección de posición para hielo** ($\varepsilon_r = 3$):
$$z_{real} = z_{aire} \cdot \frac{1}{\sqrt{\varepsilon_r}} \cdot \frac{\cos\theta_s}{\cos\theta}$$

---

## 4. Elementos Fundamentales Comunes

### 4.1 Distancia Generalizada

Todo el procesado SAR/GPR se reduce a calcular el **camino óptico** desde el transmisor hasta el píxel y vuelta al receptor, en el medio correspondiente:

$$\boxed{R_{OP}(\mathbf{r}_{tx}, \mathbf{p}, \mathbf{r}_{rx}) = \sum_{l} n_l \cdot d_l(\mathbf{r}_{tx}, \mathbf{p}, \mathbf{r}_{rx})}$$

donde $n_l = \sqrt{\varepsilon_{r,l}}$ es el índice de refracción del $l$-ésimo segmento y $d_l$ es la distancia geométrica en ese segmento.

| Caso | Expresión |
|------|-----------|
| Monoestático, aire | $2\|\mathbf{r}_{ant} - \mathbf{p}\|$ |
| Biestático, aire | $\|\mathbf{r}_t - \mathbf{p}\| + \|\mathbf{r}_r - \mathbf{p}\|$ |
| Monoestático, 2 medios | $2(R_{aire} + \sqrt{\varepsilon_r}R_{sub})$ |

---

### 4.2 Fase de la Señal

La fase de la señal SAR demodulada es siempre proporcional al camino óptico:

$$\boxed{\Phi(\mathbf{r}_{ant}, \mathbf{p}) = \frac{4\pi f}{c} \cdot R_{OP} = \frac{4\pi}{\lambda} \cdot R_{OP}}$$

Para el caso biestático con referencia al origen:

$$\Phi = \frac{2\pi}{\lambda}(R_t + R_r) = \frac{4\pi}{\lambda}\cdot\frac{R_t+R_r}{2}$$

La **señal Phase History** (Gorham & Moore, Doerry) es:

$$S(f_k, \tau_n) = A \cdot \exp\!\left(-j\frac{4\pi f_k}{c}\Delta R(\tau_n)\right)$$

donde $\Delta R(\tau_n) = d_{a_0}(\tau_n) - d_a(\tau_n)$ es el rango diferencial respecto al SRP.

---

### 4.3 Relación Fase ↔ Resolución

La resolución en cualquier dirección $\hat{u}$ es el inverso del soporte espectral en esa dirección. Se puede demostrar que:

$$\delta_u = \frac{2\pi}{\Delta k_u}$$

donde $\Delta k_u$ es el rango de números de onda adquiridos en la dirección $\hat{u}$.

**Tabla unificada fase ↔ resolución:**

| Dirección | Número de onda | Rango adquirido | Resolución |
|-----------|---------------|-----------------|------------|
| Rango (monoestático) | $k_r = 4\pi f/c$ | $\Delta k_r = 4\pi B/c$ | $\delta_r = c/(2B)$ |
| Acimut (apertura $\theta_a$) | $k_x = 2k\sin\theta$ | $\Delta k_x = 2k\theta_a$ | $\delta_x = \lambda/(2\theta_a)$ |
| Vertical (tomo. $B_\perp$) | $k_z = 4\pi B_\perp\sin\psi/(\lambda R_0)$ | $\Delta k_z = 4\pi W_z/c$ | $\delta_z \propto c/W_z$ |
| Rango biestático | $k_v = (4\pi f/c)\cos(\varphi/2)$ | $\Delta k_v = (4\pi\Delta f/c)\cos(\varphi/2)$ | $\rho_v = c/(2\Delta f\cos(\varphi/2))$ |
| GPR en profundidad | $k_z = 4\pi f\sqrt{\varepsilon_r}/c$ | $\Delta k_z = 4\pi B\sqrt{\varepsilon_r}/c$ | $\delta_z = c/(2B\sqrt{\varepsilon_r})$ |

**La relación fundamental de incertidumbre SAR:**
$$\delta_r \cdot B = \frac{c}{2}, \quad \delta_{az} \cdot \theta_a = \frac{\lambda}{2}$$

---

### 4.4 Formación de Imagen como Filtro Adaptado (Backprojection)

Todos los papers que describen formación de imagen convergen en la misma operación central:

$$\boxed{I(\mathbf{p}) = \sum_n s\!\left(\frac{R_{OP}(\mathbf{r}_n, \mathbf{p})}{c}\right) \cdot e^{+j\frac{4\pi f_0}{c}R_{OP}(\mathbf{r}_n,\mathbf{p})}}$$

donde $s(\tau)$ es el perfil de rango comprimido (pulse-compressed). Esta es la ecuación del Delay-and-Sum / Backprojection, implementada con diferentes complejidades computacionales:

| Algoritmo | Complejidad | Paper |
|-----------|-------------|-------|
| BP directo (MF) | $\mathcal{O}(N^4)$ para imagen 2D | Gorham & Moore 2010 |
| BP directo (imagen) | $\mathcal{O}(N^3)$ | Ulander 2003, Gorham 2010 |
| Fast BP (Yegulalp) | $\mathcal{O}(N^{3/2}N_{pulse})$ | Yegulalp 1999 |
| FFBP cuadtree (McCorkle) | $\mathcal{O}(N^2\log N)$ | McCorkle 1996 |
| FFBP (Ulander) | $\mathcal{O}(N^2\log N)$ | Ulander 2003 |
| 3D-FFBP (Góes) | $\mathcal{O}(P^3)$ en 3D | Góes 2020 |

---

## 5. Inconsistencias y Diferencias Importantes entre Papers

### 5.1 ★ Error en la Expresión del Desplazamiento Vertical del Número de Onda

**Paper:** Ponce et al. (2016) — citado por Góes 2022  
**Inconsistencia:** Ponce et al. usaba la sensibilidad fase-altura $\partial\varphi/\partial z$ como si fuera $\Delta k_z$:

$$\text{Incorrecto (Ponce 2016):} \quad \Delta k_z^{INCORRECTO} = \frac{4\pi b_\perp}{\lambda R_0\sin\psi}$$

$$\text{Correcto (Góes 2022):} \quad \Delta k_z = \frac{4\pi b_\perp}{\lambda R_0}\sin\psi$$

La diferencia: $\Delta k_z$ y $\partial\varphi/\partial z$ no son lo mismo:

$$\frac{\partial\varphi}{\partial z} = \frac{\Delta k_g}{\tan\psi} + \Delta k_z \neq \Delta k_z$$

**Impacto:** los valores de $\delta_z$ calculados con la expresión incorrecta son **sobrestimados** en Multi-Circular SAR y Spiral SAR.

---

### 5.2 Criterio de Resolución No Uniforme

Distintos papers usan distintas definiciones de "resolución":

| Paper | Criterio | Factor numérico |
|-------|----------|-----------------|
| Munson et al. 1983 | Primer nulo de la PSF ($2\pi/\Delta X$) | depende de geometría |
| Ishimaru (SAR circular) | Punto $e^{-1}$ del valor máximo | factor 1.12 en $\delta_{xy}$ |
| Góes 2022, SAR espiral | Half-power (−3 dB) → factor $\sqrt{\ln(2)/\pi}$ | $\sqrt{\ln(2)/\pi} \approx 0.469$ |
| SAR Tutorial Moreira | Half-power estándar | factor 0.886 para chirp rectangular |
| GPR papers | Generalmente −3 dB o primer nulo | inconsistente entre papers |

**Consecuencia:** las resoluciones numéricas no son directamente comparables entre papers sin especificar el criterio.

---

### 5.3 Notación del Ángulo Biestático

Tres convenciones distintas coexisten:

| Paper | Símbolo | Definición |
|-------|---------|------------|
| Arikan & Munson 1988 | $\beta$ | **Semi-ángulo** biestático; monoestático: $\beta=0$ |
| Horne & Yates | $\varphi$ | **Ángulo completo** entre LOS del TX y RX |
| Ding & Munson 2002 | $\beta$ | Ángulo biestático; factor de escala $2\cos(\beta/2)/c$ |

El factor de escala k-espacio se expresa como:
- Arikan: $w_{tr}/c = 2\cos(\beta)/c$
- Ding: $(2/c)\cos(\beta/2)$
- Horne: $\rho_v = c/(2\Delta f\cos(\varphi/2))$

Todas son equivalentes: $\cos\beta_{Arikan} = \cos(\varphi_{Horne}/2) = \cos(\beta_{Ding}/2)$, es decir, los tres papeles definen el mismo ángulo de forma diferente.

---

### 5.4 Tratamiento de la Refracción GPR

Hay un gradiente de sofisticación:

| Paper | Modelo de refracción | Aproximación |
|-------|---------------------|--------------|
| Johansson & Mast 1994 | Bicapa, Snell iterativo | $x_b = x_2 + \sqrt{\varepsilon_{r1}/\varepsilon_{r2}}(x_1-x_2)$ |
| Imoc/BP refractive 2023 | Interfaz plana o inclinada, Fermat | Exacto en 2D, numérico en 3D |
| Banda et al. 2015 | Interfaz plana, corrección de profundidad | $z = z_{air}\cos\theta_s/(\sqrt{\varepsilon_r}\cos\theta)$ |
| García-Fernández 2019 | Tres capas, fórmula analítica approx. | Ecs. 3–5 de González-Díaz 2019 |
| García-Fernández 2019 (UAV) | Snell con permitividad estimada | $d_{aparente} = \sqrt{\varepsilon_r}\,d_{real}$ |
| Sistema Patente 2016 | Velocidad uniforme $v = c/\sqrt{\varepsilon_r}$ | Estimación de $\varepsilon_r$ in situ |

Los modelos más simples (sin refracción) **no ubican correctamente** los objetos subsuperficiales (error $\propto\sqrt{\varepsilon_r}$ en profundidad).

---

### 5.5 Velocidad de Propagación: Constante vs. Variable

| Supuesto | Papers | Riesgo |
|---------|--------|--------|
| $v = c$ uniforme | Todos los SAR aéreos/satelitales | OK para aire; erróneo para subsuelo |
| $v = c/\sqrt{\varepsilon_r}$, $\varepsilon_r$ constante | GPR papers (la mayoría) | Error si $\varepsilon_r$ varía con profundidad |
| $v$ variable en capas | González-Díaz 2019, Johansson 1994 | Modelo más realista |
| $v$ variable (FDTD) | Irving & Knight 2006 | Modelo completo (solo simulación) |

---

### 5.6 Modelo de Señal SAR: Chirp LFM vs. Impulso vs. M-sequence

| Paper | Tipo de señal | Implicación |
|-------|--------------|-------------|
| Munson 1983, Arikan 1988, Horne/Yates | Chirp LFM continuo | Resolución = $c/(2B)$ directa |
| Doerry (Basics BP) | Chirp LFM + stretch processing (RVPE) | Corrección de fase cuadrática necesaria |
| Gorham & Moore 2010 | Stepped-frequency (datos AFRL) | Requiere IFFT para formar perfil de rango |
| García-Fernández 2019 (GPR) | M-sequence UWB (correlación cruzada) | Equivalente a impulso; BW muy grande |
| UAV multichannel 2024 | Señal UWB 1–6 GHz | Procesado por sub-bandas recomendado |

---

### 5.7 Suposición Stop-and-Go: Universal en SAR, Ignorada en GPR

Todos los papers SAR la mencionan explícitamente como suposición. Los papers GPR (García-Fernández, Johansson) la asumen implícitamente sin mencionarla. La suposición es más crítica para GPR-UAV a baja altitud y alta velocidad relativa del escáner.

---

## 6. Modelo Conceptual Unificado

### 6.1 Jerarquía de Modelos SAR/GPR

```
MODELO GENERAL
│
├── Definición de señal:  s(f) = A·exp(-j·4πf/c·R_OP)
│
├── Distancia óptica:  R_OP = Σ nₗ·dₗ  (índices de refracción por segmento)
│   ├── Monoestático en aire:  R_OP = 2‖r_ant - p‖
│   ├── Biestático en aire:   R_OP = ‖r_t - p‖ + ‖r_r - p‖
│   └── Cualquier medio:  R_OP = R_aire + √εr·R_sub  (+ Snell en interfaz)
│
├── Fase:  Φ = (4πf/c)·R_OP
│
├── Formación de imagen (Backprojection):
│   I(p) = Σₙ s(R_OP(rₙ,p)/c)·exp(+j·4πf₀/c·R_OP(rₙ,p))
│
└── Resolución:  δᵤ = 2π/Δkᵤ  en cada dirección ŷ
    ├── Rango:    Δkᵣ = 4πB/c          →  δᵣ = c/(2B)
    ├── Azimut:   Δkₓ = 2k·θₐ          →  δₓ = λ/(2θₐ)
    ├── Vertical: Δkz = 4πB⊥sinψ/(λR₀) →  δz = f(Wz, λ, B⊥, ψ₀, R₀)
    └── Biestático: factor w_tr = 2cosβ  →  δᵥ = c/(2Δf·cosβ)
```

---

### 6.2 Parámetros Clave del Sistema Unificado

| Parámetro | Símbolo | Controla |
|-----------|---------|---------|
| Velocidad en el medio | $v = c/\sqrt{\varepsilon_r}$ | Resolución en rango, propagación |
| Ancho de banda | $B$ | Resolución en rango: $\delta_r = v/(2B)$ |
| Longitud de onda | $\lambda = v/f_c$ | Resolución en azimut y vertical |
| Ángulo de apertura | $\theta_a$ | Resolución en acimut: $\delta_x = \lambda/(2\theta_a)$ |
| Apertura tomográfica efectiva | $B_\perp = B|\cos(\beta-\psi_0)|$ | Resolución vertical: $\delta_z \propto 1/W_z(B_\perp)$ |
| Ángulo de incidencia | $\psi_0$ | Mezcla rango/vertical; afecta $B_\perp$ y $W_z$ |
| Ángulo biestático | $\beta$ (semi-ángulo) | Reduce resolución: factor $\cos\beta$ en k-espacio |
| Permitividad del subsuelo | $\varepsilon_r$ | Modifica velocidad, introduce refracción, $d_{aparente}/d_{real} = \sqrt{\varepsilon_r}$ |

---

### 6.3 Tres Regímenes de Operación

**Régimen 1 — SAR aéreo/satelital (propagación en aire):**
- $\varepsilon_r = 1$, $v = c$
- Foco en resolución 2D (rango + azimut) o 3D (+ tomografía/espiral)
- Papers: SAR Tutorial, Ulander, Doerry, Munson, Arikan, Horne, Góes, Moreira DInSAR, Banda

**Régimen 2 — GPR downward-looking (propagación en dos medios):**
- $\varepsilon_r$ del suelo necesario para localizar en profundidad
- La refracción en la interfaz es la corrección más importante (sin ella: error $\propto\sqrt{\varepsilon_r}$ en posición)
- Papers: García-Fernández 2019 (x2), UAV multichannel 2024, sistema patente 2016, Johansson & Mast 1994

**Régimen 3 — Modelado numérico EM:**
- FDTD 2D (Irving & Knight 2006): simula la propagación completa de Maxwell
- Útil para validación de algoritmos de formación de imagen
- No deriva directamente resoluciones: las calcula numéricamente

---

### 6.4 Tabla de Geometrías y sus Ecuaciones de Resolución

| Geometría | $\delta_r$ (rango) | $\delta_{az}$ (azimut) | $\delta_z$ (vertical) | Papers |
|-----------|-------------------|-----------------------|----------------------|--------|
| SAR lineal monoestático | $c/(2B)$ | $d_a/2$ | — | SAR Tutorial, Doerry, Gorham |
| SAR spotlight (tomo.) | $c/(2B)$ | $\lambda/(2\theta_M)$ | — | Munson 1983 |
| SAR biestático spotlight | $c/(2\Delta f\cos\beta)$ | $c/(f_0 K)$ | — | Horne, Arikan, Ding |
| Circular SAR (1 vuelta) | $c/(2B)$ | $1.12\lambda/(2\pi\sin\psi)$ | $\sqrt{\ln2/\pi\cdot c/(W\cos\psi)}$ | Góes 2022 (cita Ishimaru) |
| Spiral/Multi-Circular SAR | $c/(2B)$ | $1.12\lambda/(2\pi\sin\psi_0)$ | $\sqrt{\ln2/\pi\cdot c/W_z}$ | Góes 2022 |
| SAR tomografía multibaseline | $c/(2B)$ | $d_a/2$ (por pasada) | $\lambda R_0/(2B_\perp)$ | SAR Tutorial, Banda 2015 |
| GPR downward-looking (aire) | $c/(2B)$ | $\lambda_c\sqrt{L^2/4+h^2}/(2L)$ | — | García-F. 2024 |
| GPR downward-looking (suelo) | $c/(2B\sqrt{\varepsilon_r})$ | idem | $c/(2B\sqrt{\varepsilon_r})$ | todos los GPR papers |

---

### 6.5 La Cadena de Procesado Unificada

```
DATOS CRUDOS (phase history)
    │
    ▼
1. CORRECCIÓN DE TRAYECTORIA (MoComp)
   ΔR(τₙ) = d_a₀(τₙ) - dₐ(τₙ)   [referencia al SRP]
    │
    ▼
2. COMPRESIÓN EN RANGO
   s(R,τₙ) = IFFT{S(f,τₙ)}   [perfil de rango]
    │
    ▼
3. CORRECCIÓN DE MEDIO [GPR solamente]
   R_OP = R_aire + √εr·R_sub   [corrección de refracción]
    │
    ▼
4. BACKPROJECTION / DELAY-AND-SUM
   I(p) = Σₙ s(R_OP(rₙ,p)/c)·exp(+j·4πf₀/c·R_OP)
    │
    ▼
5. IMAGEN SAR/GPR 3D
   → Resolución determinada por {B, θₐ, B_⊥, λ, ψ₀, εr}
```

---

### 6.6 Observación Final: La Unidad del Modelo

A pesar de la diversidad de geometrías, configuraciones y aplicaciones cubiertos por los 22 papers, **todo el procesado SAR/GPR puede describirse con tres ecuaciones fundamentales**:

1. **Fase (camino óptico):**
$$\Phi(\mathbf{r}_{ant}, \mathbf{p}) = \frac{4\pi f}{c}\,R_{OP}(\mathbf{r}_{ant}, \mathbf{p})$$

2. **Imagen (filtro adaptado/backprojection):**
$$I(\mathbf{p}) = \sum_n \int S(f, \tau_n)\,e^{+j\frac{4\pi f}{c}R_{OP}(\mathbf{r}_n,\mathbf{p})}\,df$$

3. **Resolución (dualidad k-espacio):**
$$\delta_u = \frac{2\pi}{\Delta k_u} = \frac{\pi c}{2f\Delta\theta_u}$$

La riqueza del campo no reside en ecuaciones fundamentales distintas, sino en cómo varía $R_{OP}$ (la distancia óptica) según la geometría: monoestática vs. biestática, lineal vs. circular vs. espiral, en aire vs. en subsuelo con refracción. **Cambiar la geometría es cambiar la definición de $R_{OP}$**; el resto del modelo permanece invariante.

---

## Referencias Cruzadas

| Concepto | Papers relevantes |
|----------|-------------------|
| Resolución $\delta_r = c/(2B)$ | Todos (22/22) |
| Backprojection como suma coherente | Ulander, Yegulalp, McCorkle, Doerry, Gorham, Góes, García-F., Johansson |
| Formulación tomográfica SAR | Munson 1983, Arikan 1988, SAR Tutorial |
| Factor biestático $w_{tr} = 2\cos\beta$ | Arikan 1988, Ding 2002 |
| Refracción en interfaz aire/suelo | Imoc 2023, Banda 2015, García-F. 2019×2, Johansson 1994, Patente 2016 |
| Resolución vertical espiral | Góes 2022 (expresión corregida respecto a Ponce 2016) |
| Desplazamiento de número de onda | Góes 2022, SAR Tutorial (InSAR), Drone DInSAR |
| Algoritmos FFBP | Ulander 2003, Yegulalp 1999, McCorkle 1996, Góes 2020 |
| Sistemas reales P-band dron | Góes 2022, Banda 2015, Drone DInSAR |
| GPR-SAR UAV experimental | García-F. 2019, UAV multichannel 2024, Patente 2016 |
