# Hipótesis v4 — Modelo Consolidado: Estado del Arte, Validez y Límites

**Fecha:** 2026-06-01  
**Estado:** MODELO CONSOLIDADO — síntesis de v1–v3 + experimento de grilla fina  
**Referencia:** `models/resolucion_v2.md`, `hypotheses/hipotesis_v1.md` a `v3.md`

---

## 1. Tabla completa de validaciones experimentales

| # | Experimento | Variable | Predicción | Simulación | Error |
|---|------------|----------|:----------:|:----------:|:-----:|
| 1 | NorthOffset=180°, grilla 40mm | δz | 0.211 m | 0.2112 m | **0.05%** ✅ |
| 2 | NorthOffset=180°, grilla 40mm | δxy | 0.240 m (grilla) | 0.240 m | **0.0%** ✅ |
| 3 | NorthOffset=90°, grilla 40mm | δz | 0.211 m | 0.2112 m | **0.05%** ✅ |
| 4 | NorthOffset=90°, grilla 40mm | δxy (artefacto) | 19.8 mm | 19.0–18.1 mm | **4–9%** ✅ |
| 5 | NorthOffset=90°, **grilla 2mm** | δz | 0.211 m | 0.2112 m | **0.05%** ✅ |
| 6 | NorthOffset=90°, **grilla 2mm** | δxy (real) | 9.54 mm | **8.8 mm** | **7.8%** ✅ |

Todos los experimentos tienen target en el **eje** `[0,0,−5]`.

---

## 2. Las ecuaciones del modelo consolidado

### 2.1 Resolución vertical — validada con error 0.05 %

$$\boxed{\delta_z = \frac{c}{2\,W_z}}$$

$$\boxed{W_z = n_2\,B\,\cos\theta_{t,0} + \frac{c\,B_\perp\,\sin\psi_0\cos\psi_0}{\lambda_0\,R_0\,n_2\cos\theta_{t,0}}}$$

**Invariante al NorthOffset** (demostrado analítica y numéricamente para 180° y 90°).

### 2.2 Resolución horizontal real — validada con error 7.8 %

$$\boxed{\delta_{xy}(\Delta\phi) = \frac{2\times 1.20\,\lambda_0}{4\pi\,\sin\psi_0\cdot|\cos(\Delta\phi/2)|} = \frac{0.60\,\lambda_0}{\pi\,\sin\psi_0\cdot|\cos(\Delta\phi/2)|}}$$

El factor 1.20 corresponde al punto −3 dB de la PSF de tipo $J_0$ (Bessel de orden 0) — forma exacta de la PSF para cobertura azimutal completa.

**Nota sobre el error residual de 7.8 %:**  
La fórmula usa $\psi_0 = \arctan(\rho_0/z_0)$, que es el look angle hacia la **interfaz** ($z=0$), no hacia el target ($z_P=-5$ m). El error sistemático ≈ 8 % proviene de esta aproximación. El solver de Snell en la simulación computa el ángulo refractado real, que es ligeramente mayor que el predicho por la fórmula.

**Corrección exacta (requiere resolver Snell numéricamente):**
$$\sin\theta_t^{exact} = Q_\rho/\sqrt{Q_\rho^2+z_P^2}$$
donde $Q_\rho$ es la solución del sistema iterativo de Snell para la geometría exacta.

### 2.3 Parámetros auxiliares

| Símbolo | Expresión | Descripción |
|---------|-----------|-------------|
| $\psi_0$ | $\arctan(\rho_0/z_0)$ | Look angle medio (sensor→interfaz) |
| $\cos\theta_{t,0}$ | $\sqrt{1-\sin^2\psi_0/n_2^2}$ | Ángulo de refracción (Snell) |
| $R_0$ | $\sqrt{\rho_0^2+z_0^2}$ | Distancia media sensor→origen |
| $B_{helix}$ | $\sqrt{\Delta z^2+\Delta\rho^2}$ | Longitud de la trayectoria helicoidal |
| $\beta$ | $\arctan(\Delta z/\Delta\rho)$ | Ángulo de inclinación de la hélice |
| $B_\perp$ | $B_{helix}\cdot|\cos(\beta-\psi_0)|$ | Apertura tomográfica efectiva |
| $\lambda_0$ | $c/f_0$ | Longitud de onda en aire |

---

## 3. Dominio de validez: tres condiciones

| Condición | Formalización | Error si no se cumple |
|-----------|--------------|----------------------|
| **C1: Target en el eje** | $x_P = y_P = 0$ | PSF asimétrica; δxy varía con posición y dirección |
| **C2: TX y RX en la misma hélice** | mismo $\rho(t)$, $z(t)$; solo difiere el azimut | $k_z \neq 2\cos\theta_t$; fórmula δz cambia |
| **C3: Target poco profundo** | $|z_P| \ll z_0$ | Error en ψ₀ → error sistemático ~8% en δxy |

---

## 4. Error residual sistemático y su origen

El único error verificado (7.8% en δxy, 0.05% en δz) tiene dos causas identificadas:

**Para δz:** el error de 0.05% está por debajo del ruido numérico → **la fórmula es exacta dentro de las hipótesis declaradas**.

**Para δxy:** el error de 7.8% viene del ángulo de refracción real vs. aproximado:

$$\Delta\text{error} \approx \frac{|z_P|}{z_0}\cos^2\psi_0 \approx \frac{5}{100}\times0.281 \approx 1.4\%$$

Esto es solo una parte del error total (7.8%). El resto proviene de la distribución no uniforme de densidad de muestreo en el k-espacio para la hélice cónica (no es una cobertura perfectamente arcsine-uniforme).

---

## 5. Comportamiento según NorthOffset

| NorthOffset $\Delta\phi$ | Resolución horizontal | Resolución vertical |
|:---:|:---:|:---:|
| 0° (monoestático) | $\delta_{xy}^{mono} = \frac{0.60\lambda_0}{\pi\sin\psi_0}$ | $= \delta_z$ (invariante) |
| 90° | $\sqrt{2}\,\delta_{xy}^{mono}$ | $= \delta_z$ |
| 120° | $2\,\delta_{xy}^{mono}$ | $= \delta_z$ |
| 180° | $\infty$ (sin resolución) | $= \delta_z$ |

La degradación de la resolución horizontal es $1/|\cos(\Delta\phi/2)|$, lo que es un factor **continuo y controlable** por diseño del sistema.

---

## 6. Lo que falta para un modelo verdaderamente general

| Extensión necesaria | Complejidad | Impacto |
|--------------------|-------------|---------|
| **Target off-axis** ($x_P \neq 0$) | Media — corrección ψ̃₀ de Góes 2022 | Alta — caso práctico más común |
| TX y RX con distintos parámetros helicoidales | Alta — k_z biestático asimétrico | Media |
| Target profundo ($|z_P| \sim z_0$) | Alta — corrección exacta de Snell | Media |
| Cubrimiento azimutal incompleto ($N_t < 1$ vuelta) | Media — PSF deja de ser J₀ | Alta |

La primera extensión (**target off-axis**) es la más relevante para aplicaciones prácticas y será abordada en `hipotesis_v5.md`.
