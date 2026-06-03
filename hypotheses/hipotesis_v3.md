# Hipótesis v3 — Modelo Extendido: Resolución para NorthOffset Arbitrario

**Fecha:** 2026-06-01  
**Estado:** VALIDADO por simulación (NorthOffset=90°)  
**Motivación:** hipotesis_v2 predijo incorrectamente δx=δy=0.24m para NorthOffset=90°.  
**Corrección:** el modelo horizontal necesita incluir el perfil J₀ de la PSF y su interacción con la grilla.

---

## 1. Comparación de predicciones v2 vs. simulación (NorthOffset=90°)

| Eje | v2 predijo | Simulación | Error |
|-----|:----------:|:----------:|:-----:|
| δz | 0.211 m | 0.2112 m | **0.05 %** ✅ |
| δx | 0.24 m (grilla) | **0.0190 m** | ❌ Predicción incorrecta |
| δy | 0.24 m (grilla) | **0.0181 m** | ❌ Predicción incorrecta |

**¿Por qué falló la predicción horizontal?** v2 asumió que k_x=k_y=0 para NorthOffset=180° se extendería a 90°. Error: para NorthOffset=90°, los vectores k horizontales NO se cancelan completamente.

---

## 2. Corrección: análisis del k-espacio horizontal para NorthOffset genérico

### 2.1 Fórmula general de k_x para cualquier offset azimutal Δφ

Con TX en ángulo $\alpha$ y RX en $\alpha + \Delta\phi$ (para target en el eje):

$$\hat{\mathbf{e}}_2^{TX}\cdot\hat{x} = -\sin\theta_t\cos\alpha$$
$$\hat{\mathbf{e}}_2^{RX}\cdot\hat{x} = -\sin\theta_t\cos(\alpha+\Delta\phi)$$

$$k_x(\alpha) = \frac{2\pi f n_2}{c}\cdot\sin\theta_t\cdot\big(-\cos\alpha - \cos(\alpha+\Delta\phi)\big)$$

Usando identidad trigonométrica: $\cos\alpha + \cos(\alpha+\Delta\phi) = 2\cos(\Delta\phi/2)\cos(\alpha+\Delta\phi/2)$

$$\boxed{k_x(\alpha) = -\frac{2\pi f n_2}{c}\cdot 2\sin\theta_t\cos(\Delta\phi/2)\cdot\cos(\alpha+\Delta\phi/2)}$$

La amplitud máxima varía con $\Delta\phi$:

$$k_x^{max}(\Delta\phi) = \frac{4\pi f_0 n_2}{c}\sin\theta_{t,0}\cdot|\cos(\Delta\phi/2)|$$

| NorthOffset $\Delta\phi$ | $|\cos(\Delta\phi/2)|$ | $k_x^{max}$ relativo |
|:---:|:---:|:---:|
| 0° (monoestático) | 1 | $R_c^{mono}$ |
| 90° | $1/\sqrt{2}$ | $R_c^{mono}/\sqrt{2}$ |
| 180° | 0 | 0 (sin resolución horizontal) |

### 2.2 Forma del k-espacio horizontal

Para NorthOffset $\Delta\phi$, al variar $\alpha \in [0, 2\pi N_t]$, la proyección $(k_x, k_y)$ traza un **círculo** de radio:

$$R_c(\Delta\phi) = \frac{4\pi f_0 n_2}{c}\cdot\sin\theta_{t,0}\cdot|\cos(\Delta\phi/2)|$$

Verificación:
- Δφ=0°: $R_c = (4\pi f_0 n_2/c)\sin\theta_{t,0}$ = monoestático ✓  
- Δφ=90°: $R_c = (4\pi f_0 n_2/c)\sin\theta_{t,0}/\sqrt{2} = \sqrt{2} \times (2\pi f_0 n_2/c)\sin\theta_{t,0}$ ✓  
- Δφ=180°: $R_c = 0$ → sin cobertura horizontal ✓

La PSF horizontal es:
$$\text{PSF}_{xy}(\delta r) = J_0(R_c(\Delta\phi)\cdot \delta r)$$

### 2.3 Resolución horizontal verdadera (−3 dB)

La resolución horizontal **real** (media-potencia, criterio −3 dB):

$$\boxed{\delta_{xy}^{true}(\Delta\phi) = \frac{1.20}{R_c(\Delta\phi)} = \frac{1.20\,\lambda_0}{4\pi\sin\psi_0\cdot|\cos(\Delta\phi/2)|}}$$

Para $\Delta\phi = 90°$: $\delta_{xy}^{true} = 1.20\lambda_0/(4\pi\sin\psi_0\cdot\frac{1}{\sqrt{2}}) = \sqrt{2}\cdot\delta_{xy}^{mono}$

---

## 3. Fenómeno observado en la simulación: artefacto de grilla con J₀

La simulación mide δx ≈ 19 mm, pero la resolución verdadera es ≈ 9.6 mm. Esto no es un error — es un **artefacto de submuestreo** de la PSF tipo J₀ en una grilla gruesa.

### 3.1 Mecanismo

La grilla de procesado tiene `dxy = 0.04 m = 40 mm`. La PSF real es $J_0(R_c\cdot\delta x)$ con FWHM ≈ 9.6 mm, mucho más estrecha que la grilla. Los puntos adyacentes al pico están en las zonas de **lóbulos laterales** de J₀.

Perfil x simulado (confirmado por análisis post-simulación):

| x [m] | −0.12 | −0.08 | **−0.04** | **0.00** | **+0.04** | +0.08 | +0.12 |
|-------|-------|-------|-----------|---------|-----------|-------|-------|
| dB | −32.1 | −19.1 | **−12.58** | **0.0** | **−12.62** | −19.2 | −31.8 |

El punto ±40 mm está en el **lóbulo lateral** de J₀, a −12.6 dB (no en la zona de −3 dB del lóbulo principal).

### 3.2 Predicción del artefacto

La función `report3dBResolution` interpolará linealmente entre el pico (0 dB en x=0) y el punto adyacente (−12.6 dB en x=40 mm) para encontrar el cruce en −3 dB:

$$x_{app} = d_{xy}\cdot\frac{-3\,\text{dB}}{20\log_{10}(|J_0(R_c\cdot d_{xy})|)} = d_{xy}\cdot\frac{3\,\text{dB}}{|20\log_{10}(|J_0(R_c\cdot d_{xy})|)|}$$

$$\delta_{xy}^{app} = 2\,x_{app} = \frac{2\,d_{xy}\cdot 3}{|20\log_{10}(|J_0(R_c\cdot d_{xy})|)|}$$

**Cálculo numérico para $\Delta\phi = 90°$:**

| Magnitud | Valor |
|---------|-------|
| $R_c$ | 251.2 m⁻¹ |
| $R_c \cdot d_{xy}$ | $251.2 \times 0.04 = 10.05$ |
| $J_0(10.05)$ | $-0.2477$ |
| $20\log_{10}(0.2477)$ | $-12.12$ dB |
| $x_{app} = 0.04 \times 3/12.12$ | 9.90 mm |
| $\delta_{xy}^{app} = 2 \times 9.90$ | **19.8 mm** |

**Simulación:** δx = 19.0 mm, δy = 18.1 mm → **error 4–9 %** ✅

La pequeña diferencia (19.8 mm predicho vs 19.0-18.1 mm medido) se debe a:
1. La PSF no es exactamente $J_0$ (densidad de muestreo no perfectamente uniforme en $k_x$ por la hélice cónica)
2. Leve asimetría bistática para target a z=−5 m (no en la interfaz)
3. La interpolación lineal no es exacta para una función que curva como J₀

---

## 4. Validación del modelo de δz — Invariancia respecto a NorthOffset

**Resultado clave:** la resolución vertical es independiente del NorthOffset azimutal.

**Demostración:**

Para cualquier NorthOffset $\Delta\phi$, con TX y RX en la misma hélice (mismo $\rho(t)$, $z(t)$):

$$k_z(f,t) = \frac{2\pi f n_2}{c}\!\left(\hat{\mathbf{e}}_2^{TX}\cdot\hat{z} + \hat{\mathbf{e}}_2^{RX}\cdot\hat{z}\right) = \frac{2\pi f n_2}{c}\!\left(-\cos\theta_t - \cos\theta_t\right) = -\frac{4\pi f n_2}{c}\cos\theta_t$$

La componente vertical de ambos vectores unitarios es $-\cos\theta_t$ **independientemente del azimut**, porque TX y RX están a la misma altura $z(t)$ y el target está en el eje. ∴ $k_z$ no depende de $\Delta\phi$.

**Confirmación numérica:**

| NorthOffset | δz predicho | δz simulado | Error |
|:-----------:|:-----------:|:-----------:|:-----:|
| 180° | 0.211 m | 0.2112 m | 0.05% |
| **90°** | **0.211 m** | **0.2112 m** | **0.05%** |

Idénticos hasta los 4 decimales. ✅

---

## 5. Modelo completo validado v3

### 5.1 Resolución vertical (invariante a NorthOffset)

$$\boxed{\delta_z = \frac{c}{2\,W_z}}$$

$$W_z = n_2\,B\,\cos\theta_{t,0} + \frac{c\,B_\perp\,\sin\psi_0\cos\psi_0}{\lambda_0\,R_0\,n_2\cos\theta_{t,0}}$$

**Condición:** target en el eje ($x_P=y_P=0$), TX y RX en la misma hélice con igual $\rho(t)$ y $z(t)$.

### 5.2 Resolución horizontal real

$$\boxed{\delta_{xy}^{true}(\Delta\phi) = \frac{1.20\,\lambda_0}{4\pi\,\sin\psi_0\cdot|\cos(\Delta\phi/2)|}}$$

| $\Delta\phi$ | Factored $|\cos(\Delta\phi/2)|$ | $\delta_{xy}$ relativo |
|:---:|:---:|:---:|
| 0° | 1 | $\delta^{mono}$ |
| 60° | $\cos(30°) = \sqrt{3}/2 \approx 0.866$ | $1.15\,\delta^{mono}$ |
| 90° | $1/\sqrt{2} \approx 0.707$ | $\sqrt{2}\,\delta^{mono}$ |
| 120° | $\cos(60°) = 0.5$ | $2\,\delta^{mono}$ |
| 180° | 0 | $\infty$ |

### 5.3 Resolución horizontal medida por simulación (artefacto de grilla)

Cuando la PSF verdadera es mucho más estrecha que la grilla (`dxy >> δ_true`):

$$\boxed{\delta_{xy}^{app} = \frac{2\,d_{xy}\cdot 3}{|20\log_{10}(|J_0(R_c(\Delta\phi)\cdot d_{xy})|)|}}$$

Esta expresión predice lo que reportará `report3dBResolution` en la condición de submuestreo.

**Para NorthOffset=90°:** predicción = 19.8 mm, simulación = 19.0-18.1 mm (error 4-9%).

---

## 6. Tabla de verificación completa

| Experimento | Parámetro | Predicción v3 | Simulación | Error |
|-------------|-----------|:-------------:|:----------:|:-----:|
| NorthOffset=180° | δz | 0.211 m | 0.2112 m | 0.05% ✅ |
| NorthOffset=180° | δx=δy | 0.24 m (grilla) | 0.2400 m | 0.0% ✅ |
| NorthOffset=90° | δz | 0.211 m | 0.2112 m | 0.05% ✅ |
| NorthOffset=90° | δx_app | 19.8 mm | 19.0 mm | 4.2% ✅ |
| NorthOffset=90° | δy_app | 19.8 mm | 18.1 mm | 9.4% ✅ |

---

## 7. Conclusión física

La configuración biestática helicoidal con TX y RX en la misma trayectoria a offset azimutal $\Delta\phi$ tiene las siguientes propiedades para target en el eje:

1. **La resolución vertical $\delta_z$ es invariante a $\Delta\phi$** — depende únicamente del ancho de banda, la apertura tomográfica y el ángulo de incidencia.

2. **La resolución horizontal verdadera $\delta_{xy}^{true}$ degrada con $\Delta\phi$** — empeora por el factor $1/|\cos(\Delta\phi/2)|$ respecto al monoestático, hasta desaparecer en $\Delta\phi=180°$.

3. **La resolución horizontal medida por simulación (grilla gruesa) no refleja la resolución real** para el caso de 10 GHz — es un artefacto de la interacción entre el lóbulo lateral de J₀ y el espaciado de grilla de 40 mm.

4. **Para medir la resolución horizontal verdadera** se necesita `dxy ≤ λ₀/(4sinψ₀·|cos(Δφ/2)|) / 2` ≈ 2.4 mm para el sistema actual. O bien, cambiar a P-band donde la resolución es de cm en lugar de mm.
