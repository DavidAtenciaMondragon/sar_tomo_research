# Hipótesis v2 — Modelo de Resolución Validado y Generalizado

**Fecha:** 2026-06-01  
**Estado:** VALIDADO por simulación  
**Configuración:** `_espiral` · NorthOffset TX=0°, RX=180° · n₂=2 · Target `[0,0,−5]`

---

## 1. Resultado de la verificación de v1

| Eje | Predicción v1 | Simulación (-3 dB) | Error relativo |
|-----|:-------------:|:------------------:|:--------------:|
| **δx** | 0.240 m (grilla) | **0.2400 m** | 0.0 % |
| **δy** | 0.240 m (grilla) | **0.2400 m** | 0.0 % |
| **δz** | 0.211 m | **0.2112 m** | 0.05 % |

Pico detectado exactamente en `[0, 0, −5]` = posición real del target ✓

La predicción para δz es **prácticamente exacta** (error 0.05%). La predicción para δx y δy también es exacta: el sistema no tiene resolución horizontal para target en el eje con NorthOffset=180°, y los 0.24 m reportados son el límite de la grilla de procesado.

---

## 2. Crítica de la hipótesis v1

### 2.1 Lo que v1 acertó

1. **Cancelación exacta de k_x y k_y** — El razonamiento de simetría era correcto: para target en el eje y TX/RX diametralmente opuestos, los vectores horizontales del k-espacio se cancelan exactamente, dejando cobertura solo en k_z.

2. **Fórmula de δz** — La expresión $\delta_z = c/(2W_z)$ predicó el valor -3 dB medido con error de 0.05%. Esto valida el modelo de Góes 2022 también para X-band (10 GHz).

3. **Límite de grilla para x, y** — La predicción de que `report3dBResolution` reportaría exactamente el span de la grilla (0.24 m) fue correcta.

### 2.2 Ambigüedad que debe aclararse: criterio de resolución

En v1 se derivó $\delta_z = c/(2W_z)$ presentándolo como "criterio de primer nulo". Esto es incorrecto: la fórmula en realidad produce directamente el valor **−3 dB** para la PSF del sistema helicoidal con mezcla de contribuciones frecuencial y geométrica.

Comparación de criterios para δz de este sistema:

| Criterio | Fórmula | Valor | Coincidencia con sim. |
|---------|---------|-------|----------------------|
| Primer nulo (sinc) | $2\pi/\Delta k_z = 2\pi/29.0$ | 0.217 m | Regular (2.8% error) |
| −3 dB (Gaussiana, HWHM) | $\sqrt{\ln 2/\pi} \cdot c/W_z$ | 0.199 m | Peor (5.8% error) |
| **Fórmula Góes** | $c/(2W_z)$ | **0.211 m** | **Exacto (0.05%)** |

La fórmula de Góes $\delta_z = c/(2W_z)$ es un **resultado empírico-analítico** calibrado para la PSF real del SAR helicoidal. No es el primer nulo ni la HWHM gaussiana pura. Predice el −3 dB directamente porque incorpora la densidad de muestreo real del k-espacio (no uniforme en $k_z$).

### 2.3 Condición de dominancia del término tomográfico

La contribución tomográfica dominó sobre la de BW en un factor 6.8×. Esto solo ocurre porque:
- $f_0 = 10$ GHz es alta → λ pequeño → alta sensibilidad angular
- $B_\perp = 47$ m es grande (hélice cónica óptima con $\beta = \psi_0$)

Para P-band ($f_0 = 425$ MHz), la contribución de BW es mayor y el término tomográfico es comparable o menor. El modelo se aplica igual pero los valores cambian.

---

## 3. Modelo validado — Fórmula de resolución

### 3.1 Condición de aplicabilidad

La fórmula aplica cuando:

$$\text{NorthOffset}_{RX} - \text{NorthOffset}_{TX} = 180° \quad \land \quad x_P = y_P = 0$$

(target en el eje del helicoide, TX y RX diametralmente opuestos).

### 3.2 Resolución horizontal — Sin resolución (cancela exactamente)

Para target en el eje:

$$k_x(f,t) = k_y(f,t) = 0 \quad \forall\; (f,t)$$

**Consecuencia:** No existe apertura horizontal → PSF plana en (x,y) → medición queda acotada por la grilla de procesado:

$$\delta_x^{medido} = \delta_y^{medido} = 2\,l_{xy}$$

### 3.3 Resolución vertical — Fórmula validada

$$\boxed{\delta_z = \frac{c}{2\,W_z}}$$

$$\boxed{W_z = n_2\,B\,\cos\theta_{t,0} + \frac{c\,B_\perp\,\sin\psi_0\cos\psi_0}{\lambda_0\,R_0\,n_2\cos\theta_{t,0}}}$$

**Definición de cada término:**

| Símbolo | Expresión | Descripción |
|---------|-----------|-------------|
| $n_2$ | $\sqrt{\varepsilon_r}$ | Índice de refracción del suelo |
| $B$ | $f_{max}-f_{min}$ | Ancho de banda del chirp [Hz] |
| $\psi_0$ | $\arctan(\rho_0/z_0)$ | Look angle medio (sensor→interfaz) |
| $\theta_{t,0}$ | $\arcsin(\sin\psi_0/n_2)$ | Ángulo de transmisión en suelo (Snell) |
| $\cos\theta_{t,0}$ | $\sqrt{1-\sin^2\psi_0/n_2^2}$ | Componente vertical del rayo en suelo |
| $\lambda_0$ | $c/f_0$ | Longitud de onda en aire |
| $R_0$ | $\sqrt{\rho_0^2+z_0^2}$ | Distancia media sensor→target proyectada |
| $B_{helix}$ | $\sqrt{\Delta z^2+\Delta\rho^2}$ | Longitud total de la trayectoria helicoidal |
| $\beta$ | $\arctan(\Delta z/\Delta\rho)$ | Ángulo de inclinación de la hélice |
| $B_\perp$ | $B_{helix}|\cos(\beta-\psi_0)|$ | Apertura tomográfica efectiva |

**Condición óptima:** $\beta = \psi_0$ → $B_\perp = B_{helix}$ (máxima resolución vertical).

### 3.4 Verificación numérica para la config `_espiral`

| Parámetro | Valor |
|-----------|-------|
| $f_0$ | 10 GHz |
| $B$ | 50 MHz |
| $n_2$ | 2 |
| $\rho_0$ | 160 m |
| $z_0$ | 100 m |
| $\psi_0$ | 58.0° |
| $\cos\theta_{t,0}$ | 0.906 |
| $B_{helix}$ | 47.17 m |
| $\beta$ | 58.0° = $\psi_0$ → **óptimo** |
| $B_\perp$ | 47.17 m |

$$W_z = 2 \times 50\text{ MHz} \times 0.906 + \frac{3\times10^8 \times 47.17 \times 0.848 \times 0.530}{0.03 \times 188.7 \times 2 \times 0.906}$$

$$= 90.6\text{ MHz} + 619.8\text{ MHz} = 710.4\text{ MHz}$$

$$\delta_z = \frac{3\times10^8}{2 \times 710.4\times10^6} = \mathbf{0.2111\text{ m}}$$

| | Predicho | Medido | Error |
|-|:--------:|:------:|:-----:|
| $\delta_z$ | 0.2111 m | **0.2112 m** | **0.05 %** |

---

## 4. Interpretación física

### 4.1 Por qué δx = δy → ∞ con esta configuración

La PSF es la transformada de Fourier del soporte en k-espacio. Con k_x=k_y=0 para todo $(f,t)$, el soporte es un **segmento unidimensional** sobre el eje $k_z$:

$$\mathcal{K} = \{(0,\;0,\;k_z) : k_z \in [k_z^{min},\; k_z^{max}]\}$$

La TF de una distribución sobre una línea es independiente de las otras dos coordenadas → PSF constante en (x,y).

```
k_z
  │████████████████  ← soporte k-espacio (unidimensional)
  │
──┼──────────── k_x
 k_y (⊗)
```

### 4.2 Por qué δz es excelente

La hélice cónica con $\beta = \psi_0 = 58°$ maximiza $B_\perp = B_{helix}$. La contribución tomográfica $\Delta f_z = c B_\perp \sin\psi_0\cos\psi_0/(\lambda_0 R_0 n_2\cos\theta_{t,0}) = 619.8$ MHz amplía el ancho de banda vertical efectivo $W_z$ en 6.8× respecto al BW del chirp solo. Esto da una resolución vertical equivalente a tener un chirp de 710 MHz → δz ≈ 21 cm.

### 4.3 El bistático NorthOffset=180° funciona como SAR tomográfico puro

Esta configuración es equivalente, para targets en el eje, a:
- **Monoestático** con $k_z = -(4\pi f n_2/c)\cos\theta_t$ → misma fórmula de δz
- **SAR tomográfico** donde las diferentes alturas del helicoide crean la apertura en elevación que resuelve en z

No es un sistema con utilidad para imaging 3D completo de la escena. Es óptimo para:
- Detectar y medir la profundidad de objetos enterrados en el eje
- Estudiar la resolución vertical en función de la geometría helicoidal

---

## 5. Limitaciones conocidas del modelo

| Limitación | Efecto | Magnitud |
|-----------|--------|---------|
| Target en eje (x_P=y_P=0 asumido) | Para x_P ≠ 0: cancelación parcial de k_x,k_y → δx finito pero mayor que monoestático | Significativo para x_P > λ₀ |
| Aproximación $\psi_0 = \arctan(\rho_0/z_0)$ ignorando z_P | Subestima ligeramente el look angle real | ~5% para z_P = −5m, z₀=100m |
| Modelo sin pérdidas en suelo | La atenuación dieléctrica no afecta la resolución (solo la SNR) | No afecta δz |
| Criterio de resolución: solo δz | δx, δy son indeterminados con esta geometría (no son una propiedad del sistema sino de la grilla) | — |
| Profundidad arbitraria z_P | La fórmula usa ψ₀ calculado desde la interfaz, no desde z_P | Error ∝ |z_P|/z₀ |

---

## 6. Extensión al caso bistático general (NorthOffset ≠ 180°)

Para NorthOffset $\neq 180°$ con target en el eje, los vectores k_x y k_y **no se cancelan**. En ese caso, la resolución horizontal recupera la forma monoestática ponderada por el ángulo biestático $\beta_P$ en el target:

$$|\mathbf{k}_{xy}|^{max} = \frac{4\pi f_0 n_2}{c}\sin\theta_{t,0}\,|\sin(\Delta\phi/2)|$$

donde $\Delta\phi$ es el ángulo de separación entre TX y RX. Para NorthOffset=180°: $\Delta\phi = \pi$ → $|\sin(\pi/2)| = 1$ → máxima cancelación horizontal (sin resolución).

Para NorthOffset=90°: $|\sin(\pi/4)| = \sqrt{2}/2$ → hay cobertura horizontal parcial.

---

## 7. Resumen final del modelo

> **Para SAR biestático helicoidal cónico con TX/RX diametralmente opuestos (NorthOffset=180°) y target en el eje:**

$$\boxed{
\delta_z = \frac{c}{2\,W_z}, \qquad W_z = n_2 B\cos\theta_{t,0} + \frac{c\,B_\perp\sin\psi_0\cos\psi_0}{\lambda_0 R_0 n_2\cos\theta_{t,0}}
}$$

- **δx = δy**: sin resolución — acotados por la grilla de procesado
- **δz**: completamente determinado por la fórmula anterior
- **Configuración óptima:** β = ψ₀ → B_⊥ = B_helix
- **Validado** para: f₀=10 GHz, B=50 MHz, n₂=2, ρ₀=160m, z₀=100m, B_⊥=47m → **δz=0.211 m (error 0.05%)**
