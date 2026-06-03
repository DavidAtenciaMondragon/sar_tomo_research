# Hipótesis v1 — Resolución 3D para SAR Biestático Helicoidal con Target en el Eje

**Fecha:** 2026-06-01  
**Estado:** Pendiente de verificación por simulación  
**Configuración:** `_espiral` · NorthOffset TX=0°, RX=180° · n₂=2 · Target `[0,0,−5]`  
**Referencia:** `models/resolucion_v2.md`, `models/fase_modelo.md`

---

## 1. Punto de partida: vector de número de onda instantáneo

Del modelo de fase (`models/fase_modelo.md`, §5):

$$\mathbf{k}(f,t) = \frac{2\pi f\,n_2}{c}\!\left(\hat{\mathbf{e}}_2^{TX}(t) + \hat{\mathbf{e}}_2^{RX}(t)\right)$$

Los vectores unitarios apuntan desde los puntos de refracción $\mathbf{Q}_{TX}^*$ y $\mathbf{Q}_{RX}^*$ hasta el objetivo $\mathbf{P}$.

---

## 2. Simetría exacta para target en el eje (x_P = y_P = 0)

### 2.1 Posiciones de TX y RX con NorthOffset=180°

En el instante $t$ con ángulo azimutal $\alpha(t)$:

$$\mathbf{r}_{TX}(t) = \big(\rho(t)\cos\alpha,\; \rho(t)\sin\alpha,\; z(t)\big)$$
$$\mathbf{r}_{RX}(t) = \big(-\rho(t)\cos\alpha,\; -\rho(t)\sin\alpha,\; z(t)\big)$$

(misma altura $z(t)$ y radio $\rho(t)$ en todo instante, separación azimutal exacta de 180°).

### 2.2 Puntos de refracción para target en el eje

Para $\mathbf{P} = (0,0,z_P)$ con $z_P < 0$, por simetría cilíndrica:

$$\mathbf{Q}_{TX}^*(t) = \big(Q_\rho(t)\cos\alpha,\; Q_\rho(t)\sin\alpha,\; 0\big)$$
$$\mathbf{Q}_{RX}^*(t) = \big(-Q_\rho(t)\cos\alpha,\; -Q_\rho(t)\sin\alpha,\; 0\big)$$

donde $Q_\rho(t) > 0$ es el radio del punto de refracción (mismo para TX y RX por simetría).

### 2.3 Vectores unitarios en el Medio 2

$$\hat{\mathbf{e}}_2^{TX}(t) = \frac{\mathbf{P}-\mathbf{Q}_{TX}^*}{|\mathbf{P}-\mathbf{Q}_{TX}^*|} = \frac{(-Q_\rho\cos\alpha,\,-Q_\rho\sin\alpha,\,z_P)}{R_2}$$

$$\hat{\mathbf{e}}_2^{RX}(t) = \frac{\mathbf{P}-\mathbf{Q}_{RX}^*}{|\mathbf{P}-\mathbf{Q}_{RX}^*|} = \frac{(+Q_\rho\cos\alpha,\,+Q_\rho\sin\alpha,\,z_P)}{R_2}$$

donde $R_2 = \sqrt{Q_\rho^2 + z_P^2}$.

### 2.4 Cancelación exacta de componentes horizontales

$$\hat{\mathbf{e}}_2^{TX} + \hat{\mathbf{e}}_2^{RX} = \frac{1}{R_2}\big(0,\; 0,\; 2z_P\big) = \frac{2z_P}{R_2}\hat{z}$$

$$\boxed{\mathbf{k}(f,t) = \frac{2\pi f\,n_2}{c}\cdot\frac{2z_P}{R_2}\,\hat{z} = \left(0,\; 0,\; \frac{4\pi f\,n_2\,z_P}{c\,R_2(t)}\right)}$$

**Consecuencias:**

| Componente | Valor | Significado |
|------------|-------|-------------|
| $k_x$ | **0** para todo $(f,t)$ | Sin cobertura horizontal en x |
| $k_y$ | **0** para todo $(f,t)$ | Sin cobertura horizontal en y |
| $k_z$ | $\neq 0$, varía con $t$ y $f$ | Única fuente de resolución |

---

## 3. Hipótesis de resolución

### H1: Resolución horizontal

$$\boxed{\delta_x = \delta_y \rightarrow \infty \quad \text{(sin resolución para target en eje)}}$$

En la práctica, la medición queda **acotada por el tamaño de la grilla de procesado**:

$$\delta_x^{medido} \approx \delta_y^{medido} \approx 2 \cdot l_{xy} = 2 \times 0.12 = 0.24\text{ m}$$

*Predicción para simulación: el pico en el CSV mostrará $\delta_x \approx \delta_y \approx 0.24$ m (el span completo de la grilla) o un valor limitado por la interpolación en el borde.*

### H2: Resolución vertical

Usando $z_P/R_2 = -\cos\theta_t$ (ángulo de transmisión en Medio 2, con $z_P < 0$):

$$k_z(f,t) = -\frac{4\pi f\,n_2}{c}\cos\theta_t(t)$$

La extensión espectral vertical exacta:

$$\Delta k_z = \frac{4\pi n_2}{c}\!\left[(f_0+\tfrac{B}{2})\cos\theta_t^{top} - (f_0-\tfrac{B}{2})\cos\theta_t^{bot}\right]$$

donde $\theta_t^{top}$ y $\theta_t^{bot}$ son los ángulos de transmisión en la cima y base del helicoide respectivamente, obtenidos por la ley de Snell: $\sin\theta_t = \sin\psi/n_2$.

La resolución vertical (criterio primer nulo):

$$\boxed{\delta_z = \frac{2\pi}{\Delta k_z} = \frac{c}{2\,W_z}}$$

con el ancho de banda efectivo vertical en dos medios:

$$W_z = \underbrace{n_2\,B\,\cos\theta_{t,0}}_{\text{contribución de BW}} + \underbrace{\frac{c\,B_\perp\,\sin\psi_0\cos\psi_0}{\lambda_0\,R_0\,n_2\cos\theta_{t,0}}}_{\text{contribución tomográfica}}$$

---

## 4. Cálculo numérico para la configuración `_espiral`

### 4.1 Parámetros del sistema

| Parámetro | Valor | Fuente |
|-----------|-------|--------|
| $f_0$ | 10 GHz | `radarTx_espiral.json` |
| $\lambda_0 = c/f_0$ | 0.03 m | — |
| $B$ | 50 MHz | FreqMenor=−25 MHz, FreqMayor=+25 MHz |
| $n_2$ | 2 | `GS_snell_test_script.m` línea 82 |
| $\rho_0 = (\rho_{top}+\rho_{bot})/2$ | 160 m | (147.5+172.5)/2 |
| $z_0 = (z_{top}+z_{bot})/2$ | 100 m | (120+80)/2 |
| $R_0 = \sqrt{\rho_0^2+z_0^2}$ | 188.7 m | — |
| $\psi_0 = \arctan(\rho_0/z_0)$ | 58.0° | arctan(1.60) |
| $\cos\theta_{t,0} = \sqrt{1-\sin^2\psi_0/n_2^2}$ | 0.9057 | Snell con n₂=2 |

### 4.2 Parámetros de la hélice cónica

Hélice cónica: radio aumenta de 147.5 m (cima) a 172.5 m (base):

$$B_{helix} = \sqrt{(z_{top}-z_{bot})^2+(\rho_{bot}-\rho_{top})^2} = \sqrt{40^2+25^2} = 47.17\text{ m}$$

$$\beta = \arctan\!\left(\frac{40}{25}\right) = \arctan(1.60) = 58.0°$$

$$B_\perp = B_{helix}\,|\cos(\beta-\psi_0)| = 47.17\,|\cos(58°-58°)| = 47.17\text{ m}$$

> **Observación:** $\beta = \psi_0 = 58°$ → la hélice cónica es **óptima** ($B_\perp = B_{helix}$, máxima apertura tomográfica para la longitud dada).

### 4.3 Cálculo de $W_z$

**Término de ancho de banda:**
$$n_2\,B\,\cos\theta_{t,0} = 2 \times 50\times10^6 \times 0.9057 = 90.57\text{ MHz}$$

**Término tomográfico:**
$$\frac{c\,B_\perp\,\sin\psi_0\cos\psi_0}{\lambda_0\,R_0\,n_2\cos\theta_{t,0}} = \frac{3\times10^8 \times 47.17 \times 0.8480 \times 0.5299}{0.03 \times 188.7 \times 2 \times 0.9057} = \frac{6.36\times10^9}{10.26} = 619.8\text{ MHz}$$

$$W_z = 90.57 + 619.8 = 710.4\text{ MHz}$$

> La contribución tomográfica domina sobre la de ancho de banda: razón ≈ 6.8:1. Esto es consecuencia de $f_0=10$ GHz (λ pequeña → alta sensibilidad angular) combinado con la apertura geométrica de 47 m.

### 4.4 Predicción de resolución vertical

$$\boxed{\delta_z^{pred} = \frac{c}{2\,W_z} = \frac{3\times10^8}{2 \times 7.104\times10^8} \approx 0.211\text{ m} = 21.1\text{ cm}}$$

**Verificación cruzada por k-espacio directo:**

Ángulos de transmisión en los extremos (usando $\psi = \arctan(\rho/z)$):

| Posición | $\rho$ | $z$ | $\psi$ | $\sin\theta_t = \sin\psi/2$ | $\cos\theta_t$ |
|---------|--------|-----|--------|--------------------------|---------------|
| Cima (top) | 147.5 | 120 | 50.87° | 0.388 | 0.9216 |
| Base (bot) | 172.5 | 80 | 65.10° | 0.453 | 0.8911 |

$$\Delta k_z = \frac{4\pi\times2}{3\times10^8}\!\left[10.025\times10^9\times0.9216 - 9.975\times10^9\times0.8911\right]$$
$$= \frac{8\pi}{3\times10^8}\!\left[9.234\times10^9 - 8.888\times10^9\right] = \frac{8\pi}{3\times10^8}\times 3.46\times10^8 = 8\pi\times1.153 = 29.0\text{ m}^{-1}$$

$$\delta_z^{kspace} = \frac{2\pi}{29.0} = 0.217\text{ m} = 21.7\text{ cm}$$

Ambos métodos concuerdan dentro del 3%.

---

## 5. Tabla de predicciones

| Dirección | $\Delta k$ predicho | $\delta$ predicho | Criterio | Nota |
|-----------|--------------------|--------------------|---------|------|
| **X** | $\Delta k_x = 0$ | **∞** (→ 0.24 m en grilla) | — | Sin cobertura horizontal |
| **Y** | $\Delta k_y = 0$ | **∞** (→ 0.24 m en grilla) | — | Sin cobertura horizontal |
| **Z** | $\Delta k_z \approx 29.0\text{ m}^{-1}$ | **≈ 0.211–0.217 m** | Primer nulo | Dominada por apertura tomográfica |

*Nota: el criterio "primer nulo" da $\delta = 2\pi/\Delta k$. El criterio −3 dB da un valor ~20% menor: $\delta_{-3dB} \approx 0.167\text{ m}$ para PSF gaussiana o $\approx 0.174\text{ m}$ para PSF sinc.*

---

## 6. Hipótesis de forma cerrada final

$$\boxed{
\begin{aligned}
\delta_x &= \delta_y = \text{sin resolución para target en eje biestático} \\[6pt]
\delta_z &= \frac{c}{2\,W_z} \\[4pt]
W_z &= n_2 B\cos\theta_{t,0} + \frac{c\,B_\perp\,\sin\psi_0\cos\psi_0}{\lambda_0\,R_0\,n_2\cos\theta_{t,0}} \\[4pt]
B_\perp &= B_{helix}|\cos(\beta-\psi_0)| \\[4pt]
\cos\theta_{t,0} &= \sqrt{1-\sin^2\psi_0/n_2^2} \\[4pt]
\psi_0 &= \arctan(\rho_0/z_0), \quad R_0 = \sqrt{\rho_0^2+z_0^2}
\end{aligned}
}$$

**Predicción numérica:** $\delta_z \approx 0.211$ m (21.1 cm)

**Predicción para CSV de simulación:**  
- `Resolution_m` para Z: ≈ **0.21 m** (tolerancia ±20% por aproximaciones de look angle y criterio de resolución)  
- `Resolution_m` para X e Y: ≈ **0.24 m** (limitado por grilla, no por el sistema)

---

## 7. Posibles fuentes de error en la hipótesis

1. **Aproximación del look angle:** se usó $\psi_0 = \arctan(\rho_0/z_0)$ con la altura media del sensor. El look angle efectivo varía a lo largo de la trayectoria y la profundidad del target ($z_P=-5$m) introduce una pequeña corrección que se ignoró.

2. **Criterio de resolución:** la fórmula usa el criterio de primer nulo ($2\pi/\Delta k$). El PROC mide el criterio −3 dB, que da valores ~20% menores.

3. **La grilla dz=0.08 m con 7 puntos** da una extensión total de ±0.24m desde el target. Si $\delta_z \approx 0.21$ m, el lóbulo principal se extiende ±0.105m desde el pico, lo cual está dentro de los ±0.24m de la grilla. La medición debería capturar el lóbulo principal correctamente.

4. **El modelo asume target exactamente en el eje.** El target está en [0,0,-5], que es el eje del helicoide (x=0, y=0). La simetría es exacta.
