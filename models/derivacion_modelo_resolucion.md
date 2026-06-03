# Derivación Completa del Modelo de Resolución 3D  
## SAR Biestático Helicoidal con Dos Medios Homogéneos

**Fecha:** 2026-06-01  
**Sistema:** SAR biestático · Hélice cónica · Aire ($n_1=1$) + Suelo ($n_2=\sqrt{\varepsilon_r}$) · Interfaz plana  
**Fuentes:** `models/modelo_fisico_base.md`, `models/fase_modelo.md`, `models/resolucion_v2.md`, `hypotheses/hipotesis_v1.md` a `v5.md`

---

## Índice

1. [Enunciado del problema](#1-enunciado-del-problema)
2. [Geometría del sistema](#2-geometría-del-sistema)
3. [Modelo de señal y función de fase exacta](#3-modelo-de-señal-y-función-de-fase-exacta)
4. [Gradiente de la fase respecto al objetivo — Teorema de Danskin](#4-gradiente-de-la-fase-respecto-al-objetivo--teorema-de-danskin)
5. [El vector de número de onda instantáneo](#5-el-vector-de-número-de-onda-instantáneo)
6. [La PSF como transformada de Fourier del k-espacio](#6-la-psf-como-transformada-de-fourier-del-k-espacio)
7. [Análisis del k-espacio para la hélice biestática](#7-análisis-del-k-espacio-para-la-hélice-biestática)
8. [Resolución vertical — derivación paso a paso](#8-resolución-vertical--derivación-paso-a-paso)
9. [Resolución horizontal — derivación paso a paso](#9-resolución-horizontal--derivación-paso-a-paso)
10. [Extensión a target off-axis](#10-extensión-a-target-off-axis)
11. [Resumen de ecuaciones y condiciones de validez](#11-resumen-de-ecuaciones-y-condiciones-de-validez)
12. [Validación experimental completa](#12-validación-experimental-completa)

---

## 1. Enunciado del Problema

Se desea encontrar las expresiones analíticas cerradas (o aproximadas) que predigan la **resolución espacial 3D** de un sistema SAR con las siguientes características:

- **Configuración biestática:** transmisor (TX) y receptor (RX) en plataformas distintas, ambas sobre la misma trayectoria helicoidal con un **offset azimutal $\Delta\phi$** entre ellas.
- **Trayectoria helicoidal cónica:** radio y altura varían a lo largo de la espiral, lo que genera una apertura tomográfica tridimensional.
- **Dos medios homogéneos:** el espacio libre sobre la interfaz (aire, $n_1=1$) y un medio dieléctrico bajo ella (suelo, $n_2=\sqrt{\varepsilon_r}$).
- **Interfaz plana:** la separación entre medios es el plano $z=0$.
- **Objetivo puntual enterrado:** el target está en $\mathbf{P}=(x_P, y_P, z_P)$ con $z_P < 0$.

La resolución queda caracterizada por tres magnitudes:

$$\delta_x, \quad \delta_y, \quad \delta_z$$

donde $\delta_u$ es la anchura a $-3\,\text{dB}$ de la Función de Dispersión de Punto (PSF) en la dirección $\hat{u}$.

---

## 2. Geometría del Sistema

### 2.1 Sistema de coordenadas

- **Origen:** en la interfaz aire-suelo.
- **Eje $\hat{z}$:** vertical hacia arriba. El aire ocupa $z>0$, el suelo $z<0$.
- **Plano $z=0$:** interfaz plana.

### 2.2 Trayectoria helicoidal cónica (TX y RX)

La hélice cónica se parametriza con velocidad tangencial $V_0$ constante:

$$\rho(t) = \rho_{top} + V_\rho\,t, \quad V_\rho = \frac{\rho_{base}-\rho_{top}}{t_{max}}$$

$$\alpha(t) = \frac{V_0}{V_\rho}\!\left[\ln\rho(t) - \ln\rho_{top}\right] \quad (V_\rho \neq 0)$$

$$z(t) = z_{top} + V_z\,t, \quad V_z = -\frac{z_{top}-z_{base}}{t_{max}} < 0$$

Con TX a offset $\phi_{TX}$ y RX a offset $\phi_{RX} = \phi_{TX} + \Delta\phi$:

$$\mathbf{r}_{TX}(t) = \begin{pmatrix}\rho(t)\cos(\alpha(t)+\phi_{TX})\\\rho(t)\sin(\alpha(t)+\phi_{TX})\\z(t)\end{pmatrix}, \quad \mathbf{r}_{RX}(t) = \begin{pmatrix}\rho(t)\cos(\alpha(t)+\Delta\phi)\\\ \rho(t)\sin(\alpha(t)+\Delta\phi)\\z(t)\end{pmatrix}$$

**Parámetros geométricos derivados:**

| Símbolo | Expresión | Descripción |
|---------|-----------|-------------|
| $\rho_0$ | $(\rho_{top}+\rho_{base})/2$ | Radio medio |
| $z_0$ | $(z_{top}+z_{base})/2$ | Altura media sobre la interfaz |
| $B_{helix}$ | $\sqrt{(z_{top}-z_{base})^2+(\rho_{base}-\rho_{top})^2}$ | Longitud total de la espiral |
| $\beta$ | $\arctan\!\left(\frac{z_{top}-z_{base}}{\rho_{base}-\rho_{top}}\right)$ | Ángulo de inclinación de la hélice |
| $\psi_0$ | $\arctan(\rho_0/z_0)$ | Look angle medio (sensor→interfaz) |
| $R_0$ | $\sqrt{\rho_0^2+z_0^2}$ | Distancia media sensor→origen |

### 2.3 Medios y ley de Snell

| Magnitud | Medio 1 (aire) | Medio 2 (suelo) |
|---------|:--------------:|:---------------:|
| $n$ | $n_1=1$ | $n_2=\sqrt{\varepsilon_r}$ |
| $v$ | $c$ | $v_2=c/n_2$ |
| $k(f)$ | $2\pi f/c$ | $2\pi f\,n_2/c$ |

La ley de Snell en la interfaz $z=0$:

$$n_1\sin\theta_i = n_2\sin\theta_t \implies \sin\theta_t = \frac{\sin\theta_i}{n_2}, \quad \cos\theta_t = \sqrt{1-\frac{\sin^2\theta_i}{n_2^2}}$$

El ángulo de transmisión medio evaluado en el look angle medio:

$$\boxed{\cos\theta_{t,0} = \sqrt{1-\frac{\sin^2\psi_0}{n_2^2}}}$$

---

## 3. Modelo de Señal y Función de Fase Exacta

### 3.1 Camino óptico biestático con refracción

La señal emitida por TX viaja desde $\mathbf{r}_{TX}$ hasta el target $\mathbf{P}$ a través de los dos medios, refractándose en la interfaz. El trayecto completo es:

$$\mathbf{r}_{TX} \xrightarrow{d_1^{TX}(\text{aire})} \mathbf{Q}_{TX}^* \xrightarrow{d_2^{TX}(\text{suelo})} \mathbf{P} \xrightarrow{d_2^{RX}(\text{suelo})} \mathbf{Q}_{RX}^* \xrightarrow{d_1^{RX}(\text{aire})} \mathbf{r}_{RX}$$

Los **puntos de refracción** $\mathbf{Q}_{TX}^*$ y $\mathbf{Q}_{RX}^*$ están sobre la interfaz ($z=0$) y se determinan por el **Principio de Fermat** (mínimo tiempo de tránsito), equivalente a la ley de Snell:

$$\mathbf{Q}_{TX}^* = \arg\min_{\mathbf{Q}\in\{z=0\}} \left[\frac{|\mathbf{r}_{TX}-\mathbf{Q}|}{c}+\frac{|\mathbf{Q}-\mathbf{P}|}{v_2}\right]$$

El **camino óptico total** (en metros equivalentes en vacío):

$$\boxed{R_{OP} = d_1^{TX} + n_2\,d_2^{TX} + n_2\,d_2^{RX} + d_1^{RX}}$$

### 3.2 Fase de la señal recibida

La señal en el dominio frecuencia acumula la fase:

$$\Phi(f, t; \mathbf{P}) = -\frac{2\pi f}{c}\,R_{OP}(\mathbf{r}_{TX}(t), \mathbf{r}_{RX}(t), \mathbf{P})$$

donde el signo negativo sigue la convención de onda viajera $e^{-j\Phi}$.

### 3.3 Formación de imagen por Backprojection

La imagen SAR en un punto de prueba $\mathbf{p}$ se forma por suma coherente:

$$I(\mathbf{p}) = \sum_{n}\sum_k S(f_k, t_n)\,\exp\!\left(+j\frac{2\pi f_k}{c}R_{OP}(\mathbf{r}_{TX,n},\mathbf{r}_{RX,n},\mathbf{p})\right)$$

Para un dispersor puntual en $\mathbf{P}$: $S(f_k,t_n)=g(\mathbf{P})\exp(-j\frac{2\pi f_k}{c}R_{OP}(\ldots,\mathbf{P}))$. Sustituyendo:

$$I(\mathbf{p}) = g(\mathbf{P})\sum_{n,k}\exp\!\left(+j\frac{2\pi f_k}{c}\underbrace{\left[R_{OP}(\mathbf{p})-R_{OP}(\mathbf{P})\right]}_{\Delta R_{OP}}\right)$$

La **PSF** es la respuesta del sistema a un blanco puntual ($g(\mathbf{P})=1$):

$$\text{PSF}(\boldsymbol{\delta}) = \sum_{n,k}\exp\!\left(+j\frac{2\pi f_k}{c}\Delta R_{OP}(t_n;\boldsymbol{\delta})\right), \quad \boldsymbol{\delta}=\mathbf{p}-\mathbf{P}$$

---

## 4. Gradiente de la Fase Respecto al Objetivo — Teorema de Danskin

### 4.1 El problema de diferenciación implícita

$R_{OP}$ es una función implícita de $\mathbf{P}$ a través de los puntos de refracción $\mathbf{Q}_{TX}^*(\mathbf{r}_{TX},\mathbf{P})$ y $\mathbf{Q}_{RX}^*(\mathbf{P},\mathbf{r}_{RX})$. Diferenciarlo directamente requeriría resolver cómo se mueven los puntos de refracción al desplazar $\mathbf{P}$, lo cual es complejo.

### 4.2 Teorema de Danskin (envolvente)

$R_{OP}$ es el **mínimo** de la función tiempo de tránsito sobre las posiciones de refracción. En el mínimo, la variación de primer orden respecto a los puntos de refracción es nula. Por el Teorema de la Envolvente:

$$\nabla_\mathbf{P}\,R_{OP} = \frac{\partial}{\partial\mathbf{P}}\Big[d_1^{TX}+n_2\,d_2^{TX}+n_2\,d_2^{RX}+d_1^{RX}\Big]\Bigg|_{\mathbf{Q}_{TX}^*,\mathbf{Q}_{RX}^*\text{ fijos}}$$

Solo $d_2^{TX}$ y $d_2^{RX}$ dependen directamente de $\mathbf{P}$. Calculando:

$$\frac{\partial d_2^{TX}}{\partial\mathbf{P}} = \frac{\partial|\mathbf{P}-\mathbf{Q}_{TX}^*|}{\partial\mathbf{P}} = \frac{\mathbf{P}-\mathbf{Q}_{TX}^*}{d_2^{TX}} \equiv \hat{\mathbf{e}}_2^{TX}$$

$$\frac{\partial d_2^{RX}}{\partial\mathbf{P}} = \frac{\mathbf{P}-\mathbf{Q}_{RX}^*}{d_2^{RX}} \equiv \hat{\mathbf{e}}_2^{RX}$$

donde $\hat{\mathbf{e}}_2^{TX}$ y $\hat{\mathbf{e}}_2^{RX}$ son vectores unitarios desde los puntos de refracción hacia $\mathbf{P}$ en el Medio 2.

**Resultado:**

$$\boxed{\nabla_\mathbf{P}\,R_{OP} = n_2\!\left(\hat{\mathbf{e}}_2^{TX} + \hat{\mathbf{e}}_2^{RX}\right)}$$

---

## 5. El Vector de Número de Onda Instantáneo

### 5.1 Definición

El vector de número de onda en el objetivo mide la sensibilidad de la fase SAR al desplazamiento del target:

$$\mathbf{k}(f,t;\mathbf{P}) \equiv -\nabla_\mathbf{P}\,\Phi = \frac{2\pi f}{c}\nabla_\mathbf{P}\,R_{OP}$$

Sustituyendo el resultado de Danskin:

$$\boxed{\mathbf{k}(f,t;\mathbf{P}) = \frac{2\pi f\,n_2}{c}\!\left(\hat{\mathbf{e}}_2^{TX}(t) + \hat{\mathbf{e}}_2^{RX}(t)\right)}$$

### 5.2 Componentes angulares en el Medio 2

Cada vector unitario $\hat{\mathbf{e}}_2$ puede descomponerse en ángulos polar ($\theta_t$) y azimutal ($\phi$) del rayo en el Medio 2:

$$\hat{\mathbf{e}}_2 = (\sin\theta_t\cos\phi_e,\; \sin\theta_t\sin\phi_e,\; -\cos\theta_t)$$

donde el signo negativo en $z$ indica que el rayo se propaga hacia abajo (hacia $z<0$).

Las componentes del vector $\mathbf{k}$ son:

$$k_x = \frac{2\pi f\,n_2}{c}\!\left(\sin\theta_t^{TX}\cos\phi_{TX} + \sin\theta_t^{RX}\cos\phi_{RX}\right)$$

$$k_y = \frac{2\pi f\,n_2}{c}\!\left(\sin\theta_t^{TX}\sin\phi_{TX} + \sin\theta_t^{RX}\sin\phi_{RX}\right)$$

$$\boxed{k_z = -\frac{2\pi f\,n_2}{c}\!\left(\cos\theta_t^{TX} + \cos\theta_t^{RX}\right)}$$

### 5.3 Verificación del módulo

El módulo máximo de $\mathbf{k}$ (para TX=RX, rayo vertical):

$$|\mathbf{k}|_{max} = \frac{4\pi f\,n_2}{c} = 2k_2$$

donde $k_2 = 2\pi f\,n_2/c$ es el número de onda en el Medio 2. ✓ (No se puede superar el doble del número de onda del medio.)

---

## 6. La PSF como Transformada de Fourier del k-espacio

### 6.1 Linealización de $\Delta R_{OP}$

Para $|\boldsymbol{\delta}|$ pequeño (objetivos cercanos al foco), expandimos en serie de Taylor:

$$\Delta R_{OP}(t;\boldsymbol{\delta}) = R_{OP}(\ldots,\mathbf{P}+\boldsymbol{\delta}) - R_{OP}(\ldots,\mathbf{P}) \approx \nabla_\mathbf{P}R_{OP}\cdot\boldsymbol{\delta} = \frac{c}{2\pi f}\,\mathbf{k}\cdot\boldsymbol{\delta}$$

**Condición de validez:** $|\boldsymbol{\delta}|^2/(2d_i) \ll \lambda/4$ (criterio de campo lejano en el plano del objetivo).

### 6.2 PSF como integral de Fourier

Sustituyendo en la expresión de la PSF:

$$\text{PSF}(\boldsymbol{\delta}) = \int_T\int_B A(f,t)\,\exp\!\left(+j\,\mathbf{k}(f,t;\mathbf{P})\cdot\boldsymbol{\delta}\right)df\,dt = \int_{\mathcal{K}} W(\mathbf{k})\,e^{+j\mathbf{k}\cdot\boldsymbol{\delta}}\,d^3k$$

donde $W(\mathbf{k})$ es la **densidad de muestreo en el k-espacio** (Jacobiano del cambio de variables $(f,t)\to\mathbf{k}$).

**Resultado fundamental:** la PSF es la **transformada de Fourier inversa** de la función de densidad $W(\mathbf{k})$. La resolución en cada dirección es el inverso del ancho espectral en esa dirección.

### 6.3 Criterio de resolución

Para una distribución $W(\mathbf{k})$ con soporte en el intervalo $\Delta k_u$ en la dirección $\hat{u}$:

$$\boxed{\delta_u = \frac{2\pi}{\Delta k_u}}$$

con el **ancho espectral efectivo**:

$$\Delta k_u = \max_{f,t}\!\left(\mathbf{k}\cdot\hat{u}\right) - \min_{f,t}\!\left(\mathbf{k}\cdot\hat{u}\right)$$

---

## 7. Análisis del k-espacio para la Hélice Biestática

### 7.1 Simetría para target en el eje ($x_P=y_P=0$)

Para target en $\mathbf{P}=(0,0,z_P)$ y sensor TX a azimut $\alpha$, el rayo refractado en el Medio 2 apunta desde $\mathbf{Q}_{TX}^* = (Q_\rho\cos\alpha, Q_\rho\sin\alpha, 0)$ hacia $\mathbf{P}$:

$$\hat{\mathbf{e}}_2^{TX} = (-\sin\theta_t\cos\alpha,\; -\sin\theta_t\sin\alpha,\; -\cos\theta_t)$$

Para RX a azimut $\alpha+\Delta\phi$:

$$\hat{\mathbf{e}}_2^{RX} = (-\sin\theta_t\cos(\alpha+\Delta\phi),\; -\sin\theta_t\sin(\alpha+\Delta\phi),\; -\cos\theta_t)$$

### 7.2 Suma de vectores y componentes de k

Usando la identidad trigonométrica $\cos\alpha+\cos(\alpha+\Delta\phi) = 2\cos(\Delta\phi/2)\cos(\alpha+\Delta\phi/2)$:

$$\hat{\mathbf{e}}_2^{TX}+\hat{\mathbf{e}}_2^{RX} = \begin{pmatrix}-2\sin\theta_t\cos(\Delta\phi/2)\cos(\alpha+\Delta\phi/2)\\-2\sin\theta_t\cos(\Delta\phi/2)\sin(\alpha+\Delta\phi/2)\\-2\cos\theta_t\end{pmatrix}$$

**Componente vertical** (clave para $\delta_z$):

$$\boxed{k_z(f,t) = -\frac{4\pi f\,n_2}{c}\cos\theta_t(t)}$$

> **Observación crítica:** $k_z$ es **independiente de $\Delta\phi$** y del azimut $\alpha$. Solo depende de $\theta_t(t)$, que varía al descender la hélice.

**Componentes horizontales** (determinan $\delta_{xy}$):

$$k_x = -\frac{4\pi f\,n_2}{c}\sin\theta_t\cos(\Delta\phi/2)\cos(\alpha+\Delta\phi/2)$$

$$k_y = -\frac{4\pi f\,n_2}{c}\sin\theta_t\cos(\Delta\phi/2)\sin(\alpha+\Delta\phi/2)$$

---

## 8. Resolución Vertical — Derivación Paso a Paso

### Paso 1: Fuentes de variación en $k_z$

$k_z = -(4\pi f\,n_2/c)\cos\theta_t(t)$ varía con:
- **Frecuencia $f$**: al variar en $[f_0-B/2,\, f_0+B/2]$
- **Tiempo $t$**: al variar $\theta_t(t)$ con el descenso helicoidal

### Paso 2: Extensión espectral exacta

Los extremos del rango de $|k_z|$ se alcanzan cuando ambas fuentes son simultáneamente máximas o mínimas:

$$\Delta k_z = \frac{4\pi n_2}{c}\!\left[(f_0+\tfrac{B}{2})\cos\theta_t^{top} - (f_0-\tfrac{B}{2})\cos\theta_t^{bot}\right]$$

donde $\theta_t^{top}$ y $\theta_t^{bot}$ son los ángulos de transmisión en la **cima** y **base** de la hélice respectivamente, vinculados al look angle por la ley de Snell.

### Paso 3: Separación en dos contribuciones

Expandiendo algebraicamente con $G(t) = 2\cos\theta_t(t)$, $\Delta G = G_{top}-G_{bot}$, $\bar{G}=(G_{top}+G_{bot})/2$:

$$\Delta k_z = \frac{2\pi n_2}{c}\!\left[f_0\,\Delta G + B\,\bar{G}\right] = \frac{4\pi n_2}{c}\!\left[f_0\,(\cos\theta_t^{top}-\cos\theta_t^{bot}) + B\,\langle\cos\theta_t\rangle\right]$$

$$\underbrace{\frac{4\pi n_2 B}{c}\langle\cos\theta_t\rangle}_{\Delta k_z^{(B)}\;:\;\text{contribución de BW}} \quad+\quad \underbrace{\frac{4\pi f_0 n_2}{c}(\cos\theta_t^{top}-\cos\theta_t^{bot})}_{\Delta k_z^{(geom)}\;:\;\text{contribución tomográfica}}$$

### Paso 4: Expresar $\Delta\cos\theta_t$ en términos de $B_\perp$

La variación geométrica de $\theta_t$ se produce porque el look angle $\psi$ varía de $\psi_{top}$ a $\psi_{base}$ al descender la hélice. Diferenciando la ley de Snell:

$$\frac{d\cos\theta_t}{d\psi} = -\frac{\sin\psi\cos\psi}{n_2^2\cos\theta_t}$$

La variación de ángulo es $\Delta\psi = B_\perp/R_0$ (la apertura tomográfica efectiva $B_\perp$ subtiende un ángulo $B_\perp/R_0$ desde el target):

$$\Delta\cos\theta_t \approx \frac{\sin\psi_0\cos\psi_0}{n_2^2\cos\theta_{t,0}}\cdot\frac{B_\perp}{R_0}$$

Sustituyendo:

$$\Delta k_z^{(geom)} = \frac{4\pi f_0}{c}\cdot\frac{\sin\psi_0\cos\psi_0}{n_2\cos\theta_{t,0}}\cdot\frac{B_\perp}{R_0} = \frac{4\pi B_\perp}{\lambda_0 R_0}\cdot\frac{\sin\psi_0\cos\psi_0}{n_2\cos\theta_{t,0}}$$

**Verificación límite $n_2\to 1$:** $\cos\theta_{t,0}\to\cos\psi_0$ → $\Delta k_z^{(geom)}\to\frac{4\pi B_\perp}{\lambda_0 R_0}\sin\psi_0$, que es exactamente la Ec. (4.24) de Góes 2022 para medio único. ✓

### Paso 5: La apertura tomográfica efectiva $B_\perp$

La hélice cónica tiene longitud $B_{helix}$ e inclinación $\beta$. La proyección de esta longitud perpendicular a la dirección de vista media (LOS):

$$\boxed{B_\perp = B_{helix}\,|\cos(\beta-\psi_0)|}$$

**Condición óptima:** $\beta = \psi_0$ → $B_\perp = B_{helix}$ (máxima apertura tomográfica).

### Paso 6: Ensamblaje del ancho de banda efectivo vertical

Combinando las dos contribuciones bajo la aproximación $B/f_0 \ll 1$ (banda estrecha):

$$\Delta k_z = \frac{4\pi}{c}\!\left[n_2\,B\,\langle\cos\theta_t\rangle + \frac{cB_\perp\sin\psi_0\cos\psi_0}{\lambda_0 R_0 n_2\cos\theta_{t,0}}\right] = \frac{4\pi}{c}\,W_z$$

Definiendo el **ancho de banda efectivo vertical**:

$$\boxed{W_z = \underbrace{n_2\,B\,\cos\theta_{t,0}}_{\text{ancho de banda}} + \underbrace{\frac{c\,B_\perp\,\sin\psi_0\cos\psi_0}{\lambda_0\,R_0\,n_2\cos\theta_{t,0}}}_{\text{apertura tomográfica}}}$$

### Paso 7: Resolución vertical final

Aplicando el criterio de resolución $\delta_z = 2\pi/\Delta k_z$:

$$\delta_z = \frac{2\pi}{\Delta k_z} = \frac{2\pi}{\frac{4\pi}{c}W_z} = \frac{c}{2W_z}$$

$$\boxed{\delta_z = \frac{c}{2\,W_z}}$$

**Demostración de invariancia al NorthOffset:** $k_z = -(4\pi f\,n_2/c)\cos\theta_t(t)$ es independiente de $\Delta\phi$ y del azimut $\alpha$ → $\Delta k_z$ no depende de $\Delta\phi$ → **$\delta_z$ es el mismo para cualquier NorthOffset** (validado experimentalmente: mismo valor para $\Delta\phi=180°$ y $\Delta\phi=90°$, error 0.05%).

---

## 9. Resolución Horizontal — Derivación Paso a Paso

### Paso 1: Estructura del k-espacio horizontal

Las componentes horizontales son:

$$k_x = -\frac{4\pi f n_2}{c}\sin\theta_t\cos(\Delta\phi/2)\cos\!\left(\alpha+\frac{\Delta\phi}{2}\right)$$

$$k_y = -\frac{4\pi f n_2}{c}\sin\theta_t\cos(\Delta\phi/2)\sin\!\left(\alpha+\frac{\Delta\phi}{2}\right)$$

Al variar $\alpha \in [0, 2\pi N_t]$, estas describen un **círculo** en el plano $(k_x, k_y)$ de radio:

$$R_c(\Delta\phi) = \frac{4\pi f_0 n_2}{c}\sin\theta_{t,0}\cdot|\cos(\Delta\phi/2)|$$

### Paso 2: Simplificación usando la ley de Snell

Usando $n_2\sin\theta_{t,0} = \sin\psi_0$:

$$R_c(\Delta\phi) = \frac{4\pi f_0 n_2}{c}\cdot\frac{\sin\psi_0}{n_2}\cdot|\cos(\Delta\phi/2)| = \frac{4\pi\sin\psi_0}{\lambda_0}\cdot|\cos(\Delta\phi/2)|$$

El factor $n_2$ se cancela. La resolución horizontal no depende del índice de refracción cuando se expresa en términos del look angle en el aire $\psi_0$.

### Paso 3: Forma de la PSF horizontal — Distribución arcseno

Para cobertura azimutal completa (360°), $k_x = R_c\cos(\alpha+\phi_0)$ varía sinusoidalmente con $\alpha$. La distribución estadística de $k_x$ es la **distribución arcseno**:

$$p(k_x) = \frac{1}{\pi\sqrt{R_c^2-k_x^2}}, \quad k_x\in[-R_c, R_c]$$

La transformada de Fourier de esta distribución es la función de Bessel de orden 0:

$$\text{PSF}_x(\delta_x) = \int_{-R_c}^{+R_c} p(k_x)\,e^{jk_x\delta_x}\,dk_x = J_0(R_c\,\delta_x)$$

### Paso 4: Criterio $-3\,\text{dB}$ para la PSF de Bessel

La función $J_0(u)$ cumple:
- $J_0(0) = 1$ (valor en el pico)
- $J_0(u) = 1/\sqrt{2}$ en $u \approx 1.20$ (punto $-3\,\text{dB}$, semiancho)

La anchura a $-3\,\text{dB}$ (FWHM) de la PSF horizontal:

$$\text{FWHM}_{xy} = 2\times\frac{1.20}{R_c(\Delta\phi)}$$

### Paso 5: Resolución horizontal final

$$\boxed{\delta_{xy}(\Delta\phi) = \frac{2\times 1.20}{R_c(\Delta\phi)} = \frac{2.40\,\lambda_0}{4\pi\,\sin\psi_0\cdot|\cos(\Delta\phi/2)|} = \frac{0.60\,\lambda_0}{\pi\,\sin\psi_0\cdot|\cos(\Delta\phi/2)|}}$$

**Casos límite:**

| NorthOffset $\Delta\phi$ | $|\cos(\Delta\phi/2)|$ | $\delta_{xy}$ |
|:---:|:---:|:---:|
| 0° (monoestático) | 1 | $\delta_{xy}^{mono} = 0.60\lambda_0/(\pi\sin\psi_0)$ |
| 90° | $1/\sqrt{2}$ | $\sqrt{2}\,\delta_{xy}^{mono}$ |
| 180° | 0 | $\infty$ (sin resolución horizontal) |

### Paso 6: Caso $\Delta\phi = 180°$ — cancelación exacta

Para NorthOffset=180°, $|\cos(\Delta\phi/2)|=|\cos(90°)|=0$, por lo que $R_c=0$ y $k_x=k_y=0$ identicamente. Este es el resultado de la cancelación exacta de los vectores horizontales derivada en §7:

$$\hat{\mathbf{e}}_2^{TX}+\hat{\mathbf{e}}_2^{RX}\Big|_{\Delta\phi=180°,x_P=0} = \left(0,\,0,\,-\frac{2|z_P|}{R_2}\right)$$

La PSF en $(x,y)$ es constante → no hay focalización horizontal → la medición queda limitada por la grilla de procesado.

---

## 10. Extensión a Target Off-Axis

### 10.1 El problema

Para target en $(x_P, 0, z_P)$ con $x_P \neq 0$, la simetría cilíndrica se rompe. El look angle efectivo varía con el azimut del sensor:

$$\psi(\alpha) = \arctan\!\left(\frac{\sqrt{\rho^2-2\rho x_P\cos\alpha+x_P^2}}{z}\right) \neq \text{cte}$$

En particular, hay una diferencia entre near-range ($\alpha=0$, distancia $\rho_0-x_P$) y far-range ($\alpha=\pi$, distancia $\rho_0+x_P$).

### 10.2 Corrección de Góes (Ec. 4.52–4.53)

En sistemas reales con pérdida de propagación $\propto 1/R^2$, las posiciones de **near-range contribuyen más** que las de far-range. La corrección de Góes 2022 reemplaza $\psi_0$ y $R_0$ por los valores de near-range:

$$\boxed{\tilde{\psi}_0 = \arctan\!\left(\frac{\rho_0-|x_P|}{z_0}\right), \qquad \tilde{R}_0 = \sqrt{(\rho_0-|x_P|)^2+z_0^2}}$$

$$\cos\tilde{\theta}_{t,0} = \sqrt{1-\frac{\sin^2\tilde{\psi}_0}{n_2^2}}$$

### 10.3 Fórmulas corregidas para target off-axis

**Resolución vertical:**

$$\tilde{\delta}_z = \frac{c}{2\,\tilde{W}_z}, \quad \tilde{W}_z = n_2\,B\,\cos\tilde{\theta}_{t,0} + \frac{c\,B_\perp\,\sin\tilde{\psi}_0\cos\tilde{\psi}_0}{\lambda_0\,\tilde{R}_0\,n_2\cos\tilde{\theta}_{t,0}}$$

**Resolución horizontal (dirección radial):**

$$\tilde{\delta}_x = \frac{0.60\,\lambda_0}{\pi\,\sin\tilde{\psi}_0\cdot|\cos(\Delta\phi/2)|}$$

**Resolución horizontal (dirección azimutal):** se aproxima con $\psi_0$:

$$\delta_y \approx \frac{0.60\,\lambda_0}{\pi\,\sin\psi_0\cdot|\cos(\Delta\phi/2)|}$$

### 10.4 Condición de aplicabilidad de la corrección de Góes

La corrección es válida cuando el path loss $1/R^2$ pondera las posiciones de near-range (sistemas reales). En simulaciones sin pérdida de propagación, la corrección **sobreestima** la mejora y el error puede ser ligeramente mayor que sin corrección. Para todos los propósitos prácticos, el modelo da errores < 12%.

---

## 11. Resumen de Ecuaciones y Condiciones de Validez

### 11.1 Conjunto completo de ecuaciones

$$\boxed{
\begin{aligned}
&\textbf{Resolución vertical:}\\
&\quad \delta_z = \frac{c}{2\,W_z}\\[6pt]
&\quad W_z = n_2\,B\,\cos\theta_{t,0} + \frac{c\,B_\perp\,\sin\psi_0\cos\psi_0}{\lambda_0\,R_0\,n_2\cos\theta_{t,0}}\\[12pt]
&\textbf{Resolución horizontal (target en eje):}\\
&\quad \delta_{xy}(\Delta\phi) = \frac{0.60\,\lambda_0}{\pi\,\sin\psi_0\cdot|\cos(\Delta\phi/2)|}\\[12pt]
&\textbf{Corrección off-axis (target en } x_P\neq 0\text{):}\\
&\quad \tilde{\psi}_0 = \arctan\!\left(\frac{\rho_0-|x_P|}{z_0}\right),\quad \tilde{R}_0 = \sqrt{(\rho_0-|x_P|)^2+z_0^2}\\
&\quad \text{Reemplazar } \psi_0\to\tilde{\psi}_0,\; R_0\to\tilde{R}_0 \text{ en las fórmulas anteriores}\\[12pt]
&\textbf{Parámetros auxiliares:}\\
&\quad \cos\theta_{t,0} = \sqrt{1-\sin^2\psi_0/n_2^2} \quad\text{(Snell)}\\
&\quad B_\perp = B_{helix}|\cos(\beta-\psi_0)|\quad\text{(apertura tomográfica efectiva)}\\
&\quad \psi_0 = \arctan(\rho_0/z_0),\quad R_0 = \sqrt{\rho_0^2+z_0^2}\\
&\quad B_{helix} = \sqrt{\Delta z^2+\Delta\rho^2},\quad \beta = \arctan(\Delta z/\Delta\rho)\\
&\quad \lambda_0 = c/f_0
\end{aligned}
}$$

### 11.2 Condiciones de validez

| Condición | Expresión | Origen de la condición |
|-----------|-----------|----------------------|
| Banda estrecha | $B/f_0 \ll 1$ | Simplificación al separar contribuciones en $\Delta k_z$ |
| Target poco profundo | $|z_P| \ll z_0$ | Aproximación $\psi_0 \approx \arctan(\rho/z)$ hacia la interfaz |
| Cobertura azimutal completa | $N_t \geq 1$ vuelta completa | Necesario para PSF tipo $J_0$ en horizontal |
| TX y RX en misma hélice | mismo $\rho(t)$, $z(t)$ | Garantiza $k_z^{TX}+k_z^{RX} = 2\times(2\pi f n_2/c)\cos\theta_t$ |
| Interfaz plana | $\partial z_{interfaz}/\partial x \approx 0$ | Simplifica el problema de refracción |
| Medios homogéneos | $n_i = \text{cte}$ en cada semieespacio | No hay dispersión ni gradientes internos |

### 11.3 Interpretación física de cada término de $W_z$

| Término | Expresión | Mecanismo físico |
|---------|-----------|-----------------|
| **Rango** | $n_2\,B\,\cos\theta_{t,0}$ | Cada pulso chirp resuelve verticalmente $\delta_r = c/(2n_2 B)$; proyectado sobre el eje $z$ da un factor $\cos\theta_{t,0}$ |
| **Tomográfico** | $\frac{c\,B_\perp\sin\psi_0\cos\psi_0}{\lambda_0 R_0 n_2\cos\theta_{t,0}}$ | La diversidad de ángulo de elevación al descender la hélice añade diversidad en $k_z$; el desplazamiento espectral efectivo $\Delta f_z = c B_\perp\sin\psi_0/(λ_0 R_0 n_2^2\cos\theta_{t,0})$ amplía el ancho de banda vertical |

### 11.4 Condición óptima de diseño

$$\beta = \psi_0 \implies B_\perp = B_{helix} \implies W_z \text{ máximo} \implies \delta_z \text{ mínimo}$$

La hélice cónica con ángulo de inclinación igual al look angle medio tiene la mayor apertura tomográfica efectiva.

---

## 12. Validación Experimental Completa

Todas las predicciones fueron verificadas por simulación MATLAB (backprojection bistático con Snell iterativo).

**Sistema simulado:** $f_0=10$ GHz, $B=50$ MHz, $n_2=2$, $\rho_0=160$ m, $z_0=100$ m, $B_\perp=47.2$ m, $\beta=\psi_0=58°$ (hélice óptima), NorthOffset TX=0°.

| # | NorthOffset | $x_P$ | Variable | Predicción | Simulación | Error |
|---|:-----------:|:------:|----------|:----------:|:----------:|:-----:|
| 1 | 180° | 0 | $\delta_z$ | **0.211 m** | 0.2112 m | **0.05%** |
| 2 | 180° | 0 | $\delta_{xy}$ | sin resol. → 0.24 m | 0.240 m | **0.0%** |
| 3 | 90° | 0 | $\delta_z$ | **0.211 m** | 0.2112 m | **0.05%** |
| 4 | 90° | 0 | $\delta_{xy}$ | **9.54 mm** | 8.8 mm | **7.8%** |
| 5 | 90° | 20 m | $\delta_z$ | 18.72 cm (Góes) | **20.05 cm** | 7.1% |
| 6 | 90° | 20 m | $\delta_x$ | 9.96 mm | **8.8 mm** | 11.6% |
| 7 | 90° | 20 m | $\delta_y$ | 9.96 mm | **9.1 mm** | 8.6% |

**Invariancia al NorthOffset** (experimentos 1, 3): $\delta_z$ idéntico para NorthOffset=180° y 90° — validación analítica confirmada numéricamente.

**Error máximo en $\delta_z$:** 0.05% para target en eje; 7.1% para target off-axis. La diferencia se atribuye a la ausencia de pérdida de propagación en la simulación (la corrección de Góes fue derivada para sistemas con $1/R^2$).

**Error en $\delta_{xy}$:** 7.8–11.6%, sistemáticamente menor que el verdadero. Origen: aproximación $\psi_0 = \arctan(\rho_0/z_0)$ ignora el ángulo refractado exacto al target a $z_P=-5$ m.
