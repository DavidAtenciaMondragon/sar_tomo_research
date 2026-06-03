# Resolución Espacial del SAR Biestático Helicoidal con Dos Medios

**Versión:** 1.0  
**Fecha:** 2026-05-31  
**Referencias:** `models/fase_modelo.md`, `models/modelo_fisico_base.md`  
**Configuración:** Biestático · Hélice cónica · Aire + Suelo ($n_2 = \sqrt{\varepsilon_r}$) · Interfaz plana

---

## Índice

1. [Función de Dispersión de Punto (PSF) desde la integral SAR](#1-función-de-dispersión-de-punto-psf-desde-la-integral-sar)
2. [Representación en el Espacio de Número de Onda](#2-representación-en-el-espacio-de-número-de-onda)
3. [Criterio General de Resolución](#3-criterio-general-de-resolución)
4. [Cálculo del Gradiente de la Fase respecto a $(x,y,z)$](#4-cálculo-del-gradiente-de-la-fase-respecto-a-xyz)
5. [Cobertura del k-espacio: componentes horizontal y vertical](#5-cobertura-del-k-espacio-componentes-horizontal-y-vertical)
6. [Extensión espectral $\Delta k_z$: contribuciones y expresión exacta](#6-extensión-espectral-delta-k_z-contribuciones-y-expresión-exacta)
7. [Extensión espectral $\Delta k_{xy}$: cobertura horizontal](#7-extensión-espectral-delta-k_xy-cobertura-horizontal)
8. [Aproximaciones justificadas hacia forma cerrada](#8-aproximaciones-justificadas-hacia-forma-cerrada)
9. [Expresiones de resolución final](#9-expresiones-de-resolución-final)
10. [Casos especiales y verificación](#10-casos-especiales-y-verificación)
11. [Tabla resumen y mapa de dependencias](#11-tabla-resumen-y-mapa-de-dependencias)

---

## 1. Función de Dispersión de Punto (PSF) desde la Integral SAR

### 1.1 Imagen SAR como suma coherente

La imagen SAR formada por backprojection para un píxel de prueba $\mathbf{p}$ es (de `fase_modelo.md`, §9):

$$I(\mathbf{p}) = \sum_{n=1}^{N_p}\sum_{k=1}^{K} S(f_k, t_n)\,e^{+j\frac{2\pi f_k}{c}R_{OP}(\mathbf{r}_{TX,n},\mathbf{r}_{RX,n},\mathbf{p})}$$

Para un **dispersor puntual** en $\mathbf{P}$ con reflectividad compleja $g(\mathbf{P})$:

$$S(f_k, t_n) = g(\mathbf{P})\,A(f_k, t_n)\,e^{-j\frac{2\pi f_k}{c}R_{OP}(\mathbf{r}_{TX,n},\mathbf{r}_{RX,n},\mathbf{P})}$$

Sustituyendo en la imagen:

$$I(\mathbf{p}) = g(\mathbf{P})\sum_{n,k} A(f_k,t_n)\,e^{+j\frac{2\pi f_k}{c}\left[R_{OP}(\ldots,\mathbf{p}) - R_{OP}(\ldots,\mathbf{P})\right]}$$

### 1.2 Definición de la PSF

La **Función de Dispersión de Punto (PSF)** describe la respuesta del sistema a un blanco puntual ideal. Sea $\boldsymbol{\delta} = \mathbf{p} - \mathbf{P}$ el desplazamiento desde el foco verdadero:

$$\text{PSF}(\boldsymbol{\delta}) = \sum_{n,k} A(f_k,t_n)\,e^{+j\frac{2\pi f_k}{c}\Delta R_{OP}(t_n, f_k;\,\boldsymbol{\delta})}$$

donde:

$$\Delta R_{OP}(t_n, f_k;\,\boldsymbol{\delta}) \equiv R_{OP}(\mathbf{r}_{TX,n},\mathbf{r}_{RX,n},\mathbf{P}+\boldsymbol{\delta}) - R_{OP}(\mathbf{r}_{TX,n},\mathbf{r}_{RX,n},\mathbf{P})$$

La **resolución espacial** del sistema queda determinada por la anchura de $|\text{PSF}(\boldsymbol{\delta})|$ en cada dirección.

---

## 2. Representación en el Espacio de Número de Onda

### 2.1 Linealización de $\Delta R_{OP}$ — Approximación de Campo Próximo Linealizado

Para $|\boldsymbol{\delta}|$ pequeño respecto a las distancias de propagación ($|\boldsymbol{\delta}| \ll d_i$), expandimos $R_{OP}$ en serie de Taylor a **primer orden** alrededor de $\mathbf{P}$:

$$R_{OP}(\ldots,\mathbf{P}+\boldsymbol{\delta}) = R_{OP}(\ldots,\mathbf{P}) + \nabla_\mathbf{P}R_{OP}\cdot\boldsymbol{\delta} + \mathcal{O}(|\boldsymbol{\delta}|^2)$$

De donde:

$$\Delta R_{OP} \approx \nabla_\mathbf{P}R_{OP}(t_n)\cdot\boldsymbol{\delta}$$

> **Justificación:** Esta linealización es válida cuando la variación cuadrática es despreciable frente a $\lambda/4$ (criterio de campo lejano estándar en SAR): $|\boldsymbol{\delta}|^2/(2d_i) \ll \lambda/4$, es decir $|\boldsymbol{\delta}| \ll \sqrt{d_i\lambda/2}$.

### 2.2 Sustitución en la PSF

Usando el resultado de `fase_modelo.md` (§3, Danskin):

$$\nabla_\mathbf{P}R_{OP}(t_n) = n_2\!\left(\hat{\mathbf{e}}_2^{TX}(t_n) + \hat{\mathbf{e}}_2^{RX}(t_n)\right)$$

la PSF en primer orden se convierte en:

$$\text{PSF}(\boldsymbol{\delta}) = \sum_{n,k} A(f_k,t_n)\,\exp\!\left(+j\,\frac{2\pi f_k\,n_2}{c}\!\left[\hat{\mathbf{e}}_2^{TX}(t_n)+\hat{\mathbf{e}}_2^{RX}(t_n)\right]\cdot\boldsymbol{\delta}\right)$$

### 2.3 Cambio de variable al k-espacio

Definiendo el **vector de número de onda instantáneo** (de `fase_modelo.md`, §5):

$$\mathbf{k}(f_k, t_n) = \frac{2\pi f_k\,n_2}{c}\!\left(\hat{\mathbf{e}}_2^{TX}(t_n) + \hat{\mathbf{e}}_2^{RX}(t_n)\right)$$

la PSF adopta la forma de una **transformada de Fourier**:

$$\boxed{\text{PSF}(\boldsymbol{\delta}) = \sum_{n,k} A(f_k,t_n)\,e^{+j\,\mathbf{k}(f_k,t_n)\cdot\boldsymbol{\delta}} \approx \int_{\mathcal{K}} W(\mathbf{k})\,e^{+j\,\mathbf{k}\cdot\boldsymbol{\delta}}\,d^3k}$$

donde $W(\mathbf{k})$ es la **densidad de muestreo** en el k-espacio (Jacobiano de la transformación $(f_k, t_n) \to \mathbf{k}$, ponderada por $A$).

> **Conclusión fundamental:** La imagen SAR en el entorno del foco es la **transformada de Fourier inversa** de la función de densidad del k-espacio $W(\mathbf{k})$. La PSF es ancha donde el soporte espectral es estrecho y angosta donde el soporte es amplio.

---

## 3. Criterio General de Resolución

### 3.1 Principio de incertidumbre de Fourier

Para una distribución $W(\mathbf{k})$ con soporte en el intervalo $[k_{u,min},\, k_{u,max}]$ a lo largo de la dirección $\hat{u}$, la anchura de la PSF proyectada satisface:

$$\delta_u \cdot \Delta k_u \sim 2\pi$$

con el **ancho espectral efectivo** en dirección $\hat{u}$:

$$\Delta k_u \equiv \max_{f_k, t_n}\!\left(\mathbf{k}\cdot\hat{u}\right) - \min_{f_k, t_n}\!\left(\mathbf{k}\cdot\hat{u}\right)$$

### 3.2 Criterio de resolución operacional

$$\boxed{\delta_u = \frac{C_w \cdot 2\pi}{\Delta k_u}}$$

donde $C_w$ es un factor numérico que depende de la **forma de la ventana** de ponderación y de la **forma del soporte espectral**:

| Tipo de soporte | Tipo de PSF | $C_w$ | Criterio |
|-----------------|-------------|--------|----------|
| Rectangular 1D | $\mathrm{sinc}$ | $0.886$ | $-3\,\text{dB}$ (FWHM) |
| Hamming 1D | $\mathrm{sinc}$ apodizado | $1.30$ | $-3\,\text{dB}$ |
| Anular completo (circular SAR) | $J_0$ (Bessel orden 0) | $\approx 1.12\cdot(2\pi/4\pi) \to $ ver §7 | $e^{-1}$ radio |
| Gaussiano (Spiral tomografía) | Gaussiana | $\sqrt{\ln 2/\pi} \approx 0.470$ | HWHM |

> **Convención de este documento:** se usa el criterio de **primer nulo** ($C_w = 1$) como definición de resolución nominal, anotando explícitamente las conversiones cuando se compare con la literatura.

### 3.3 Resolución en cada dirección

Para las tres direcciones Cartesianas:

$$\delta_x = \frac{2\pi}{\Delta k_x}, \qquad \delta_y = \frac{2\pi}{\Delta k_y}, \qquad \delta_z = \frac{2\pi}{\Delta k_z}$$

La tarea central es calcular $\Delta k_x$, $\Delta k_y$, $\Delta k_z$ exactamente y luego aproximar.

---

## 4. Cálculo del Gradiente de la Fase respecto a $(x,y,z)$

### 4.1 Componentes del vector $\mathbf{k}$

El vector de número de onda (de `fase_modelo.md`, §5) con índices explícitos es:

$$\mathbf{k}(f_k, t_n) = \frac{2\pi f_k n_2}{c}\begin{pmatrix} (\hat{\mathbf{e}}_2^{TX})_x + (\hat{\mathbf{e}}_2^{RX})_x \\ (\hat{\mathbf{e}}_2^{TX})_y + (\hat{\mathbf{e}}_2^{RX})_y \\ (\hat{\mathbf{e}}_2^{TX})_z + (\hat{\mathbf{e}}_2^{RX})_z \end{pmatrix}$$

### 4.2 Vectores unitarios en el Medio 2 — expresión angular

Para un objetivo en $\mathbf{P} = (x_P, y_P, z_P)$ con $z_P < 0$, y puntos de refracción $\mathbf{Q}_{TX}^* = (Q_x^T, Q_y^T, 0)$ y $\mathbf{Q}_{RX}^* = (Q_x^R, Q_y^R, 0)$, los vectores unitarios son:

$$\hat{\mathbf{e}}_2^{TX} = \frac{\mathbf{P}-\mathbf{Q}_{TX}^*}{d_2^{TX}} = \frac{1}{d_2^{TX}}\begin{pmatrix}x_P - Q_x^T \\ y_P - Q_y^T \\ z_P\end{pmatrix}$$

Introduciendo los **ángulos de transmisión** $\theta_t^{TX}$, $\phi_{TX}$ (ángulos polar y azimutal del rayo en Medio 2 medidos desde $-\hat{z}$):

$$\hat{\mathbf{e}}_2^{TX} = \begin{pmatrix}\sin\theta_t^{TX}\cos\phi_{TX} \\ \sin\theta_t^{TX}\sin\phi_{TX} \\ -\cos\theta_t^{TX}\end{pmatrix}$$

donde:
- $\cos\theta_t^{TX} = |z_P|/d_2^{TX}$ (componente vertical hacia abajo, positiva)
- $\sin\theta_t^{TX} = \sqrt{(x_P-Q_x^T)^2+(y_P-Q_y^T)^2}\,/\,d_2^{TX}$
- $\phi_{TX}$: ángulo azimutal del rayo en Medio 2 (igual al ángulo del TX en el plano horizontal, para interfaz plana)

Análogamente para $\hat{\mathbf{e}}_2^{RX}$ con ángulos $\theta_t^{RX}$ y $\phi_{RX}$.

### 4.3 Componentes cartesianas del vector $\mathbf{k}$

Sustituyendo:

$$k_x(f_k,t_n) = \frac{2\pi f_k n_2}{c}\!\left(\sin\theta_t^{TX}\cos\phi_{TX} + \sin\theta_t^{RX}\cos\phi_{RX}\right) \tag{k$_x$}$$

$$k_y(f_k,t_n) = \frac{2\pi f_k n_2}{c}\!\left(\sin\theta_t^{TX}\sin\phi_{TX} + \sin\theta_t^{RX}\sin\phi_{RX}\right) \tag{k$_y$}$$

$$k_z(f_k,t_n) = -\frac{2\pi f_k n_2}{c}\!\left(\cos\theta_t^{TX} + \cos\theta_t^{RX}\right) \tag{k$_z$}$$

> **Signo de $k_z$:** negativo porque $\hat{\mathbf{e}}_2^{TX}\cdot\hat{z} = -\cos\theta_t^{TX} < 0$. El número de onda vertical apunta en la dirección $-\hat{z}$ (hacia el suelo), consistente con la propagación hacia el objetivo.

### 4.4 Dependencias explícitas en las variables del sistema

| Componente | Depende de | Varía con |
|------------|-----------|-----------|
| $k_x$, $k_y$ | $\phi_{TX}(t_n)$, $\phi_{RX}(t_n)$: azimut del sensor | Rotación de la hélice (período $T_{vuelta}$) |
| $k_x$, $k_y$ | $\theta_t^{TX}(t_n)$, $\theta_t^{RX}(t_n)$: ángulo polar | Cambio de $\rho$ y $z$ en la hélice cónica |
| $k_z$ | $\theta_t^{TX}(t_n)$, $\theta_t^{RX}(t_n)$ | Descenso helicoidal (cambio de $\psi$) |
| Todos | $f_k$ | Variación en la banda $[f_0-B/2,\, f_0+B/2]$ |

---

## 5. Cobertura del k-espacio: Componentes Horizontal y Vertical

### 5.1 Estructura del k-espacio para SAR helicoidal

Al variar $(f_k, t_n)$, el vector $\mathbf{k}$ traza una superficie en $\mathbb{R}^3$. Para el SAR helicoidal con cobertura azimutal completa ($\alpha \in [0, 2\pi N_t]$), la cobertura tiene la estructura siguiente:

**En el plano horizontal $(k_x, k_y)$:** Para ángulo fijo $f_k = f_0$ y ángulo de incidencia fijo $\theta_t^{TX}$, el vector $(\phi_{TX}, \phi_{RX})$ barre todos los azimuts → la proyección $(k_x, k_y)$ traza una **curva cerrada** (círculo para monoestático, elipse para biestático con separación fija).

**En la componente $k_z$:** No depende del azimut $\phi$, sino solo de los ángulos de transmisión polar $\theta_t$ y de la frecuencia $f_k$. Varía **lentamente** con el número de vuelta (descenso helicoidal).

### 5.2 Conexión con los ángulos en el Medio 1 (ley de Snell)

Los ángulos de transmisión $\theta_t$ en el Medio 2 están relacionados con los ángulos de incidencia $\psi$ en el Medio 1 (ángulo de incidencia desde la vertical = look angle) por:

$$\sin\theta_t = \frac{n_1}{n_2}\sin\psi = \frac{\sin\psi}{n_2}, \qquad \cos\theta_t = \sqrt{1 - \frac{\sin^2\psi}{n_2^2}}$$

El look angle $\psi$ de la hélice en función del tiempo (para el objetivo en $\mathbf{P} = (x_P, y_P, z_P)$):

$$\psi(t) \approx \arctan\!\left(\frac{\sqrt{(\rho(t)-\rho_P)^2 + \rho_P^2 - 2\rho(t)\rho_P\cos(\alpha(t)-\alpha_P)}}{z(t) - z_P}\right)$$

donde $\rho_P = \sqrt{x_P^2+y_P^2}$ y $\alpha_P = \arctan(y_P/x_P)$.

Para el objetivo **en el eje de la hélice** ($x_P = y_P = 0$, caso de referencia):

$$\psi(t) = \arctan\!\left(\frac{\rho(t)}{z(t) - z_P}\right) \approx \arctan\!\left(\frac{\rho(t)}{z(t)}\right) \quad \text{si } |z_P| \ll z(t)$$

---

## 6. Extensión Espectral $\Delta k_z$: Contribuciones y Expresión Exacta

### 6.1 Expresión exacta de $k_z$

$$k_z(f_k, t_n) = -\frac{2\pi f_k n_2}{c}\!\left(\cos\theta_t^{TX}(t_n) + \cos\theta_t^{RX}(t_n)\right)$$

con $\cos\theta_t = \sqrt{1 - \sin^2\psi/n_2^2}$ y $\psi = \psi(t_n)$.

### 6.2 Dos fuentes independientes de variación en $k_z$

**Fuente 1 — Variación de frecuencia** (a geometría fija $t_n = t_0$):

$$\frac{\partial k_z}{\partial f_k}\bigg|_{t_0} = -\frac{2\pi n_2}{c}\!\left(\cos\theta_t^{TX}(t_0) + \cos\theta_t^{RX}(t_0)\right)$$

Contribución al rango de $k_z$ de la variación en $f$:

$$\Delta k_z^{(f)}(t_0) = \left|\frac{\partial k_z}{\partial f_k}\right| B = \frac{2\pi n_2 B}{c}\!\left(\cos\theta_t^{TX}(t_0) + \cos\theta_t^{RX}(t_0)\right) \tag{A}$$

> Esto es la **contribución de rango** a la resolución vertical. Es equivalente a la resolución en rango proyectada sobre el eje $z$.

**Fuente 2 — Variación de geometría** (a frecuencia fija $f_k = f_0$, la hélice desciende):

$$\frac{\partial k_z}{\partial \psi}\bigg|_{f_0} = \frac{2\pi f_0 n_2}{c}\frac{\partial}{\partial\psi}\!\left(\cos\theta_t^{TX}+\cos\theta_t^{RX}\right)$$

Diferenciando $\cos\theta_t = \sqrt{1 - \sin^2\psi/n_2^2}$ respecto a $\psi$:

$$\frac{d\cos\theta_t}{d\psi} = -\frac{\sin\psi\cos\psi}{n_2^2\cos\theta_t}$$

entonces:

$$\frac{\partial k_z}{\partial \psi}\bigg|_{f_0} = \frac{2\pi f_0 n_2}{c}\cdot\frac{\sin\psi\cos\psi}{n_2^2\cos\theta_t}\cdot 2 = \frac{4\pi f_0 \sin\psi\cos\psi}{c\,n_2\cos\theta_t}$$

(monostático: $\theta_t^{TX} = \theta_t^{RX} = \theta_t$)

La variación del ángulo $\psi$ entre la cima y la base de la hélice es $\Delta\psi = \psi_{base} - \psi_{top}$. La contribución al rango de $k_z$:

$$\Delta k_z^{(geom)} = \frac{4\pi f_0 \sin\psi_0\cos\psi_0}{c\,n_2\cos\theta_{t,0}}\cdot\Delta\psi \tag{B}$$

### 6.3 Expresión exacta del rango total $\Delta k_z$

El rango completo de $k_z$ en toda la apertura $(f_k, t_n)$ es:

$$k_z^{max} = \frac{2\pi f_{max} n_2}{c}\!\left(\cos\theta_t^{TX,min} + \cos\theta_t^{RX,min}\right)$$

$$k_z^{min} = \frac{2\pi f_{min} n_2}{c}\!\left(\cos\theta_t^{TX,max} + \cos\theta_t^{RX,max}\right)$$

donde $\theta_t^{min/max}$ se alcanzan en los extremos de la hélice.

$$\boxed{\Delta k_z = \frac{2\pi n_2}{c}\!\left[f_{max}\!\left(\cos\theta_t^{TX,min}+\cos\theta_t^{RX,min}\right) - f_{min}\!\left(\cos\theta_t^{TX,max}+\cos\theta_t^{RX,max}\right)\right]} \tag{exacta}$$

> Esta expresión no tiene forma cerrada en términos de parámetros del sistema sin especificar la geometría de la hélice, porque los ángulos $\theta_t$ dependen implícitamente de $\rho(t)$ y $z(t)$ a través de Snell.

---

## 7. Extensión Espectral $\Delta k_{xy}$: Cobertura Horizontal

### 7.1 Estructura del k-espacio horizontal

Para un objetivo en el eje de la hélice ($x_P = y_P = 0$) y cobertura azimutal completa ($\alpha_{TX} \in [0, 2\pi N_t]$):

**Caso monoestático** ($\mathbf{r}_{TX} = \mathbf{r}_{RX}$, $\phi_{TX} = \phi_{RX} = \alpha(t)$):

$$(k_x, k_y) = \frac{4\pi f_k n_2}{c}\sin\theta_t(t)\cdot(-\cos\alpha(t),\; -\sin\alpha(t))$$

A frecuencia $f_k$ y ángulo $\theta_t$ fijos, al variar $\alpha \in [0, 2\pi]$, se traza un **círculo** de radio:

$$R_k(f_k, t) = \frac{4\pi f_k n_2}{c}\sin\theta_t(t) \tag{radio del círculo en k-espacio}$$

Al variar $f_k \in [f_0-B/2, f_0+B/2]$ y $\theta_t \in [\theta_t^{min}, \theta_t^{max}]$, la cobertura $(k_x, k_y)$ es un **anillo** (corona circular).

**Radio mínimo del anillo:**

$$R_{in} = \frac{4\pi f_{min}\,n_2}{c}\sin\theta_t^{min}$$

**Radio máximo del anillo:**

$$R_{out} = \frac{4\pi f_{max}\,n_2}{c}\sin\theta_t^{max}$$

### 7.2 Anchura del anillo en k-espacio horizontal

$$\Delta k_{ring} = R_{out} - R_{in} = \frac{4\pi n_2}{c}\!\left(f_{max}\sin\theta_t^{max} - f_{min}\sin\theta_t^{min}\right)$$

### 7.3 Cobertura total $\Delta k_x = \Delta k_y$

Para cobertura azimutal completa (360°):

$$\max k_x - \min k_x = 2R_{out} = \frac{8\pi f_{max}\,n_2}{c}\sin\theta_t^{max}$$

$$\boxed{\Delta k_{xy} \equiv \Delta k_x = \Delta k_y = 2R_{out} = \frac{8\pi f_{max}\,n_2}{c}\sin\theta_t^{max}} \tag{exacta}$$

> Esta es la extensión del **diámetro** del anillo de k-espacio. Sin embargo, la PSF en el plano horizontal no es una función sinc sino una función de Bessel $J_0$, lo que cambia el factor numérico de resolución (ver §8.3).

---

## 8. Aproximaciones Justificadas hacia Forma Cerrada

A continuación se introducen aproximaciones en orden creciente de restricción, con justificación física en cada paso.

### 8.1 Aproximación 1 — Objetivo en el eje de la hélice

**Enunciado:** $x_P = y_P = 0$ (objetivo exactamente bajo el centro de la hélice).

**Justificación:** El caso más simple y la geometría de referencia para SAR circular/espiral. Los ángulos $\psi$ y $\theta_t$ dependen solo de $\rho(t)/|z(t)|$, no del azimut.

**Consecuencias:**
- Los puntos de refracción $\mathbf{Q}_{TX}^*$ y $\mathbf{Q}_{RX}^*$ se sitúan exactamente bajo TX y RX respectivamente (azimut idéntico al sensor).
- Los ángulos $\phi_{TX} = \alpha_{TX}+\pi$ y $\phi_{RX} = \alpha_{RX}+\pi$ (en dirección opuesta al sensor).
- Los ángulos $\theta_t$ dependen solo del tiempo $t$ (a través de $\rho(t)$ y $z(t)$), no del azimut.

### 8.2 Aproximación 2 — Banda estrecha ($B \ll f_0$)

**Enunciado:** La variación fraccional de frecuencia es pequeña: $B/f_0 \ll 1$.

**Justificación:** Típica en sistemas SAR microondas. Para P-band ($f_0 = 435$ MHz, $B = 50$–85 MHz): $B/f_0 \approx 0.11$–0.20. Para sistemas de menor BW fraccional, la aproximación es excelente.

**Consecuencias:**
- Se reemplaza $f_{max} \approx f_{min} \approx f_0$ donde no sea la **diferencia** la que importa.
- La contribución de frecuencia a $\Delta k_z$:

$$\Delta k_z^{(f)} \approx \frac{4\pi n_2 B}{c}\cos\theta_{t,0} \tag{A'}$$

- El radio medio del anillo: $R_c \approx \frac{4\pi f_0 n_2}{c}\sin\theta_{t,0}$

### 8.3 Aproximación 3 — Variación angular lenta ($\Delta\psi \ll \psi_0$)

**Enunciado:** La variación del ángulo de incidencia entre la cima y la base de la hélice es pequeña: $|\psi_{base} - \psi_{top}| = \Delta\psi \ll \psi_0$.

**Justificación:** Para hélices con $N_t \gg 1$ vueltas y relación de aspecto $\Delta z/\rho_0$ moderada, $\Delta\psi$ es pequeño comparado con $\psi_0$.

**Consecuencias:** Linealizar $\cos\theta_t$ respecto a $\psi$:

$$\Delta\cos\theta_t \approx -\frac{d\cos\theta_t}{d\psi}\bigg|_{\psi_0}\Delta\psi = \frac{\sin\psi_0\cos\psi_0}{n_2^2\cos\theta_{t,0}}\Delta\psi$$

La contribución geométrica (monostática) a $\Delta k_z$:

$$\Delta k_z^{(geom)} \approx \frac{4\pi f_0 n_2}{c}\cdot\frac{\sin\psi_0\cos\psi_0}{n_2^2\cos\theta_{t,0}}\cdot\Delta\psi = \frac{4\pi f_0\sin\psi_0\cos\psi_0}{c\,n_2\cos\theta_{t,0}}\cdot\Delta\psi \tag{B'}$$

### 8.4 Conexión de $\Delta\psi$ con la apertura tomográfica efectiva $B_\perp$

La variación angular del look angle $\Delta\psi$ está relacionada con la longitud de la **baseline perpendicular a la LOS** (apertura tomográfica efectiva $B_\perp$) por:

$$\Delta\psi \approx \frac{B_\perp}{R_0}$$

donde $R_0 = \sqrt{\rho_0^2 + z_0^2}$ es la distancia media radar-objetivo.

> **Derivación:** Un desplazamiento $\Delta\mathbf{r}_{ant} = B_\perp\hat{\mathbf{n}}_\perp$ (perpendicular a la LOS) produce un cambio en el ángulo de incidencia de $\Delta\psi = B_\perp/R_0$ en la aproximación de campo lejano.

Sustituyendo:

$$\Delta k_z^{(geom)} \approx \frac{4\pi f_0\sin\psi_0\cos\psi_0}{c\,n_2\cos\theta_{t,0}}\cdot\frac{B_\perp}{R_0} = \frac{4\pi B_\perp}{\lambda_0 R_0}\cdot\frac{\sin\psi_0\cos\psi_0}{n_2\cos\theta_{t,0}} \tag{B''}$$

### 8.5 Aproximación 4 — Incidencia moderada con dos medios ($n_2 > 1$)

Para ángulos de incidencia $\psi_0$ tales que $\sin\psi_0 \ll n_2$ (rayo casi vertical en el suelo), $\cos\theta_{t,0} \approx 1$ y:

$$\Delta k_z^{(geom)} \approx \frac{4\pi B_\perp}{\lambda_0 R_0}\cdot\sin\psi_0\cos\psi_0 = \frac{2\pi B_\perp}{\lambda_0 R_0}\sin(2\psi_0) \tag{B'''}$$

Para el caso **un solo medio** ($n_2 = 1$, aire), usando $\cos\theta_t = \cos\psi$:

$$\frac{\sin\psi_0\cos\psi_0}{n_2\cos\theta_{t,0}}\bigg|_{n_2=1} = \frac{\sin\psi_0\cos\psi_0}{1\cdot\cos\psi_0} = \sin\psi_0$$

$$\Delta k_z^{(geom)}\bigg|_{n_2=1} = \frac{4\pi B_\perp}{\lambda_0 R_0}\sin\psi_0$$

Que es exactamente la **Ec. (4.24) de Góes 2022**: $\Delta k_z(B_\perp) \approx \frac{4\pi B_\perp}{\lambda_0 R_0}\sin\psi$. ✓

---

## 9. Expresiones de Resolución Final

### 9.1 Resolución Vertical $\delta_z$ — Dos medios, forma semi-cerrada

Combinando las contribuciones (A') y (B'') bajo las aproximaciones 1–3:

$$\Delta k_z \approx \frac{4\pi n_2 B}{c}\cos\theta_{t,0} + \frac{4\pi B_\perp}{\lambda_0 R_0}\cdot\frac{\sin\psi_0\cos\psi_0}{n_2\cos\theta_{t,0}}$$

$$= \frac{4\pi}{c}\!\left[n_2 B\cos\theta_{t,0} + \frac{cB_\perp}{\lambda_0 R_0}\cdot\frac{\sin\psi_0\cos\psi_0}{n_2\cos\theta_{t,0}}\right]$$

Definiendo el **ancho de banda efectivo vertical en dos medios**:

$$\boxed{W_z^{(2\,medios)} = \underbrace{n_2 B\cos\theta_{t,0}}_{\text{contribución de rango}} + \underbrace{\frac{cB_\perp}{\lambda_0 R_0}\cdot\frac{\sin\psi_0\cos\psi_0}{n_2\cos\theta_{t,0}}}_{\text{contribución tomográfica}}}$$

La resolución vertical (criterio primer nulo):

$$\boxed{\delta_z^{(2\,medios)} = \frac{c}{2\,W_z^{(2\,medios)}}}$$

> **Nota sobre el factor 2:** El factor $2$ en el denominador proviene de la relación $\Delta k_z = 4\pi W_z^{(2\,medios)}/c$, que al aplicar $\delta_z = 2\pi/\Delta k_z$ da $\delta_z = c/(2W_z)$. ✓

**Verificación dimensional:**
- $n_2 B\cos\theta_{t,0}$: $[\text{adimensional}][\text{Hz}][\text{adimensional}] = [\text{Hz}]$ ✓
- $cB_\perp/(\lambda_0 R_0)$: $[\text{m/s}][\text{m}]/([\text{m}][\text{m}]) = [\text{Hz}]$ ✓
- $\delta_z = c/(2W_z)$: $[\text{m/s}]/[\text{Hz}] = [\text{m}]$ ✓

### 9.2 Caso Límite Medio Único ($n_2 = 1$, sin refracción)

Con $n_2 = 1$: $\theta_{t,0} = \psi_0$ y $\cos\theta_{t,0} = \cos\psi_0$:

$$W_z^{(n_2=1)} = B\cos\psi_0 + \frac{cB_\perp}{\lambda_0 R_0}\sin\psi_0$$

Comparando con Góes 2022 (Ec. 4.25–4.27): $W_z = W\cos\psi_0 + \Delta f_z(B_\perp)$ con $\Delta f_z = \frac{cB_\perp}{\lambda_0 R_0}\sin\psi_0$. ✓

### 9.3 Resolución Horizontal $\delta_x = \delta_y$ — Dos medios

Para la PSF de tipo Bessel (cobertura anular completa en azimut), la resolución horizontal se define como el **radio del primer nulo** de $J_0(R_c\,r)$, que ocurre en $R_c\,r = 2.405$:

$$\delta_{xy}^{(nulo)} = \frac{2.405}{R_c}$$

con radio medio del anillo:

$$R_c = \frac{4\pi f_0 n_2}{c}\sin\theta_{t,0} = \frac{4\pi n_2}{\lambda_0}\sin\theta_{t,0}$$

luego:

$$\boxed{\delta_{xy}^{(nulo)} = \frac{2.405\,\lambda_0}{4\pi\,n_2\sin\theta_{t,0}}}$$

En términos del ángulo de incidencia en el aire $\psi_0$ (usando $\sin\theta_{t,0} = \sin\psi_0/n_2$):

$$\boxed{\delta_{xy}^{(nulo)} = \frac{2.405\,\lambda_0}{4\pi\sin\psi_0}}$$

> **Observación importante:** La resolución horizontal **no depende de $n_2$** cuando se expresa en términos del look angle en aire $\psi_0$. El índice de refracción $n_2$ afecta el ángulo de transmisión $\theta_{t,0}$, pero estos efectos se cancelan al usar el look angle $\psi_0$ medido en aire.

**Para criterio $e^{-1}$** de $J_0$ (Ishimaru 1998, $J_0(u)=e^{-1}$ en $u \approx 1.83$):

$$\delta_{xy}^{(e^{-1})} = \frac{1.83}{R_c} = \frac{1.83\,\lambda_0}{4\pi\sin\psi_0}$$

**Para criterio −3 dB** ($J_0(u)=1/\sqrt{2}$ en $u \approx 1.20$):

$$\delta_{xy}^{(-3\text{dB})} = \frac{1.20\,\lambda_0}{4\pi\sin\psi_0}$$

> La fórmula de Góes/Ishimaru $\delta_{xy} = 1.12\lambda/(2\pi\sin\psi)$ usa el punto donde $J_0(u) = e^{-1}$ (con $u \approx 1.12$ en la literatura) y el factor $2\pi$ en el denominador. Esta discrepancia con $4\pi$ sugiere que la literatura usa el número de onda **one-way** ($k = 2\pi/\lambda$) mientras que el k-space de la fase SAR usa el **two-way** ($k_{SAR} = 4\pi/\lambda$). Esto se clarifica en §10.4.

### 9.4 Resolución en Rango Slant ($\delta_r$, referencia)

La resolución en rango slant (en el Medio 2, proyectada a lo largo del rayo) es:

$$\delta_r = \frac{v_2}{2B} = \frac{c}{2n_2 B}$$

> Nótese que $\delta_r$ disminuye (mejora) al aumentar $n_2$: el suelo tiene longitudes de onda físicas más cortas, mejorando la resolución por unidad de ancho de banda.

---

## 10. Casos Especiales y Verificación

### 10.1 Medio único, monoestático ($n_2 = 1$, TX = RX, objetivo en superficie)

Con $n_2 = 1$, $z_P = 0$:

$$W_z^{(n_2=1)} = B\cos\psi_0 + \frac{cB_\perp}{\lambda_0 R_0}\sin\psi_0, \qquad \delta_z = \frac{c}{2W_z}$$

$$\delta_{xy} = \frac{2.405\,\lambda_0}{4\pi\sin\psi_0}, \qquad \delta_r = \frac{c}{2B}$$

Coincide con los resultados de SAR circular/espiral en literatura (Góes 2022). ✓

### 10.2 SAR circular puro (una sola vuelta, $B_\perp = 0$)

Para una sola vuelta sin cambio de altura: $B_\perp = 0$, por lo que:

$$W_z^{circ} = n_2 B\cos\theta_{t,0}, \qquad \delta_z^{circ} = \frac{c}{2n_2 B\cos\theta_{t,0}}$$

En el caso de un solo medio: $\delta_z^{circ} = c/(2B\cos\psi_0)$. Coincide con Góes Ec. (4.5) con $\psi$ como look angle. ✓

La resolución vertical de SAR circular es **mala** porque solo la frecuencia $B$ contribuye a $\Delta k_z$; no hay apertura tomográfica.

### 10.3 SAR helicoidal con $B \to 0$ (pulso monofrecuencia, solo geometría)

Con $B = 0$ (solo una frecuencia $f_0$), toda la resolución vertical proviene de la geometría:

$$\delta_z^{(B=0)} = \frac{c}{2} \cdot \frac{c/\lambda_0}{4\pi B_\perp/(R_0)} \cdot \frac{n_2\cos\theta_{t,0}}{\sin\psi_0\cos\psi_0}$$

$$= \frac{\lambda_0 R_0}{2B_\perp} \cdot \frac{n_2\cos\theta_{t,0}}{\sin\psi_0\cos\psi_0}$$

Para $n_2 = 1$: $\delta_z^{(B=0)} = \frac{\lambda_0 R_0}{2B_\perp\sin\psi_0}$, que es la fórmula de SAR Tomografía (Ec. 4.2 de Góes). ✓

### 10.4 Aclaración del factor 2π vs 4π (Ishimaru)

La fórmula de Ishimaru $\delta_{xy} = 1.12\lambda/(2\pi\sin\psi)$ usa la definición de k-espacio **one-way** ($k = 2\pi/\lambda$), donde la derivada de la fase respecto a $r$ da $k$, no $2k$.

En el modelo de este documento (two-way phase $\Phi = -(4\pi/\lambda)R$), el k-space tiene el factor doble. Reconciliación:

$$\delta_{xy}^{Ishimaru} = \frac{1.12\lambda}{2\pi\sin\psi} = \frac{1.12 \times 2.405^{-1} \times 2.405\lambda}{2\pi\sin\psi} = \frac{2.405\lambda}{2\pi \times 2.405/1.12 \cdot \sin\psi}$$

Usando $R_c^{one-way} = 2\pi\sin\psi/\lambda$ (one-way): $\delta = 1.12/R_c^{one-way}$.
Usando $R_c^{two-way} = 4\pi\sin\psi/\lambda$ (two-way): $\delta = 2.405/R_c^{two-way}$.

Ambas expresiones dan $\delta_{xy} = 2.405\lambda/(4\pi\sin\psi)$ con el criterio de primer nulo. La "diferencia" en los factores numéricos es solo una diferencia en el criterio de resolución:

$$\frac{1.12}{R_c^{one-way}} = \frac{1.12\lambda}{2\pi\sin\psi} = \frac{2.24\lambda}{4\pi\sin\psi}$$

No es exactamente igual a $2.405\lambda/(4\pi\sin\psi)$ (primer nulo), lo cual confirma que Ishimaru usa el punto $e^{-1}$ de $J_0$ (aproximadamente $u_{e^{-1}} \approx 1.83$, no 2.24). La literatura usa diferentes criterios; la tabla de §3 establece la correspondencia.

---

## 11. Tabla Resumen y Mapa de Dependencias

### 11.1 Expresiones de resolución (bajo aproximaciones 1–4)

| Dirección | Expresión | Criterio | Depende de |
|-----------|-----------|----------|-----------|
| **Rango** | $\delta_r = c/(2n_2 B)$ | Rayleigh sinc | $n_2$, $B$ |
| **Horizontal** $\delta_{xy}$ | $2.405\lambda_0/(4\pi\sin\psi_0)$ | Primer nulo $J_0$ | $\lambda_0$, $\psi_0$ |
| **Vertical** $\delta_z$ | $c/(2W_z)$ | Rayleigh sinc | $n_2$, $B$, $B_\perp$, $\psi_0$, $R_0$, $\lambda_0$ |

con:

$$W_z = n_2 B\cos\theta_{t,0} + \frac{cB_\perp}{\lambda_0 R_0}\cdot\frac{\sin\psi_0\cos\psi_0}{n_2\cos\theta_{t,0}}, \qquad \cos\theta_{t,0} = \sqrt{1 - \frac{\sin^2\psi_0}{n_2^2}}$$

### 11.2 Factores de mejora por presencia del suelo ($n_2 > 1$)

| Magnitud | Efecto de $n_2$ |
|---------|----------------|
| $\delta_r = c/(2n_2 B)$ | Mejora: $\delta_r$ disminuye con $n_2$ |
| $\delta_{xy} = 2.405\lambda_0/(4\pi\sin\psi_0)$ | **No cambia** con $n_2$ (en look angle $\psi_0$) |
| $\delta_z$ (término rango) $\sim c/(n_2 B\cos\theta_{t,0})$ | Mejora moderada: $n_2\cos\theta_{t,0}$ aumenta con $n_2$ |
| $\delta_z$ (término tomog.) $\sim \lambda_0 R_0 n_2\cos\theta_{t,0}/(B_\perp\sin\psi_0\cos\psi_0)$ | Empeora ligeramente con $n_2$ (ángulo refractado más pequeño → menor diversidad en $k_z$) |

### 11.3 Mapa de dependencias

```
Parámetros de diseño
│
├── Trayectoria helicoidal: {ρ_top, ρ_base, z_top, z_base, N_t, V₀}
│   └──→ ψ₀ = arctan(ρ₀/z₀),  B⊥ = f(β, ψ₀, B)
│
├── Radar: {f₀, B}
│   └──→ λ₀ = c/f₀
│
├── Geometría: R₀ = √(ρ₀²+z₀²)
│
└── Suelo: {ε_r → n₂ = √ε_r}
    └──→ θ_{t,0} = arcsin(sinψ₀/n₂)

         ↓ Resolución ↓

δ_r = c/(2n₂B)                     [rango]
δ_xy = 2.405λ₀/(4π sinψ₀)          [horizontal, independ. de n₂]
δ_z = c / [2(n₂B cosθ_{t,0} + cB⊥sinψ₀cosψ₀/(λ₀R₀n₂cosθ_{t,0}))]  [vertical]
```

---

## Limitaciones del Modelo en v1

Las siguientes hipótesis han sido utilizadas y deben revisarse para aplicaciones específicas:

1. **Linealización de $R_{OP}$ (primer orden):** válida para objetivos cuyo desplazamiento del foco es $|\boldsymbol{\delta}| \ll \sqrt{d_i\lambda/2}$. Para objetivos fuera del eje de la hélice se necesita el término cuadrático (de Hessiano).

2. **Objetivo en el eje de la hélice:** la resolución horizontal cambia para objetivos off-axis (el k-vector no traza círculos perfectos). Se requiere corrección del ángulo equivalente $\tilde{\psi}_0$ (ver Góes 2022 Ec. 4.52–4.53).

3. **Independencia de las contribuciones rango y tomográfica:** en la expresión $\Delta k_z = \Delta k_z^{(f)} + \Delta k_z^{(geom)}$, se ha asumido que el máximo total es la suma. Esto es válido cuando los rangos de $f$ y $\psi$ no se correlacionan en la cobertura del k-espacio (que en general no es exacto).

4. **Denominador del criterio de resolución:** se usa el criterio de primer nulo ($C_w = 1$). Para comparación con datos experimentales usar $C_w = 0.886$ (−3 dB, sinc) o la relación de Bessel apropiada.

5. **Factor numérico en $\delta_{xy}$:** el valor 2.405 (primer nulo de $J_0$) puede diferir del valor experimental debido a: (a) cobertura azimutal incompleta, (b) ponderación no uniforme, (c) target no en el eje.
