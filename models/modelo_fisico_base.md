# Modelo Físico Base — SAR Biestático Helicoidal con Dos Medios

**Versión:** 1.0  
**Fecha:** 2026-05-31  
**Referencia:** `notes/modelo_conceptual.md`  
**Configuración:** Biestático · Trayectoria helicoidal · Interfaz plana · Dos medios (aire + suelo) · Refracción (Snell)

---

## Índice

1. [Sistema de coordenadas y convenciones](#1-sistema-de-coordenadas-y-convenciones)
2. [Definición de los dos medios homogéneos](#2-definición-de-los-dos-medios-homogéneos)
3. [Parametrización de la trayectoria helicoidal](#3-parametrización-de-la-trayectoria-helicoidal)
4. [Definición del punto objetivo en el medio 2](#4-definición-del-punto-objetivo-en-el-medio-2)
5. [Geometría de propagación biestática en dos medios](#5-geometría-de-propagación-biestática-en-dos-medios)
6. [Puntos de refracción en la interfaz plana](#6-puntos-de-refracción-en-la-interfaz-plana)
7. [Ley de Snell en la interfaz](#7-ley-de-snell-en-la-interfaz)
8. [Distancia óptica total — expresión exacta](#8-distancia-óptica-total--expresión-exacta)
9. [Modelo de la señal SAR biestática](#9-modelo-de-la-señal-sar-biestática)
10. [Resumen de la cadena causal del modelo](#10-resumen-de-la-cadena-causal-del-modelo)

---

## 1. Sistema de Coordenadas y Convenciones

Se adopta un sistema de coordenadas **Cartesiano 3D** orientado de la siguiente manera:

$$\mathcal{F} = \{O;\; \hat{x},\; \hat{y},\; \hat{z}\}$$

- El **origen** $O$ es el punto de referencia en la interfaz entre los dos medios (nivel del suelo).
- El eje $\hat{z}$ apunta **verticalmente hacia arriba**: $z > 0$ en el aire, $z < 0$ en el subsuelo.
- El plano $z = 0$ coincide con la **interfaz plana** entre los dos medios.
- Los ejes $\hat{x}$ e $\hat{y}$ son horizontales y definen el plano de la apertura de adquisición.

**Convenio de vectores:**

Un punto genérico en $\mathbb{R}^3$ se escribe como:

$$\mathbf{r} = x\,\hat{x} + y\,\hat{y} + z\,\hat{z} \equiv (x, y, z)^T$$

La norma euclidiana de un vector $\mathbf{a}$ se denota $|\mathbf{a}| \equiv \|\mathbf{a}\|_2$.

---

## 2. Definición de los Dos Medios Homogéneos

### Medio 1 — Aire (semieespacio superior, $z > 0$)

| Magnitud | Símbolo | Valor |
|----------|---------|-------|
| Permitividad relativa | $\varepsilon_{r,1}$ | $1$ |
| Permeabilidad relativa | $\mu_{r,1}$ | $1$ |
| Índice de refracción | $n_1 = \sqrt{\varepsilon_{r,1}\,\mu_{r,1}}$ | $1$ |
| Velocidad de fase | $v_1 = c / n_1$ | $c$ |
| Número de onda | $k_1(f) = 2\pi f\,n_1 / c$ | $2\pi f/c$ |

### Medio 2 — Suelo (semieespacio inferior, $z < 0$)

| Magnitud | Símbolo | Expresión general | Valor típico |
|----------|---------|-------------------|-------------|
| Permitividad relativa | $\varepsilon_{r,2}$ | libre | $\varepsilon_r$ |
| Permeabilidad relativa | $\mu_{r,2}$ | $\approx 1$ (suelo no magnético) | $1$ |
| Índice de refracción | $n_2 = \sqrt{\varepsilon_{r,2}}$ | $\sqrt{\varepsilon_r}$ | $2$ (para $\varepsilon_r = 4$) |
| Velocidad de fase | $v_2 = c / n_2$ | $c/\sqrt{\varepsilon_r}$ | $c/2$ |
| Número de onda | $k_2(f) = 2\pi f\,n_2 / c$ | $2\pi f\sqrt{\varepsilon_r}/c$ | $4\pi f/c$ |

> **Nota sobre el valor típico:** en muchos suelos (arena húmeda, tierra arcillosa seca), $\varepsilon_r \approx 4$, lo que da $n_2 = 2$ y $v_2 = c/2$.

### Hipótesis de los medios

1. **Homogeneidad:** $n_1$ y $n_2$ son constantes (no dependen de la posición dentro de su semieespacio ni de la frecuencia en el rango de trabajo).
2. **Isotrópía:** la respuesta EM es igual en todas las direcciones.
3. **No dispersividad:** $n_i$ no depende de $f$ (equivalente a decir que la velocidad de fase = velocidad de grupo dentro del rango de trabajo).
4. **Sin pérdidas de propagación** (solo se modela el cambio de velocidad y la refracción, no la atenuación $\propto e^{-\alpha z}$).

---

## 3. Parametrización de la Trayectoria Helicoidal

### 3.1 Geometría General de la Hélice Cónica

Una hélice cónica en 3D se define por los siguientes parámetros geométricos:

| Símbolo | Descripción |
|---------|-------------|
| $\rho_{top}$ | Radio en la cima de la espiral (mayor altura) |
| $\rho_{base}$ | Radio en la base de la espiral (menor altura) |
| $z_{top}$ | Altura de la cima ($z_{top} > 0$) |
| $z_{base}$ | Altura de la base ($z_{base} > 0$, $z_{base} < z_{top}$) |
| $N_t$ | Número total de vueltas |
| $V_0$ | Velocidad tangencial del sensor (constante) |

La **velocidad radial** (tasa de cambio del radio) y la **velocidad vertical** son:

$$V_\rho = \frac{\rho_{base} - \rho_{top}}{t_{max}}, \qquad V_z = -\frac{z_{top} - z_{base}}{t_{max}} \quad (V_z < 0)$$

donde $t_{max}$ es el tiempo total de vuelo de la espiral completa (ver Ec. 3.6).

### 3.2 Radio en Función del Tiempo

El radio varía linealmente desde la cima hasta la base:

$$\boxed{\rho(t) = \rho_{top} + V_\rho\,t, \qquad 0 \leq t \leq t_{max}}$$

Para la **espiral cilíndrica** (caso particular, $\rho_{base} = \rho_{top} = \rho_0$): $V_\rho = 0$ y $\rho(t) = \rho_0$.

### 3.3 Ángulo de Azimut en Función del Tiempo

Imponiendo velocidad tangencial constante $V_0$, la velocidad angular es $\dot{\alpha}(t) = V_0/\rho(t)$, por lo que:

$$\boxed{\alpha(t) = \int_0^t \frac{V_0}{\rho(\tau)}\,d\tau = \frac{V_0}{V_\rho}\Big[\ln\!\big(\rho(t)\big) - \ln(\rho_{top})\Big], \qquad V_\rho \neq 0}$$

Para la espiral **cilíndrica** ($V_\rho = 0$, límite $V_\rho \to 0$):

$$\alpha(t) = \frac{V_0}{\rho_0}\,t$$

La condición $\alpha(t_{max}) = 2\pi N_t$ define:

$$\boxed{t_{max} = \frac{2\pi N_t}{V_0} \cdot \frac{\rho_{base} - \rho_{top}}{\ln(\rho_{base}) - \ln(\rho_{top})}} \qquad (V_\rho \neq 0)$$

$$t_{max} = \frac{2\pi N_t\,\rho_0}{V_0} \qquad \text{(cilíndrica)}$$

### 3.4 Altura en Función del Tiempo

$$\boxed{z(t) = z_{top} + V_z\,t = z_{top} - \frac{z_{top} - z_{base}}{t_{max}}\,t}$$

### 3.5 Posición Cartesiana de la Hélice

Las coordenadas Cartesianas en función del tiempo son:

$$\boxed{
\begin{aligned}
x(t) &= \rho(t)\cos\!\big(\alpha(t)\big) \\
y(t) &= \rho(t)\sin\!\big(\alpha(t)\big) \\
z(t) &= z_{top} + V_z\,t
\end{aligned}
}$$

### 3.6 Parametrización de TX y RX (Sistema Biestático)

En la configuración biestática, el **transmisor (TX)** y el **receptor (RX)** vuelan en trayectorias helicoidales **independientes**, posiblemente con parámetros distintos. Se definen:

**Transmisor TX:**

$$\mathbf{r}_{TX}(t_{TX}) = \begin{pmatrix} \rho_{TX}(t_{TX})\cos\!\big(\alpha_{TX}(t_{TX})\big) \\ \rho_{TX}(t_{TX})\sin\!\big(\alpha_{TX}(t_{TX})\big) \\ z_{TX}(t_{TX}) \end{pmatrix}$$

con parámetros $\{\rho_{top}^{TX},\, \rho_{base}^{TX},\, z_{top}^{TX},\, z_{base}^{TX},\, N_t^{TX},\, V_0^{TX}\}$.

**Receptor RX:**

$$\mathbf{r}_{RX}(t_{RX}) = \begin{pmatrix} \rho_{RX}(t_{RX})\cos\!\big(\alpha_{RX}(t_{RX})\big) \\ \rho_{RX}(t_{RX})\sin\!\big(\alpha_{RX}(t_{RX})\big) \\ z_{RX}(t_{RX}) \end{pmatrix}$$

con parámetros $\{\rho_{top}^{RX},\, \rho_{base}^{RX},\, z_{top}^{RX},\, z_{base}^{RX},\, N_t^{RX},\, V_0^{RX}\}$.

> **Condición de existencia en el Medio 1:** se exige $z_{TX}(t) > 0$ y $z_{RX}(t) > 0$ para todo $t$ (ambas antenas permanecen siempre en el aire).

> **Hipótesis stop-and-go:** la plataforma se asume estacionaria durante la transmisión y recepción de un pulso. Así, $t_{TX}$ y $t_{RX}$ son el **índice de slow-time** (pulso $n$) y se toman como parámetros fijos durante el cálculo de cada eco. Se escribe:

$$\mathbf{r}_{TX,n} \equiv \mathbf{r}_{TX}(t_n), \qquad \mathbf{r}_{RX,n} \equiv \mathbf{r}_{RX}(t_n)$$

---

## 4. Definición del Punto Objetivo en el Medio 2

El punto objetivo (reflector subsuperficial) está ubicado en el semieespacio del suelo:

$$\boxed{\mathbf{P} = (x_P,\; y_P,\; z_P)^T, \qquad z_P < 0}$$

**Parámetros:**
- $(x_P, y_P)$: posición horizontal del objetivo
- $z_P < 0$: profundidad bajo la interfaz (negativa en la convención adoptada)
- $|z_P|$: profundidad real en el medio 2

> **El objetivo es un dispersor puntual** isótropo (responde igual en todas las direcciones y frecuencias dentro del rango de trabajo). Sin esta suposición, la reflectividad sería una función $g(\mathbf{P}, f, \hat{\mathbf{r}}_{inc}, \hat{\mathbf{r}}_{obs})$ del ángulo de incidencia, observación y frecuencia.

---

## 5. Geometría de Propagación Biestática en Dos Medios

### 5.1 Descripción del Camino de Propagación

En la configuración biestática con dos medios, la señal recorre el siguiente camino (ver diagrama):

```
        r_TX(n)              r_RX(n)
          ·                    ·
           \                  /
   Medio 1  \  (aire, n₁=1)  /
  (z > 0)    \              /
              ·            ·
           Q_TX(n)      Q_RX(n)          ← Interfaz z = 0
   ─────────────────────────────────────────────────────
   Medio 2   \            /
  (z < 0)     \  (suelo, n₂)
               \          /
                \        /
                    P             ← Objetivo en z = z_P < 0
```

El camino completo de la señal es:

$$\mathbf{r}_{TX,n} \;\xrightarrow{\text{Medio 1}}\; \mathbf{Q}_{TX,n} \;\xrightarrow{\text{Medio 2}}\; \mathbf{P} \;\xrightarrow{\text{Medio 2}}\; \mathbf{Q}_{RX,n} \;\xrightarrow{\text{Medio 1}}\; \mathbf{r}_{RX,n}$$

### 5.2 Cuatro Segmentos de Propagación

| Segmento | Trayecto | Medio | Velocidad | Longitud geométrica |
|----------|----------|-------|-----------|---------------------|
| $\ell_1^{TX}$ | $\mathbf{r}_{TX,n} \to \mathbf{Q}_{TX,n}$ | 1 (aire) | $c$ | $d_1^{TX} = |\mathbf{r}_{TX,n} - \mathbf{Q}_{TX,n}|$ |
| $\ell_2^{TX}$ | $\mathbf{Q}_{TX,n} \to \mathbf{P}$ | 2 (suelo) | $v_2$ | $d_2^{TX} = |\mathbf{Q}_{TX,n} - \mathbf{P}|$ |
| $\ell_2^{RX}$ | $\mathbf{P} \to \mathbf{Q}_{RX,n}$ | 2 (suelo) | $v_2$ | $d_2^{RX} = |\mathbf{P} - \mathbf{Q}_{RX,n}|$ |
| $\ell_1^{RX}$ | $\mathbf{Q}_{RX,n} \to \mathbf{r}_{RX,n}$ | 1 (aire) | $c$ | $d_1^{RX} = |\mathbf{Q}_{RX,n} - \mathbf{r}_{RX,n}|$ |

### 5.3 Tiempo de Tránsito Total

El tiempo total que tarda la señal en recorrer todo el camino es:

$$\tau_n(\mathbf{P}, \mathbf{Q}_{TX,n}, \mathbf{Q}_{RX,n}) = \frac{d_1^{TX}}{c} + \frac{d_2^{TX}}{v_2} + \frac{d_2^{RX}}{v_2} + \frac{d_1^{RX}}{c}$$

Expresado en función del índice de refracción $n_2 = c/v_2$:

$$\boxed{\tau_n = \frac{1}{c}\!\left[\,d_1^{TX} + n_2\,d_2^{TX} + n_2\,d_2^{RX} + d_1^{RX}\,\right]}$$

---

## 6. Puntos de Refracción en la Interfaz Plana

### 6.1 Restricción Geométrica — Interfaz Plana

La interfaz es el plano $z = 0$. Por tanto, los puntos de refracción son de la forma:

$$\mathbf{Q}_{TX,n} = (Q_x^{TX},\; Q_y^{TX},\; 0)^T \in \mathbb{R}^2 \times \{0\}$$

$$\mathbf{Q}_{RX,n} = (Q_x^{RX},\; Q_y^{RX},\; 0)^T \in \mathbb{R}^2 \times \{0\}$$

Cada punto $\mathbf{Q}$ depende de la posición del sensor y de la posición del objetivo $\mathbf{P}$.

### 6.2 Principio de Fermat — Formulación Variacional

Los puntos de refracción son los que **minimizan el tiempo de tránsito** en sus respectivos sub-caminos. Esto es equivalente a la ley de Snell (ambos son consecuencia del principio de mínima acción óptica):

**Para el lado TX** (camino $\mathbf{r}_{TX,n} \to \mathbf{Q}_{TX,n} \to \mathbf{P}$):

$$\boxed{\mathbf{Q}_{TX,n}^* = \arg\min_{\mathbf{Q} \in \{z=0\}} \left[\frac{|\mathbf{r}_{TX,n} - \mathbf{Q}|}{c} + \frac{|\mathbf{Q} - \mathbf{P}|}{v_2}\right]}$$

**Para el lado RX** (camino $\mathbf{P} \to \mathbf{Q}_{RX,n} \to \mathbf{r}_{RX,n}$):

$$\boxed{\mathbf{Q}_{RX,n}^* = \arg\min_{\mathbf{Q} \in \{z=0\}} \left[\frac{|\mathbf{P} - \mathbf{Q}|}{v_2} + \frac{|\mathbf{Q} - \mathbf{r}_{RX,n}|}{c}\right]}$$

> **Observación clave:** los dos problemas de refracción (TX y RX) son **independientes entre sí**. El punto $\mathbf{Q}_{TX,n}^*$ depende únicamente de $(\mathbf{r}_{TX,n}, \mathbf{P})$, y el punto $\mathbf{Q}_{RX,n}^*$ depende únicamente de $(\mathbf{P}, \mathbf{r}_{RX,n})$.

### 6.3 Longitudes Explícitas de los Cuatro Segmentos

Dado que $\mathbf{r}_{TX,n} = (x_{TX}, y_{TX}, z_{TX})^T$ con $z_{TX} > 0$, $\mathbf{P} = (x_P, y_P, z_P)^T$ con $z_P < 0$, y $\mathbf{Q}_{TX} = (Q_x^{TX}, Q_y^{TX}, 0)^T$:

$$d_1^{TX} = \sqrt{(x_{TX} - Q_x^{TX})^2 + (y_{TX} - Q_y^{TX})^2 + z_{TX}^2}$$

$$d_2^{TX} = \sqrt{(Q_x^{TX} - x_P)^2 + (Q_y^{TX} - y_P)^2 + z_P^2}$$

$$d_2^{RX} = \sqrt{(x_P - Q_x^{RX})^2 + (y_P - Q_y^{RX})^2 + z_P^2}$$

$$d_1^{RX} = \sqrt{(Q_x^{RX} - x_{RX})^2 + (Q_y^{RX} - y_{RX})^2 + z_{RX}^2}$$

---

## 7. Ley de Snell en la Interfaz

### 7.1 Derivación a partir de Fermat

Para el lado TX, la condición de mínimo sobre $\mathbf{Q}_{TX} = (Q_x, Q_y, 0)$ es $\nabla_{\mathbf{Q}}\,\tau = \mathbf{0}$, lo que produce dos ecuaciones (una por componente horizontal):

**Componente $x$:**

$$\frac{\partial\tau}{\partial Q_x} = \frac{Q_x - x_{TX}}{c\,d_1^{TX}} + \frac{Q_x - x_P}{v_2\,d_2^{TX}} = 0$$

**Componente $y$:**

$$\frac{\partial\tau}{\partial Q_y} = \frac{Q_y - y_{TX}}{c\,d_1^{TX}} + \frac{Q_y - y_P}{v_2\,d_2^{TX}} = 0$$

Estas dos ecuaciones se pueden reescribir vectorialmente. Definamos:

$$\hat{\ell}_1^{TX} = \frac{\mathbf{Q}_{TX} - \mathbf{r}_{TX,n}}{|\mathbf{Q}_{TX} - \mathbf{r}_{TX,n}|} \quad \text{(dirección del rayo en Medio 1, lado TX)}$$

$$\hat{\ell}_2^{TX} = \frac{\mathbf{P} - \mathbf{Q}_{TX}}{|\mathbf{P} - \mathbf{Q}_{TX}|} \quad \text{(dirección del rayo en Medio 2, lado TX)}$$

Entonces las condiciones de mínimo equivalen a:

$$\frac{(\hat{\ell}_1^{TX})_x}{c} + \frac{(\hat{\ell}_2^{TX})_x}{v_2} = 0, \qquad \frac{(\hat{\ell}_1^{TX})_y}{c} + \frac{(\hat{\ell}_2^{TX})_y}{v_2} = 0$$

### 7.2 Ley de Snell en Forma Vectorial

Sea $\hat{\mathbf{n}} = \hat{z} = (0,0,1)^T$ la normal a la interfaz apuntando hacia el Medio 1. La **Ley de Snell vectorial** en la interfaz $z = 0$ para el rayo TX es:

$$\boxed{n_1\!\left(\hat{\ell}_1^{TX} \times \hat{\mathbf{n}}\right) = n_2\!\left(\hat{\ell}_2^{TX} \times \hat{\mathbf{n}}\right)}$$

O equivalentemente en términos de ángulos con la normal:

$$\boxed{n_1\sin\theta_i^{TX} = n_2\sin\theta_t^{TX}}$$

donde:
- $\theta_i^{TX}$: ángulo de **incidencia** del rayo TX con respecto a la normal $\hat{\mathbf{n}}$ (medido en el Medio 1)
- $\theta_t^{TX}$: ángulo de **transmisión/refracción** del rayo TX en el Medio 2

Con $n_1 = 1$:

$$\sin\theta_t^{TX} = \frac{\sin\theta_i^{TX}}{n_2} = \frac{\sin\theta_i^{TX}}{\sqrt{\varepsilon_r}}$$

Análogamente para el rayo RX:

$$\boxed{n_1\sin\theta_i^{RX} = n_2\sin\theta_t^{RX}} \quad \Longrightarrow \quad \sin\theta_t^{RX} = \frac{\sin\theta_i^{RX}}{\sqrt{\varepsilon_r}}$$

### 7.3 Ángulos en Términos de Coordenadas

Los ángulos de incidencia y transmisión se expresan exactamente como:

**Lado TX:**

$$\cos\theta_i^{TX} = \frac{z_{TX}}{d_1^{TX}} = \frac{z_{TX}}{\sqrt{(x_{TX}-Q_x^{TX})^2 + (y_{TX}-Q_y^{TX})^2 + z_{TX}^2}}$$

$$\cos\theta_t^{TX} = \frac{|z_P|}{d_2^{TX}} = \frac{|z_P|}{\sqrt{(Q_x^{TX}-x_P)^2 + (Q_y^{TX}-y_P)^2 + z_P^2}}$$

**Lado RX:**

$$\cos\theta_i^{RX} = \frac{z_{RX}}{d_1^{RX}} = \frac{z_{RX}}{\sqrt{(x_{RX}-Q_x^{RX})^2 + (y_{RX}-Q_y^{RX})^2 + z_{RX}^2}}$$

$$\cos\theta_t^{RX} = \frac{|z_P|}{d_2^{RX}} = \frac{|z_P|}{\sqrt{(Q_x^{RX}-x_P)^2 + (Q_y^{RX}-y_P)^2 + z_P^2}}$$

### 7.4 Sistema de Ecuaciones para $\mathbf{Q}_{TX}^*$ — Forma Implícita

Sustituyendo la ley de Snell en las condiciones de mínimo, se obtiene el sistema implícito para $(Q_x^{TX}, Q_y^{TX})$:

$$\frac{Q_x^{TX} - x_{TX}}{d_1^{TX}(Q_x^{TX}, Q_y^{TX})} = -\frac{n_2\,(Q_x^{TX} - x_P)}{d_2^{TX}(Q_x^{TX}, Q_y^{TX})}$$

$$\frac{Q_y^{TX} - y_{TX}}{d_1^{TX}(Q_x^{TX}, Q_y^{TX})} = -\frac{n_2\,(Q_y^{TX} - y_P)}{d_2^{TX}(Q_x^{TX}, Q_y^{TX})}$$

Este sistema de dos ecuaciones no lineales en dos incógnitas $(Q_x^{TX}, Q_y^{TX})$ **no tiene solución analítica cerrada en el caso 3D general**. Requiere resolución numérica (e.g., método de Newton-Raphson, bisección en 1D si el problema es planar, o algoritmos de minimización directa).

> **Caso especial 2D:** si $y_{TX} = y_P = y_{RX} = Q_y = 0$ (problema planar en el plano $xz$), el sistema se reduce a una sola ecuación trascendental en $Q_x$.

Análogamente para $(Q_x^{RX}, Q_y^{RX})$:

$$\frac{Q_x^{RX} - x_P}{d_2^{RX}} \cdot n_2 = -\frac{Q_x^{RX} - x_{RX}}{d_1^{RX}}$$

$$\frac{Q_y^{RX} - y_P}{d_2^{RX}} \cdot n_2 = -\frac{Q_y^{RX} - y_{RX}}{d_1^{RX}}$$

---

## 8. Distancia Óptica Total — Expresión Exacta

### 8.1 Definición del Camino Óptico

La **longitud de camino óptico** (Optical Path Length, OPL) es la distancia equivalente en el vacío que produciría el mismo desfase de fase que el camino real en los dos medios. Está definida como:

$$\text{OPL} = \sum_i n_i \cdot d_i$$

### 8.2 Expresión Exacta para el Sistema Biestático

Evaluando en los puntos de refracción óptimos $\mathbf{Q}_{TX,n}^*$ y $\mathbf{Q}_{RX,n}^*$:

$$\boxed{R_{OP}(\mathbf{r}_{TX,n}, \mathbf{r}_{RX,n}, \mathbf{P}) = d_1^{TX} + n_2\,d_2^{TX} + n_2\,d_2^{RX} + d_1^{RX}}$$

Expandido:

$$R_{OP} = \underbrace{|\mathbf{r}_{TX,n} - \mathbf{Q}_{TX,n}^*|}_{\text{TX en aire}} + n_2\,\underbrace{|\mathbf{Q}_{TX,n}^* - \mathbf{P}|}_{\text{TX en suelo}} + n_2\,\underbrace{|\mathbf{P} - \mathbf{Q}_{RX,n}^*|}_{\text{RX en suelo}} + \underbrace{|\mathbf{Q}_{RX,n}^* - \mathbf{r}_{RX,n}|}_{\text{RX en aire}}$$

### 8.3 Tiempo de Retraso del Eco

$$\boxed{\tau_n(\mathbf{P}) = \frac{R_{OP}(\mathbf{r}_{TX,n}, \mathbf{r}_{RX,n}, \mathbf{P})}{c}}$$

### 8.4 Caso Límite Monoestático ($\mathbf{r}_{TX} = \mathbf{r}_{RX}$ y $\mathbf{Q}_{TX}^* = \mathbf{Q}_{RX}^* = \mathbf{Q}^*$)

Cuando TX y RX coinciden (misma antena), y por simetría $\mathbf{Q}_{TX}^* = \mathbf{Q}_{RX}^* \equiv \mathbf{Q}^*$:

$$R_{OP}^{mono} = 2\,|\mathbf{r}_{ant} - \mathbf{Q}^*| + 2\,n_2\,|\mathbf{Q}^* - \mathbf{P}|$$

que coincide con la expresión del modelo BP refractive de Imoc 2023 y García-Fernández 2019:
$$R_{OP}^{mono} = 2(d_1 + n_2\,d_2)$$

### 8.5 Relación con el Factor Biestático $w_{tr}$

En el límite de campo lejano (lejos del eje de la espiral y con $R_0 \gg \Delta\rho$), el tiempo de eco biestático puede escribirse aproximadamente como:

$$\tau_n \approx \frac{|\mathbf{r}_{TX,n}| + |\mathbf{r}_{RX,n}|}{c} + \frac{w_{tr}}{c}\,(\hat{\mathbf{r}}_{tx} \cdot \Delta\mathbf{P}_\perp + n_2\,|z_P|)$$

donde $w_{tr} = 2\cos\beta$ es el factor biestático (semi-ángulo $\beta$) de Arikan & Munson 1988. Esta aproximación **no se usará aquí**; se mantiene la expresión exacta.

---

## 9. Modelo de la Señal SAR Biestática

### 9.1 Pulso Transmitido

El TX emite en el pulso $n$ un pulso chirp LFM complejo:

$$s_{TX}(\tau;\,n) = \mathrm{rect}\!\left(\frac{\tau}{T}\right)\exp\!\left[j\!\left(2\pi f_0\,\tau + \pi\,\gamma\,\tau^2\right)\right]$$

**Variables:**
- $\tau$ — fast-time (tiempo dentro del pulso)
- $f_0$ — frecuencia portadora central [Hz]
- $T$ — duración del pulso [s]
- $\gamma$ — tasa de chirp [Hz/s]; ancho de banda $B = \gamma T$

### 9.2 Señal Recibida (Eco de un Dispersor Puntual en P)

La señal recibida en el receptor RX, proveniente del dispersor puntual en $\mathbf{P}$ con reflectividad $g(\mathbf{P})$, es una versión retrasada y atenuada del pulso transmitido:

$$s_{RX}(\tau, t_n;\,\mathbf{P}) = \frac{g(\mathbf{P})\,A(\mathbf{r}_{TX,n}, \mathbf{P}, \mathbf{r}_{RX,n})}{R_{OP}^2} \cdot s_{TX}\!\left(\tau - \tau_n(\mathbf{P});\,n\right)$$

donde:
- $A(\cdot)$ — factor de amplitud (incluye patrón de antena, pérdidas de reflexión en la interfaz)
- $R_{OP}^2$ — pérdida de propagación de dos vías (ley de potencias)
- $\tau_n(\mathbf{P}) = R_{OP}/c$ — retraso de eco exacto

### 9.3 Demodulación (Baseband IQ)

Multiplicando por la portadora de referencia $\exp(-j2\pi f_0 \tau)$ y filtrando paso-bajo:

$$s_{bb}(\tau, t_n;\,\mathbf{P}) = \frac{g(\mathbf{P})\,A}{R_{OP}^2}\,\mathrm{rect}\!\left(\frac{\tau - \tau_n}{T}\right)\exp\!\left[j\pi\gamma(\tau-\tau_n)^2\right]\exp\!\left[-j2\pi f_0\,\tau_n(\mathbf{P})\right]$$

### 9.4 Phase History en Frecuencia (Stepped-Frequency o Tras Compresión)

En el dominio de la frecuencia, la contribución de un solo dispersor puntual a la Phase History es:

$$\boxed{S(f_k, t_n;\,\mathbf{P}) = g(\mathbf{P})\,A(f_k,\,t_n)\,\exp\!\left(-j\,\frac{2\pi\,f_k}{c}\,R_{OP}(\mathbf{r}_{TX,n}, \mathbf{r}_{RX,n}, \mathbf{P})\right)}$$

donde $f_k = f_0 + k\Delta f$ es la k-ésima frecuencia de la banda procesada.

La fase es por tanto:

$$\boxed{\Phi_n(\mathbf{P}) = -\frac{2\pi f_k}{c}\,R_{OP}(\mathbf{r}_{TX,n}, \mathbf{r}_{RX,n}, \mathbf{P})}$$

> Esta expresión es **exacta** dado el modelo de dos medios con Snell. Toda la información de la geometría (posiciones TX, RX, ubicación de $\mathbf{P}$, índices de refracción) está contenida en $R_{OP}$.

### 9.5 Imagen SAR por Backprojection (Suma Coherente)

Para reconstruir la imagen SAR $I(\mathbf{p})$ en un píxel de prueba $\mathbf{p}$ (que puede estar en cualquier posición del espacio):

$$\boxed{I(\mathbf{p}) = \sum_{n=1}^{N_p} \sum_{k=1}^{K} S(f_k, t_n) \cdot \exp\!\left(+j\,\frac{2\pi f_k}{c}\,R_{OP}(\mathbf{r}_{TX,n}, \mathbf{r}_{RX,n}, \mathbf{p})\right)}$$

donde $R_{OP}(\mathbf{r}_{TX,n}, \mathbf{r}_{RX,n}, \mathbf{p})$ se calcula para cada par $(n, \mathbf{p})$ resolviendo el problema de refracción descrito en §6 y §7.

---

## 10. Resumen de la Cadena Causal del Modelo

```
SISTEMA BIESTÁTICO HELICOIDAL — DOS MEDIOS (MODELO EXACTO, SIN APROXIMACIONES)
│
├── 1. TRAYECTORIAS (slow-time index n)
│     r_TX(n) = (ρ_TX·cosα_TX, ρ_TX·sinα_TX, z_TX)    [hélice TX, z > 0]
│     r_RX(n) = (ρ_RX·cosα_RX, ρ_RX·sinα_RX, z_RX)    [hélice RX, z > 0]
│     → ρ_i(t), α_i(t), z_i(t) dados por Ecs. (3.2)–(3.4)
│
├── 2. OBJETIVO
│     P = (x_P, y_P, z_P),   z_P < 0  [en el suelo]
│
├── 3. REFRACCIÓN (para cada n y cada p de imagen)
│     Q_TX* = argmin_{Q∈z=0} [|r_TX - Q|/c + |Q - P|/v₂]
│     Q_RX* = argmin_{Q∈z=0} [|P - Q|/v₂ + |Q - r_RX|/c]
│     → Condición: ley de Snell  sin θᵢ = n₂ sin θₜ  en cada lado
│     → Solución: sistema no lineal (resolución numérica en 3D general)
│
├── 4. CAMINO ÓPTICO (exacto)
│     R_OP = |r_TX - Q_TX*| + n₂|Q_TX* - P| + n₂|P - Q_RX*| + |Q_RX* - r_RX|
│     τ_n = R_OP / c   [tiempo de retraso del eco]
│
├── 5. FASE DE LA SEÑAL
│     Φ_n(P) = -(2πf/c) · R_OP
│
└── 6. FORMACIÓN DE IMAGEN (Backprojection)
      I(p) = ΣₙΣₖ S(fₖ,tₙ) · exp(+j·2πfₖ/c · R_OP(r_TX,n, r_RX,n, p))
      → Máximo coherente cuando p = P (foco perfecto)
```

---

## Tabla de Símbolos del Modelo

| Símbolo | Descripción | Unidad | Dominio |
|---------|-------------|--------|---------|
| $c$ | Velocidad de la luz en el vacío | m/s | constante |
| $\varepsilon_r$ | Permitividad relativa del suelo | — | $\varepsilon_r \geq 1$ |
| $n_2 = \sqrt{\varepsilon_r}$ | Índice de refracción del suelo | — | $n_2 \geq 1$ |
| $v_2 = c/n_2$ | Velocidad de propagación en suelo | m/s | $v_2 \leq c$ |
| $\mathbf{r}_{TX,n}$ | Posición del TX en el pulso $n$ | m | $z > 0$ |
| $\mathbf{r}_{RX,n}$ | Posición del RX en el pulso $n$ | m | $z > 0$ |
| $\rho_i(t)$ | Radio de la hélice $i$ | m | $\rho_i > 0$ |
| $\alpha_i(t)$ | Ángulo azimutal de la hélice $i$ | rad | $[0, 2\pi N_t]$ |
| $z_i(t)$ | Altura de la hélice $i$ | m | $z_i > 0$ |
| $V_0^i$ | Velocidad tangencial de la plataforma $i$ | m/s | — |
| $\beta$ | Semi-ángulo biestático | rad | $[0, \pi/2)$ |
| $\mathbf{P}$ | Posición del objetivo | m | $z_P < 0$ |
| $\mathbf{Q}_{TX,n}^*$ | Punto de refracción lado TX | m | $z=0$ |
| $\mathbf{Q}_{RX,n}^*$ | Punto de refracción lado RX | m | $z=0$ |
| $\theta_i^{TX}$ | Ángulo de incidencia lado TX | rad | $[0, \pi/2)$ |
| $\theta_t^{TX}$ | Ángulo de transmisión lado TX | rad | $[0, \theta_{crit})$ |
| $d_1^{TX}$ | Distancia TX $\to$ $\mathbf{Q}_{TX}^*$ (en aire) | m | $d_1 > 0$ |
| $d_2^{TX}$ | Distancia $\mathbf{Q}_{TX}^*$ $\to$ $\mathbf{P}$ (en suelo) | m | $d_2 > 0$ |
| $d_2^{RX}$ | Distancia $\mathbf{P}$ $\to$ $\mathbf{Q}_{RX}^*$ (en suelo) | m | $d_2 > 0$ |
| $d_1^{RX}$ | Distancia $\mathbf{Q}_{RX}^*$ $\to$ RX (en aire) | m | $d_1 > 0$ |
| $R_{OP}$ | Camino óptico total | m | $R_{OP} > 0$ |
| $\tau_n$ | Tiempo de retraso del eco | s | $\tau > 0$ |
| $f_0$ | Frecuencia portadora | Hz | — |
| $f_k$ | Frecuencias del chirp/stepped | Hz | $[f_0 - B/2, f_0 + B/2]$ |
| $B$ | Ancho de banda | Hz | $B > 0$ |
| $g(\mathbf{P})$ | Reflectividad compleja del objetivo | — | — |
| $I(\mathbf{p})$ | Imagen SAR reconstruida | — | — |

---

## Consideraciones para Extensión del Modelo

> Las siguientes extensiones se dejan para pasos posteriores (sin aproximaciones aquí):

1. **Reflexión en interfaz:** el coeficiente de transmisión de Fresnel ($T_{12}$) modifica $A(\cdot)$; para incidencia normal $T_{12} = 2n_1/(n_1+n_2)$.

2. **Ángulo crítico de reflexión total:** para $\sin\theta_i > n_1/n_2$, no hay onda transmitida al Medio 2. El ángulo crítico es $\theta_{crit} = \arcsin(n_1/n_2) = \arcsin(1/n_2)$.

3. **Dependencia en frecuencia de $n_2$:** para suelos reales, $\varepsilon_r(f)$ sigue modelos como Debye o Cole-Cole. Se requiere corrección de dispersión.

4. **Ruido de fase biestático:** la fase efectiva medida incluye ruido de los dos osciladores TX y RX independientes, aumentando la incertidumbre ~3 dB respecto al monoestático.

5. **Interfaz no plana (DTM):** $z = h(x,y)$ reemplaza el plano $z=0$; el cálculo de $\mathbf{Q}^*$ se vuelve un problema de optimización en una superficie curva.
