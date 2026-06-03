# Función de Fase del Sistema SAR Biestático Helicoidal con Refracción

**Versión:** 1.0  
**Fecha:** 2026-05-31  
**Referencia base:** `models/modelo_fisico_base.md`  
**Configuración:** Biestático · Hélice cónica TX y RX · Dos medios · Refracción (Snell) · Sin aproximaciones

---

## Índice

1. [Definición exacta de la función de fase](#1-definición-exacta-de-la-función-de-fase)
2. [Estructura interna de la fase: cuatro contribuciones](#2-estructura-interna-de-la-fase-cuatro-contribuciones)
3. [Gradiente de R_OP respecto al objetivo (teorema de Danskin)](#3-gradiente-de-r_op-respecto-al-objetivo-teorema-de-danskin)
4. [Gradiente de R_OP respecto a las posiciones de los sensores](#4-gradiente-de-r_op-respecto-a-las-posiciones-de-los-sensores)
5. [Vector de número de onda instantáneo en P](#5-vector-de-número-de-onda-instantáneo-en-p)
6. [Velocidad de la trayectoria helicoidal](#6-velocidad-de-la-trayectoria-helicoidal)
7. [Historia de fase: variación con el movimiento helicoidal](#7-historia-de-fase-variación-con-el-movimiento-helicoidal)
8. [Frecuencia Doppler instantánea](#8-frecuencia-doppler-instantánea)
9. [Tasa FM en azimut (segunda derivada de la fase)](#9-tasa-fm-en-azimut-segunda-derivada-de-la-fase)
10. [Propiedades de la historia de fase para SAR helicoidal](#10-propiedades-de-la-historia-de-fase-para-sar-helicoidal)
11. [Casos especiales — verificación de consistencia](#11-casos-especiales--verificación-de-consistencia)
12. [Resumen de relaciones clave](#12-resumen-de-relaciones-clave)

---

## 1. Definición Exacta de la Función de Fase

### 1.1 Fase de la Señal Recibida

A partir del modelo físico base (§9 de `modelo_fisico_base.md`), la fase acumulada por la señal en el pulso $n$ a la frecuencia $f_k$ es:

$$\boxed{\Phi(f_k, t_n;\, \mathbf{P}) = -\frac{2\pi f_k}{c}\,R_{OP}\!\left(\mathbf{r}_{TX,n},\, \mathbf{r}_{RX,n},\, \mathbf{P}\right)}$$

donde el signo negativo indica que la fase disminuye en la dirección de propagación (convenio de onda viajera $e^{-j\Phi}$), y el **camino óptico total exacto** es:

$$R_{OP} = d_1^{TX} + n_2\,d_2^{TX} + n_2\,d_2^{RX} + d_1^{RX}$$

con $d_i^{TX/RX}$ los cuatro segmentos geométricos y $n_2 = \sqrt{\varepsilon_r}$ el índice de refracción del suelo.

### 1.2 Relación con la Longitud de Onda

Expresada en términos de longitud de onda $\lambda_k = c/f_k$:

$$\Phi(f_k, t_n;\,\mathbf{P}) = -\frac{2\pi}{\lambda_k}\,R_{OP}$$

Para la frecuencia portadora central $f_0$ y longitud de onda $\lambda_0 = c/f_0$:

$$\Phi(f_0, t_n;\,\mathbf{P}) = -\frac{2\pi}{\lambda_0}\,R_{OP}$$

### 1.3 Argumentos Completos de la Función de Fase

La fase es una función de **seis grupos de variables**:

$$\Phi = \Phi\!\underbrace{\!\left(\vphantom{A^A}f_k,\right.}_{\text{frecuencia}}\underbrace{t_n,}_{\text{slow-time}}\underbrace{x_P, y_P, z_P,}_{\text{objetivo}}\underbrace{n_2,}_{\text{medio}}\underbrace{\left.\mathbf{r}_{TX,n},\, \mathbf{r}_{RX,n}\right)}_{\text{sensores}}$$

Sin embargo, $\mathbf{r}_{TX,n} = \mathbf{r}_{TX}(t_n)$ y $\mathbf{r}_{RX,n} = \mathbf{r}_{RX}(t_n)$ están completamente determinados por $t_n$ y los parámetros helicoidales, por lo que la dependencia funcional efectiva es:

$$\Phi = \Phi(f_k,\; t_n;\; \mathbf{P},\; n_2)$$

donde los parámetros de la trayectoria son constantes del diseño del sistema.

---

## 2. Estructura Interna de la Fase: Cuatro Contribuciones

### 2.1 Descomposición por Segmentos

Usando la linealidad de la fase respecto a $R_{OP}$:

$$\Phi = \Phi_1^{TX} + \Phi_2^{TX} + \Phi_2^{RX} + \Phi_1^{RX}$$

| Componente | Expresión | Medio | Descripción |
|------------|-----------|-------|-------------|
| $\Phi_1^{TX}$ | $-\dfrac{2\pi f_k}{c}\,d_1^{TX}$ | Aire | Propagación TX → $\mathbf{Q}_{TX}^*$ |
| $\Phi_2^{TX}$ | $-\dfrac{2\pi f_k\,n_2}{c}\,d_2^{TX}$ | Suelo | Propagación $\mathbf{Q}_{TX}^*$ → $\mathbf{P}$ |
| $\Phi_2^{RX}$ | $-\dfrac{2\pi f_k\,n_2}{c}\,d_2^{RX}$ | Suelo | Propagación $\mathbf{P}$ → $\mathbf{Q}_{RX}^*$ |
| $\Phi_1^{RX}$ | $-\dfrac{2\pi f_k}{c}\,d_1^{RX}$ | Aire | Propagación $\mathbf{Q}_{RX}^*$ → RX |

### 2.2 Expresión Geométrica Completa (Forma Explícita)

Sustituyendo las distancias euclidianas exactas (con $\mathbf{r}_{TX} = (x_T, y_T, z_T)^T$, $\mathbf{r}_{RX} = (x_R, y_R, z_R)^T$, $\mathbf{P} = (x_P, y_P, z_P)^T$, $\mathbf{Q}_{TX}^* = (Q_x^T, Q_y^T, 0)^T$, $\mathbf{Q}_{RX}^* = (Q_x^R, Q_y^R, 0)^T$):

$$R_{OP} =
\underbrace{\sqrt{(x_T-Q_x^T)^2+(y_T-Q_y^T)^2+z_T^2}}_{d_1^{TX}}
+ n_2\underbrace{\sqrt{(Q_x^T-x_P)^2+(Q_y^T-y_P)^2+z_P^2}}_{d_2^{TX}}$$

$$+\; n_2\underbrace{\sqrt{(x_P-Q_x^R)^2+(y_P-Q_y^R)^2+z_P^2}}_{d_2^{RX}}
+ \underbrace{\sqrt{(Q_x^R-x_R)^2+(Q_y^R-y_R)^2+z_R^2}}_{d_1^{RX}}$$

con los puntos de refracción $(Q_x^T, Q_y^T)$ y $(Q_x^R, Q_y^R)$ determinados por la ley de Snell (ver §6 y §7 del modelo físico base).

> **Dependencia implícita:** $\mathbf{Q}_{TX}^*$ y $\mathbf{Q}_{RX}^*$ dependen de $(\mathbf{r}_{TX}, \mathbf{P})$ y $(\mathbf{P}, \mathbf{r}_{RX})$ respectivamente. Esta dependencia hace que $R_{OP}$ sea una función **no analíticamente cerrada** en el caso 3D general.

---

## 3. Gradiente de $R_{OP}$ Respecto al Objetivo (Teorema de Danskin)

### 3.1 Teorema de Danskin — Enunciado Aplicado

Sea $R_{OP} = \min_{\mathbf{Q}_{TX}, \mathbf{Q}_{RX}} F(\mathbf{r}_{TX}, \mathbf{Q}_{TX}, \mathbf{P}, \mathbf{Q}_{RX}, \mathbf{r}_{RX})$ donde el mínimo se alcanza en los puntos óptimos $\mathbf{Q}_{TX}^*(\mathbf{r}_{TX}, \mathbf{P})$ y $\mathbf{Q}_{RX}^*(\mathbf{P}, \mathbf{r}_{RX})$.

El **Teorema de Danskin** (o Teorema de la Envolvente) establece que si los puntos de mínimo son únicos y la función es diferenciable en el entorno del mínimo, entonces:

$$\frac{\partial R_{OP}}{\partial \mathbf{P}} = \frac{\partial F}{\partial \mathbf{P}}\bigg|_{\mathbf{Q}_{TX}^*, \mathbf{Q}_{RX}^* \text{ fijos}}$$

es decir, la derivada de la envolvente respecto al parámetro externo $\mathbf{P}$ es igual a la derivada parcial de $F$ respecto a $\mathbf{P}$ evaluada en los puntos óptimos, **sin derivar** los propios puntos óptimos.

> **Intuición:** en el mínimo, un pequeño cambio en $\mathbf{Q}$ no cambia $F$ (en primer orden). Solo importa el efecto directo de $\mathbf{P}$ sobre $F$.

### 3.2 Cálculo del Gradiente respecto a P

Aplicando Danskin con $\mathbf{Q}_{TX}^*$ y $\mathbf{Q}_{RX}^*$ fijos:

$$\frac{\partial R_{OP}}{\partial \mathbf{P}} = n_2\frac{\partial d_2^{TX}}{\partial \mathbf{P}} + n_2\frac{\partial d_2^{RX}}{\partial \mathbf{P}}$$

Calculando cada término (la derivada de una norma euclidiana respecto a $\mathbf{P}$):

$$\frac{\partial d_2^{TX}}{\partial \mathbf{P}} = \frac{\partial}{\partial \mathbf{P}}|\mathbf{P} - \mathbf{Q}_{TX}^*| = \frac{\mathbf{P} - \mathbf{Q}_{TX}^*}{|\mathbf{P} - \mathbf{Q}_{TX}^*|} \equiv \hat{\mathbf{e}}_2^{TX}$$

$$\frac{\partial d_2^{RX}}{\partial \mathbf{P}} = \frac{\partial}{\partial \mathbf{P}}|\mathbf{P} - \mathbf{Q}_{RX}^*| = \frac{\mathbf{P} - \mathbf{Q}_{RX}^*}{|\mathbf{P} - \mathbf{Q}_{RX}^*|} \equiv \hat{\mathbf{e}}_2^{RX}$$

**Definición de vectores unitarios en el Medio 2:**

| Vector | Expresión | Descripción |
|--------|-----------|-------------|
| $\hat{\mathbf{e}}_2^{TX}$ | $\dfrac{\mathbf{P} - \mathbf{Q}_{TX}^*}{d_2^{TX}}$ | Dirección del rayo TX en el suelo (de interfaz a P) |
| $\hat{\mathbf{e}}_2^{RX}$ | $\dfrac{\mathbf{P} - \mathbf{Q}_{RX}^*}{d_2^{RX}}$ | Dirección opuesta al rayo RX en el suelo (de interfaz a P en sentido inverso) |

Entonces:

$$\boxed{\nabla_\mathbf{P}\,R_{OP} = n_2\!\left(\hat{\mathbf{e}}_2^{TX} + \hat{\mathbf{e}}_2^{RX}\right)}$$

### 3.3 Interpretación Geométrica del Gradiente en P

El vector $\nabla_\mathbf{P}\,R_{OP}$ es la **suma de los vectores unitarios de los dos rayos que llegan a P en el Medio 2**, escalada por $n_2$.

- En el caso **monoestático** ($\mathbf{Q}_{TX}^* = \mathbf{Q}_{RX}^* = \mathbf{Q}^*$, mismo punto): $\hat{\mathbf{e}}_2^{TX} = \hat{\mathbf{e}}_2^{RX}$ → $\nabla_\mathbf{P}\,R_{OP} = 2n_2\hat{\mathbf{e}}_2$.
- En el caso **biestático general**: los dos rayos llegan a P desde ángulos diferentes → la suma $\hat{\mathbf{e}}_2^{TX} + \hat{\mathbf{e}}_2^{RX}$ apunta en la **dirección bisectriz** de los dos rayos en el Medio 2.

### 3.4 Gradiente de la Fase respecto a P

Dado que $\Phi = -(2\pi f_k/c)\,R_{OP}$:

$$\boxed{\nabla_\mathbf{P}\,\Phi(f_k, t_n;\,\mathbf{P}) = -\frac{2\pi f_k\,n_2}{c}\!\left(\hat{\mathbf{e}}_2^{TX}(t_n) + \hat{\mathbf{e}}_2^{RX}(t_n)\right)}$$

Este gradiente es el **vector de número de onda instantáneo** en el espacio del objetivo (ver §5).

---

## 4. Gradiente de $R_{OP}$ Respecto a las Posiciones de los Sensores

### 4.1 Respecto a la Posición del Transmisor

Por Danskin (solo $d_1^{TX}$ depende de $\mathbf{r}_{TX}$ directamente, con $\mathbf{Q}_{TX}^*$ fijo):

$$\frac{\partial R_{OP}}{\partial \mathbf{r}_{TX}} = \frac{\partial d_1^{TX}}{\partial \mathbf{r}_{TX}} = \frac{\partial}{\partial \mathbf{r}_{TX}}|\mathbf{r}_{TX} - \mathbf{Q}_{TX}^*| = \frac{\mathbf{r}_{TX} - \mathbf{Q}_{TX}^*}{d_1^{TX}} \equiv \hat{\mathbf{e}}_1^{TX}$$

### 4.2 Respecto a la Posición del Receptor

Análogamente:

$$\frac{\partial R_{OP}}{\partial \mathbf{r}_{RX}} = \frac{\partial d_1^{RX}}{\partial \mathbf{r}_{RX}} = \frac{\mathbf{r}_{RX} - \mathbf{Q}_{RX}^*}{d_1^{RX}} \equiv \hat{\mathbf{e}}_1^{RX}$$

**Definición de vectores unitarios en el Medio 1:**

| Vector | Expresión | Descripción |
|--------|-----------|-------------|
| $\hat{\mathbf{e}}_1^{TX}$ | $\dfrac{\mathbf{r}_{TX} - \mathbf{Q}_{TX}^*}{d_1^{TX}}$ | Dirección desde $\mathbf{Q}_{TX}^*$ hacia TX (rayo TX en aire, sentido ascendente) |
| $\hat{\mathbf{e}}_1^{RX}$ | $\dfrac{\mathbf{r}_{RX} - \mathbf{Q}_{RX}^*}{d_1^{RX}}$ | Dirección desde $\mathbf{Q}_{RX}^*$ hacia RX (rayo RX en aire, sentido ascendente) |

### 4.3 Resumen de todos los gradientes de $R_{OP}$

$$\boxed{
\nabla_{\mathbf{r}_{TX}}\,R_{OP} = \hat{\mathbf{e}}_1^{TX}, \qquad
\nabla_{\mathbf{r}_{RX}}\,R_{OP} = \hat{\mathbf{e}}_1^{RX}, \qquad
\nabla_{\mathbf{P}}\,R_{OP} = n_2\!\left(\hat{\mathbf{e}}_2^{TX} + \hat{\mathbf{e}}_2^{RX}\right)
}$$

**Verificación de consistencia:** La relación entre los vectores unitarios en los dos medios está gobernada por la ley de Snell (§7 del modelo físico base):

$$\hat{\mathbf{e}}_1^{TX} \times \hat{\mathbf{n}} = \sin\theta_i^{TX}\,\hat{\mathbf{t}}^{TX}, \qquad \hat{\mathbf{e}}_2^{TX} \times \hat{\mathbf{n}} = \sin\theta_t^{TX}\,\hat{\mathbf{t}}^{TX}$$

$$\Rightarrow \quad n_1\sin\theta_i^{TX} = n_2\sin\theta_t^{TX} \quad \checkmark$$

donde $\hat{\mathbf{n}} = \hat{z}$ y $\hat{\mathbf{t}}^{TX}$ es el vector tangente a la interfaz en el plano de incidencia TX.

---

## 5. Vector de Número de Onda Instantáneo en P

### 5.1 Definición

El **vector de número de onda instantáneo en el objetivo** $\mathbf{P}$ es el gradiente de la fase de la señal (con signo invertido para seguir el convenio de física) respecto a la posición del objetivo:

$$\mathbf{k}(f_k, t_n;\,\mathbf{P}) \equiv -\nabla_\mathbf{P}\,\Phi = \frac{2\pi f_k}{c}\,\nabla_\mathbf{P}\,R_{OP}$$

$$\boxed{\mathbf{k}(f_k, t_n;\,\mathbf{P}) = \frac{2\pi f_k\,n_2}{c}\!\left(\hat{\mathbf{e}}_2^{TX}(t_n) + \hat{\mathbf{e}}_2^{RX}(t_n)\right)}$$

### 5.2 Interpretación Física

El vector $\mathbf{k}$ tiene módulo:

$$|\mathbf{k}| = \frac{2\pi f_k\,n_2}{c}\,\left|\hat{\mathbf{e}}_2^{TX} + \hat{\mathbf{e}}_2^{RX}\right| = \frac{2\pi f_k\,n_2}{c}\cdot 2\cos\!\left(\frac{\beta_{P}}{2}\right)$$

donde $\beta_P$ es el **ángulo biestático en el objetivo** (ángulo entre los dos rayos en el Medio 2 que llegan a P).

Este módulo es máximo cuando $\beta_P = 0$ (monoestático en el suelo, misma dirección de incidencia y reflexión) y nulo cuando $\beta_P = \pi$ (geometría ortogonal imposible).

### 5.3 Cobertura del Espacio de Número de Onda (k-space)

Cuando $t_n$ varía a lo largo de la trayectoria helicoidal (de $0$ a $t_{max}$) y $f_k$ varía en el ancho de banda $[f_0 - B/2,\; f_0 + B/2]$, el vector $\mathbf{k}(f_k, t_n;\,\mathbf{P})$ traza una **superficie en el k-espacio 3D**.

El volumen de k-espacio cubierto determina la resolución 3D del sistema:

$$\delta_u = \frac{2\pi}{\Delta k_u} \quad \text{en la dirección } \hat{u}$$

donde $\Delta k_u = \max_{f_k, t_n} \mathbf{k} \cdot \hat{u} - \min_{f_k, t_n} \mathbf{k} \cdot \hat{u}$ es la extensión de la cobertura en esa dirección.

---

## 6. Velocidad de la Trayectoria Helicoidal

### 6.1 Derivación de la Velocidad Exacta

Para la hélice cónica con parámetros $\{\rho(t), \alpha(t), z(t)\}$, la velocidad cartesiana se obtiene diferenciando:

$$\mathbf{r}(t) = \begin{pmatrix} \rho(t)\cos\alpha(t) \\ \rho(t)\sin\alpha(t) \\ z(t) \end{pmatrix}$$

$$\dot{\mathbf{r}}(t) = \frac{d\mathbf{r}}{dt} = \begin{pmatrix} \dot{\rho}\cos\alpha - \rho\dot{\alpha}\sin\alpha \\ \dot{\rho}\sin\alpha + \rho\dot{\alpha}\cos\alpha \\ \dot{z} \end{pmatrix}$$

Sustituyendo $\dot{\rho} = V_\rho$, $\dot{\alpha} = V_0/\rho(t)$ y $\dot{z} = V_z$:

$$\boxed{\dot{\mathbf{r}}(t) = \begin{pmatrix} V_\rho\cos\alpha(t) - V_0\sin\alpha(t) \\ V_\rho\sin\alpha(t) + V_0\cos\alpha(t) \\ V_z \end{pmatrix}}$$

### 6.2 Módulo de la Velocidad

$$\left|\dot{\mathbf{r}}(t)\right|^2 = V_\rho^2\cos^2\alpha - 2V_\rho V_0\cos\alpha\sin\alpha + V_0^2\sin^2\alpha + V_\rho^2\sin^2\alpha + 2V_\rho V_0\sin\alpha\cos\alpha + V_0^2\cos^2\alpha + V_z^2$$

$$\boxed{\left|\dot{\mathbf{r}}(t)\right| = \sqrt{V_\rho^2 + V_0^2 + V_z^2} = V_{total}}$$

El módulo de la velocidad es **constante** a lo largo de toda la hélice — resultado esperado dado que $V_0$, $V_\rho$ y $V_z$ son constantes.

### 6.3 Descomposición de la Velocidad en Componentes Físicas

$$\dot{\mathbf{r}} = \underbrace{V_0\hat{\boldsymbol{\tau}}_\alpha(t)}_{\text{tangencial}} + \underbrace{V_\rho\hat{\boldsymbol{\rho}}(t)}_{\text{radial}} + \underbrace{V_z\hat{z}}_{\text{vertical}}$$

donde los vectores locales del sistema cilíndrico son:

$$\hat{\boldsymbol{\rho}}(t) = \begin{pmatrix}\cos\alpha(t) \\ \sin\alpha(t) \\ 0\end{pmatrix}, \qquad \hat{\boldsymbol{\tau}}_\alpha(t) = \begin{pmatrix}-\sin\alpha(t) \\ \cos\alpha(t) \\ 0\end{pmatrix}, \qquad \hat{z} = \begin{pmatrix}0 \\ 0 \\ 1\end{pmatrix}$$

### 6.4 Velocidades TX y RX

Para el sistema biestático:

$$\dot{\mathbf{r}}_{TX}(t) = V_\rho^{TX}\hat{\boldsymbol{\rho}}_{TX}(t) + V_0^{TX}\hat{\boldsymbol{\tau}}_{TX}(t) + V_z^{TX}\hat{z}$$

$$\dot{\mathbf{r}}_{RX}(t) = V_\rho^{RX}\hat{\boldsymbol{\rho}}_{RX}(t) + V_0^{RX}\hat{\boldsymbol{\tau}}_{RX}(t) + V_z^{RX}\hat{z}$$

---

## 7. Historia de Fase: Variación con el Movimiento Helicoidal

### 7.1 Definición de la Historia de Fase

La **historia de fase** es la función de $t_n$ que describe cómo evoluciona la fase del eco de un objetivo $\mathbf{P}$ fijo a lo largo de la apertura sintética, a frecuencia portadora $f_0$:

$$\Psi(t;\,\mathbf{P}) \equiv \Phi(f_0, t;\,\mathbf{P}) = -\frac{2\pi f_0}{c}\,R_{OP}(\mathbf{r}_{TX}(t), \mathbf{r}_{RX}(t), \mathbf{P})$$

La historia de fase contiene toda la información sobre la **posición, velocidad y trayectoria** del objetivo relativa a los sensores.

### 7.2 Expresión Expandida de $R_{OP}(t)$

Sustituyendo las coordenadas helicoidales:

$$R_{OP}(t) = \sqrt{(\rho_{TX}\cos\alpha_{TX} - Q_x^T)^2 + (\rho_{TX}\sin\alpha_{TX} - Q_y^T)^2 + z_{TX}^2}$$

$$+ n_2\sqrt{(Q_x^T - x_P)^2 + (Q_y^T - y_P)^2 + z_P^2}$$

$$+ n_2\sqrt{(x_P - Q_x^R)^2 + (y_P - Q_y^R)^2 + z_P^2}$$

$$+ \sqrt{(Q_x^R - \rho_{RX}\cos\alpha_{RX})^2 + (Q_y^R - \rho_{RX}\sin\alpha_{RX})^2 + z_{RX}^2}$$

donde $\rho_{TX}(t) = \rho_{top}^{TX} + V_\rho^{TX} t$, $\alpha_{TX}(t) = (V_0^{TX}/V_\rho^{TX})\ln(\rho_{TX}/\rho_{top}^{TX})$, $z_{TX}(t) = z_{top}^{TX} + V_z^{TX} t$ (y análogo para RX), con los puntos de refracción $\mathbf{Q}_{TX}^*(t)$ y $\mathbf{Q}_{RX}^*(t)$ dependiendo implícitamente de $t$.

### 7.3 Periodicidad de la Historia de Fase

La historia de fase tiene una **estructura aproximadamente periódica** con el período de la hélice $T_{vuelta} = 2\pi\rho_0/V_0$. En cada vuelta completa, el ángulo $\alpha$ aumenta en $2\pi$, la posición $(x, y)$ se aproxima a la del inicio, pero la altura $z$ y el radio $\rho$ cambian en:

$$\Delta z_{vuelta} = V_z \cdot T_{vuelta}, \qquad \Delta\rho_{vuelta} = V_\rho \cdot T_{vuelta}$$

Para la hélice cilíndrica ($V_\rho = 0$, $\Delta\rho_{vuelta} = 0$), el movimiento es estrictamente periódico en $(x,y)$ pero no en $z$. Esto hace que la historia de fase sea **quasi-periódica**: igual periodo en azimut pero con deriva lenta en la componente vertical.

---

## 8. Frecuencia Doppler Instantánea

### 8.1 Derivada Temporal de $R_{OP}$ — Aplicación de Danskin

Dado que $R_{OP}$ es mínimo sobre $(\mathbf{Q}_{TX}^*, \mathbf{Q}_{RX}^*)$ para cada $t$, y estas posiciones varían con $t$ (a medida que cambian $\mathbf{r}_{TX}(t)$ y $\mathbf{r}_{RX}(t)$), aplicamos Danskin respecto a $t$:

$$\frac{dR_{OP}}{dt} = \frac{\partial R_{OP}}{\partial \mathbf{r}_{TX}}\cdot\dot{\mathbf{r}}_{TX} + \frac{\partial R_{OP}}{\partial \mathbf{r}_{RX}}\cdot\dot{\mathbf{r}}_{RX}$$

Sustituyendo los gradientes de §4:

$$\boxed{\frac{dR_{OP}}{dt} = \hat{\mathbf{e}}_1^{TX}(t)\cdot\dot{\mathbf{r}}_{TX}(t) + \hat{\mathbf{e}}_1^{RX}(t)\cdot\dot{\mathbf{r}}_{RX}(t)}$$

### 8.2 Frecuencia Doppler Instantánea

La frecuencia Doppler instantánea es la tasa de cambio de fase dividida por $-2\pi$ (en Hz):

$$\boxed{f_D(t) = -\frac{1}{\lambda_0}\frac{dR_{OP}}{dt} = -\frac{1}{\lambda_0}\!\left[\hat{\mathbf{e}}_1^{TX}\cdot\dot{\mathbf{r}}_{TX} + \hat{\mathbf{e}}_1^{RX}\cdot\dot{\mathbf{r}}_{RX}\right]}$$

### 8.3 Descomposición del Doppler Helicoidal

Usando la descomposición de la velocidad en componentes cilíndricas:

**Contribución TX:**

$$\hat{\mathbf{e}}_1^{TX}\cdot\dot{\mathbf{r}}_{TX} = \underbrace{V_0^{TX}\!\left(\hat{\mathbf{e}}_1^{TX}\cdot\hat{\boldsymbol{\tau}}_{TX}\right)}_{f_D^{az,TX}\;\text{(componente azimutal)}} + \underbrace{V_\rho^{TX}\!\left(\hat{\mathbf{e}}_1^{TX}\cdot\hat{\boldsymbol{\rho}}_{TX}\right)}_{f_D^{rad,TX}\;\text{(componente radial)}} + \underbrace{V_z^{TX}\!\left(\hat{\mathbf{e}}_1^{TX}\cdot\hat{z}\right)}_{f_D^{z,TX}\;\text{(componente vertical)}}$$

**Contribución RX:** análogo.

Por tanto la frecuencia Doppler total es:

$$f_D(t) = -\frac{1}{\lambda_0}\!\left[f_D^{az,TX} + f_D^{rad,TX} + f_D^{z,TX} + f_D^{az,RX} + f_D^{rad,RX} + f_D^{z,RX}\right]$$

### 8.4 Interpretación Geométrica de Cada Componente

| Componente | Magnitud | Descripción física |
|------------|----------|--------------------|
| $f_D^{az}$ | $V_0\,\hat{\mathbf{e}}_1\cdot\hat{\boldsymbol{\tau}}$ | Doppler azimutal: cruce del haz angular → genera la apertura sintética en $xy$ |
| $f_D^{rad}$ | $V_\rho\,\hat{\mathbf{e}}_1\cdot\hat{\boldsymbol{\rho}}$ | Doppler radial: acercamiento/alejamiento radial → **solo en espiral cónica**; cero en cilíndrica |
| $f_D^{z}$ | $V_z\,(\hat{\mathbf{e}}_1\cdot\hat{z})$ | Doppler vertical: descenso helicoidal → contribuye a la resolución vertical (apertura tomográfica) |

> La componente $f_D^z$ es especialmente importante: es la componente que genera la **diversidad de ángulo de elevación** que permite la resolución 3D. Su magnitud es $V_z\cos(\pi/2 - \theta_i) = V_z\sin\theta_i$ (donde $\theta_i$ es el ángulo de incidencia en la interfaz).

### 8.5 Ancho de Banda Doppler Total

El ancho de banda Doppler de la señal SAR helicoidal es:

$$B_D(t) = \max_{t \in [0, t_{max}]} f_D(t) - \min_{t \in [0, t_{max}]} f_D(t)$$

Para la hélice cilíndrica, la componente azimutal varía periódicamente con amplitud $\pm V_0/\lambda_0$, por lo que:

$$B_D \approx \frac{2V_0}{\lambda_0} \cdot \frac{2\rho_0\sin(\Delta\alpha_{max}/2)}{R_0}$$

donde $\Delta\alpha_{max}$ es el ángulo de apertura total iluminado.

---

## 9. Tasa FM en Azimut (Segunda Derivada de la Fase)

### 9.1 Derivada Segunda de $R_{OP}$

Diferenciando $\dot{R}_{OP} = \hat{\mathbf{e}}_1^{TX}\cdot\dot{\mathbf{r}}_{TX} + \hat{\mathbf{e}}_1^{RX}\cdot\dot{\mathbf{r}}_{RX}$ respecto a $t$:

$$\ddot{R}_{OP} = \dot{\hat{\mathbf{e}}}_1^{TX}\cdot\dot{\mathbf{r}}_{TX} + \hat{\mathbf{e}}_1^{TX}\cdot\ddot{\mathbf{r}}_{TX} + \dot{\hat{\mathbf{e}}}_1^{RX}\cdot\dot{\mathbf{r}}_{RX} + \hat{\mathbf{e}}_1^{RX}\cdot\ddot{\mathbf{r}}_{RX}$$

### 9.2 Aceleración de la Trayectoria Helicoidal

Para velocidades $V_\rho$, $V_0$, $V_z$ constantes, la aceleración es:

$$\ddot{\mathbf{r}}(t) = \frac{d}{dt}\dot{\mathbf{r}}(t) = \begin{pmatrix} -V_\rho\dot{\alpha}\sin\alpha - V_0\dot{\alpha}\cos\alpha - V_\rho\dot{\alpha}\sin\alpha \\ V_\rho\dot{\alpha}\cos\alpha - V_0\dot{\alpha}\sin\alpha + V_\rho\dot{\alpha}\cos\alpha \\ 0 \end{pmatrix}$$

Simplificando con $\dot{\alpha} = V_0/\rho$:

$$\ddot{\mathbf{r}}(t) = \frac{d}{dt}\!\begin{pmatrix}V_\rho\cos\alpha - V_0\sin\alpha \\ V_\rho\sin\alpha + V_0\cos\alpha \\ V_z\end{pmatrix} = \dot{\alpha}\begin{pmatrix}-V_\rho\sin\alpha - V_0\cos\alpha \\ V_\rho\cos\alpha - V_0\sin\alpha \\ 0\end{pmatrix} = \frac{V_0}{\rho}\begin{pmatrix}-V_\rho\sin\alpha - V_0\cos\alpha \\ V_\rho\cos\alpha - V_0\sin\alpha \\ 0\end{pmatrix}$$

### 9.3 Derivada Temporal del Vector Unitario $\hat{\mathbf{e}}_1^{TX}$

Sea $\hat{\mathbf{e}}_1^{TX} = \boldsymbol{\delta}^{TX}/d_1^{TX}$ con $\boldsymbol{\delta}^{TX} = \mathbf{r}_{TX} - \mathbf{Q}_{TX}^*$.

$$\dot{\hat{\mathbf{e}}}_1^{TX} = \frac{\dot{\boldsymbol{\delta}}^{TX}}{d_1^{TX}} - \frac{\boldsymbol{\delta}^{TX}\,\dot{d}_1^{TX}}{(d_1^{TX})^2} = \frac{1}{d_1^{TX}}\!\left(\dot{\boldsymbol{\delta}}^{TX} - \hat{\mathbf{e}}_1^{TX}\,(\hat{\mathbf{e}}_1^{TX}\cdot\dot{\boldsymbol{\delta}}^{TX})\right)$$

donde $\dot{\boldsymbol{\delta}}^{TX} = \dot{\mathbf{r}}_{TX} - \dot{\mathbf{Q}}_{TX}^*$ depende de cómo se mueve el punto de refracción con el tiempo (obtenido diferenciando el sistema Snell implícito del modelo físico base).

La tasa FM de azimut es entonces:

$$\boxed{\gamma_{az} = -\frac{2\pi f_0}{c}\ddot{R}_{OP} = -\frac{1}{\lambda_0}\frac{d^2 R_{OP}}{dt^2}}$$

Para el caso práctico de apertura pequeña (ángulo subtendido $\ll 1$ rad), se puede aproximar $\dot{\mathbf{Q}}^* \approx 0$ (el punto de refracción se mueve lentamente). En ese caso:

$$\gamma_{az} \approx -\frac{1}{\lambda_0}\!\left(\frac{|\dot{\mathbf{r}}_{TX}|^2 - (\hat{\mathbf{e}}_1^{TX}\cdot\dot{\mathbf{r}}_{TX})^2}{d_1^{TX}} + \frac{|\dot{\mathbf{r}}_{RX}|^2 - (\hat{\mathbf{e}}_1^{RX}\cdot\dot{\mathbf{r}}_{RX})^2}{d_1^{RX}}\right)$$

> La expresión exacta sin esta aproximación requiere el Jacobiano del sistema de refracción Snell, que se derivará en un paso posterior.

---

## 10. Propiedades de la Historia de Fase para SAR Helicoidal

### 10.1 Variación del Doppler con el Azimut

Para una vuelta completa de la hélice ($\alpha: 0 \to 2\pi$), la componente azimutal del Doppler varía sinusoidalmente:

$$f_D^{az}(t) \propto -V_0\,\frac{(\hat{\mathbf{e}}_1^{TX}\cdot\hat{\boldsymbol{\tau}}_{TX})}{|\lambda_0|}$$

La proyección $\hat{\mathbf{e}}_1^{TX}\cdot\hat{\boldsymbol{\tau}}_{TX}$ varía de $-\sin\psi_H$ a $+\sin\psi_H$ donde $\psi_H$ es el ángulo horizontal subtendido (relacionado con el ángulo de incidencia proyectado).

### 10.2 Efecto del Índice de Refracción en la Historia de Fase

El índice de refracción $n_2$ afecta a la historia de fase **a través de los puntos de refracción**: cuando cambia $t_n$, cambia la posición del sensor, lo que cambia $\mathbf{Q}_{TX}^*(t)$ y $\mathbf{Q}_{RX}^*(t)$ y por tanto $d_2^{TX}(t)$ y $d_2^{RX}(t)$.

La fase total puede descomponerse como:

$$\Psi(t) = -\frac{2\pi f_0}{c}\!\left[\underbrace{d_1^{TX}(t) + d_1^{RX}(t)}_{\text{varía con } t} + n_2\underbrace{\left(d_2^{TX}(t) + d_2^{RX}(t)\right)}_{\text{varía con } t \text{ (lentamente)}}\right]$$

El término $d_2^{TX} + d_2^{RX}$ varía más lentamente que $d_1^{TX} + d_1^{RX}$ porque las distancias en el suelo no cambian con la rotación azimutal (el objetivo está fijo y el punto de refracción se mueve poco con el azimut comparado con el movimiento del sensor).

### 10.3 Separación de la Fase en Componentes Lenta y Rápida

$$\Psi(t) = \underbrace{\Psi_0}_{\text{fase constante}} + \underbrace{\Psi_{az}(t)}_{\text{rápida: azimut}} + \underbrace{\Psi_z(t)}_{\text{lenta: elevación}}$$

donde:

- $\Psi_0 = -(2\pi f_0/c)(R_{0,TX} + n_2 d_{2,0}^{TX} + n_2 d_{2,0}^{RX} + R_{0,RX})$ es el valor en el centro de la apertura
- $\Psi_{az}(t)$: variación rápida debida al movimiento azimutal ($\Delta\alpha$ vuelta a vuelta) — genera la apertura en el plano $xy$
- $\Psi_z(t)$: variación lenta debida al descenso vertical ($V_z \cdot t$) — genera la **apertura tomográfica vertical** que permite resolver en $z$

### 10.4 Estructura Diferencial de la Fase entre Vueltas

Sea $\Delta\Psi_{vuelta}(t)$ la diferencia de fase entre dos posiciones de la espiral separadas verticalmente en $\Delta z_{vuelta} = V_z\,T_{vuelta}$ al mismo ángulo azimutal $\alpha$:

$$\Delta\Psi_{vuelta} \approx -\frac{2\pi f_0}{c}\,\Delta R_{OP}\bigg|_{\Delta z_{vuelta}}$$

Esta diferencia de fase es equivalente a la que se usaría en SAR Interferometría entre dos pasadas paralelas. Contiene información sobre la **altura del objetivo $z_P$** a través del número de onda vertical (ver §5 del modelo conceptual).

---

## 11. Casos Especiales — Verificación de Consistencia

### 11.1 Caso Monoestático ($\mathbf{r}_{TX} = \mathbf{r}_{RX}$, $\mathbf{Q}_{TX}^* = \mathbf{Q}_{RX}^* = \mathbf{Q}^*$)

$$R_{OP}^{mono} = 2\left(d_1^* + n_2\,d_2^*\right)$$

$$\Phi^{mono} = -\frac{4\pi f}{c}\left(d_1^* + n_2\,d_2^*\right)$$

$$\nabla_\mathbf{P}\,R_{OP}^{mono} = 2n_2\hat{\mathbf{e}}_2^*, \qquad \frac{dR_{OP}^{mono}}{dt} = 2\hat{\mathbf{e}}_1^*\cdot\dot{\mathbf{r}}$$

$$f_D^{mono} = -\frac{2}{\lambda_0}(\hat{\mathbf{e}}_1^*\cdot\dot{\mathbf{r}})$$

El factor 2 es el esperado para el monoestático (Moreira et al. 2013, Doerry 2016). ✓

### 11.2 Caso Sin Refracción ($n_2 = 1$, objetivo en $z_P = 0^-$)

Cuando $n_2 = 1$, los puntos de refracción se ubican en la proyección directa (el rayo no se dobla) y:

$$R_{OP}^{n_2=1} = d_1^{TX} + d_2^{TX} + d_2^{RX} + d_1^{RX} = |\mathbf{r}_{TX} - \mathbf{P}| + |\mathbf{P} - \mathbf{r}_{RX}|$$

que es la fórmula biestática estándar en un solo medio. ✓

### 11.3 Caso Biestático en el Vacío ($n_2 = 1$, $z_P \geq 0$)

$$R_{OP} = |\mathbf{r}_{TX} - \mathbf{P}| + |\mathbf{P} - \mathbf{r}_{RX}|$$

$$\nabla_\mathbf{P}\,R_{OP} = \hat{\mathbf{e}}_{TX\to P} + \hat{\mathbf{e}}_{RX\to P}$$

donde $\hat{\mathbf{e}}_{TX\to P} = (\mathbf{P} - \mathbf{r}_{TX})/|\mathbf{P} - \mathbf{r}_{TX}|$ y $\hat{\mathbf{e}}_{RX\to P} = (\mathbf{P} - \mathbf{r}_{RX})/|\mathbf{P} - \mathbf{r}_{RX}|$.

Este resultado coincide con la expresión del vector de scattering biestático ($\hat{\mathbf{k}}_i + \hat{\mathbf{k}}_s$ en la literatura de antenas). ✓

### 11.4 Caso SAR Circular Puro ($V_\rho = 0$, $V_z = 0$)

Para una sola vuelta circular ($V_z = 0$ → $z_{TX} = z_{RX} = z_0$ constante, $\rho = \rho_0$ constante):

- Solo la componente azimutal de la velocidad es no nula: $\dot{\mathbf{r}} = V_0\hat{\boldsymbol{\tau}}(t)$
- El Doppler es puramente azimutal: $f_D = -V_0(\hat{\mathbf{e}}_1\cdot\hat{\boldsymbol{\tau}})/\lambda_0$
- No hay componente vertical del Doppler → sin apertura en elevación → sin resolución en $z$

Esto es consistente con la literatura: SAR circular de una sola vuelta produce imágenes 2D, no 3D (Ishimaru et al. 1998, Góes 2022). ✓

---

## 12. Resumen de Relaciones Clave

### 12.1 Tabla de Gradientes y Derivadas

| Magnitud | Expresión exacta | Herramienta |
|---------|-----------------|-------------|
| $\nabla_\mathbf{P}\,R_{OP}$ | $n_2(\hat{\mathbf{e}}_2^{TX} + \hat{\mathbf{e}}_2^{RX})$ | Danskin |
| $\nabla_{\mathbf{r}_{TX}}\,R_{OP}$ | $\hat{\mathbf{e}}_1^{TX}$ | Danskin |
| $\nabla_{\mathbf{r}_{RX}}\,R_{OP}$ | $\hat{\mathbf{e}}_1^{RX}$ | Danskin |
| $dR_{OP}/dt$ | $\hat{\mathbf{e}}_1^{TX}\cdot\dot{\mathbf{r}}_{TX} + \hat{\mathbf{e}}_1^{RX}\cdot\dot{\mathbf{r}}_{RX}$ | Danskin + cadena |
| $|\dot{\mathbf{r}}|$ | $\sqrt{V_\rho^2 + V_0^2 + V_z^2}$ (constante) | Diferenciación directa |

### 12.2 Diagrama de Dependencias

```
Parámetros de la hélice: {ρ_top, ρ_base, z_top, z_base, N_t, V₀}
        │
        ▼
Trayectorias: r_TX(t), r_RX(t)   y   ṙ_TX(t), ṙ_RX(t)
        │
        │ + Objetivo P = (x_P, y_P, z_P) e índice n₂
        ▼
Puntos de refracción: Q*_TX(t), Q*_RX(t)   ← Snell / Fermat
        │
        ▼
Distancias: d₁ᵀˣ(t), d₂ᵀˣ(t), d₂ᴿˣ(t), d₁ᴿˣ(t)
        │
        ▼
Camino óptico: R_OP(t) = d₁ᵀˣ + n₂d₂ᵀˣ + n₂d₂ᴿˣ + d₁ᴿˣ
        │
        ├──→  Fase: Φ(f,t;P) = -(2πf/c)·R_OP
        │
        ├──→  Vector k en P: k = (2πf·n₂/c)(ê₂ᵀˣ + ê₂ᴿˣ)
        │
        ├──→  Doppler: f_D = -(1/λ)(ê₁ᵀˣ·ṙ_TX + ê₁ᴿˣ·ṙ_RX)
        │
        └──→  FM rate: γ_az = -(1/λ)·d²R_OP/dt²
```

### 12.3 Ecuaciones Fundamentales del Modelo de Fase

$$\boxed{
\begin{aligned}
&\textbf{Fase:} && \Phi(f_k, t;\mathbf{P}) = -\frac{2\pi f_k}{c}\,R_{OP}(\mathbf{r}_{TX}(t), \mathbf{r}_{RX}(t), \mathbf{P}) \\[6pt]
&\textbf{Camino óptico:} && R_{OP} = d_1^{TX} + n_2\,d_2^{TX} + n_2\,d_2^{RX} + d_1^{RX} \\[6pt]
&\textbf{Vector k en P:} && \mathbf{k} = \frac{2\pi f_k\,n_2}{c}\!\left(\hat{\mathbf{e}}_2^{TX} + \hat{\mathbf{e}}_2^{RX}\right) \\[6pt]
&\textbf{Doppler:} && f_D = -\frac{1}{\lambda_0}\!\left(\hat{\mathbf{e}}_1^{TX}\cdot\dot{\mathbf{r}}_{TX} + \hat{\mathbf{e}}_1^{RX}\cdot\dot{\mathbf{r}}_{RX}\right) \\[6pt]
&\textbf{Velocidad helicoidal:} && \dot{\mathbf{r}}(t) = V_\rho\hat{\boldsymbol{\rho}}(t) + V_0\hat{\boldsymbol{\tau}}(t) + V_z\hat{z}
\end{aligned}
}$$

---

## Notas para Pasos Siguientes

> Las siguientes derivaciones quedan pendientes para el próximo modelo:

1. **Jacobiano del sistema Snell** $\partial\mathbf{Q}^*/\partial t$: necesario para la expresión exacta de $\ddot{R}_{OP}$ (tasa FM).

2. **Expansión de Taylor de $R_{OP}$ alrededor del centro de apertura**: conduce a la historia de fase cuadrática y la relación con la resolución SAR.

3. **Mapa de cobertura del k-espacio 3D**: integrar el vector $\mathbf{k}(f_k, t)$ a lo largo de la apertura helicoidal para determinar $\delta_x$, $\delta_y$, $\delta_z$.

4. **Condición de foco perfecto**: $R_{OP}(\mathbf{r}_{TX}, \mathbf{r}_{RX}, \mathbf{p}) = R_{OP}(\mathbf{r}_{TX}, \mathbf{r}_{RX}, \mathbf{P})$ solo para $\mathbf{p} = \mathbf{P}$ (unicidad de la solución de la imagen SAR).
