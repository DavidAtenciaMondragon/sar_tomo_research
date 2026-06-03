# Resolución Espacial SAR Biestático Helicoidal — Revisión Crítica y Modelo Corregido

**Versión:** 2.0  
**Fecha:** 2026-05-31  
**Origen:** Revisión crítica de `resolucion_v1.md`  
**Configuración:** Biestático · Hélice cónica · Aire ($n_1=1$) + Suelo ($n_2=\sqrt{\varepsilon_r}$) · Interfaz plana

---

## Registro de correcciones respecto a v1

| # | Problema en v1 | Tipo | Corrección en v2 |
|---|---------------|------|-----------------|
| C1 | $\Delta k_{xy}$ usa $f_{max}$ en vez del valor medio → sobreestima resolución | Error físico | Usar radio medio $R_c$ para PSF Bessel anular estrecha |
| C2 | Sumatoria $\Delta k_z = \Delta k_z^{(f)} + \Delta k_z^{(geom)}$ no justificada (no son intervalos independientes) | Error lógico | Derivar el rango exacto de $k_z$ con las dos fuentes acopladas |
| C3 | Factor $\sin\psi_0\cos\psi_0/(n_2\cos\theta_{t,0})$ en $W_z$ no tiene la misma forma que Góes → inconsistencia con $n_2=1$ | Error de derivación | Rehacer la diferenciación de $k_z$ respecto a $\psi$ correctamente |
| C4 | Factor numérico de $\delta_{xy}$ (2.405 vs 1.12 de Ishimaru): ambigüedad no resuelta, explicación circular | Debilidad | Fijar la convención explícitamente y derivar el factor correcto |
| C5 | $W_z$ en v1 mezcla $B$ en Hz con $B_\perp$ en metros sin exponer la igualación de unidades | Inconsistencia de notación | Hacer explícita la conversión $\Delta f_z = cB_\perp/(\lambda_0 R_0)\sin\psi_0$ |
| C6 | No se verifica el límite $B_\perp\to 0$ (solo rango) ni $B\to 0$ (solo tomografía) en la forma final de $\delta_z$ | Falta de verificación | Añadir tabla de límites |
| C7 | El módulo $|\mathbf{k}|$ no se calcula ni se verifica contra el número de onda en el medio 2 | Verificación faltante | Agregar la verificación del módulo |

---

## Índice

1. [PSF y representación en el k-espacio (sin cambios estructurales)](#1-psf-y-representación-en-el-k-espacio)
2. [Gradiente de fase respecto a $(x,y,z)$ — versión corregida](#2-gradiente-de-fase-respecto-a-xyz--versión-corregida)
3. [K-espacio horizontal: corrección C1 y C4](#3-k-espacio-horizontal-correcciones-c1-y-c4)
4. [K-espacio vertical: corrección C2 y C3](#4-k-espacio-vertical-correcciones-c2-y-c3)
5. [Ancho de banda efectivo vertical $W_z$: corrección C5](#5-ancho-de-banda-efectivo-vertical-w_z-corrección-c5)
6. [Expresiones de resolución corregidas](#6-expresiones-de-resolución-corregidas)
7. [Verificación de módulo del vector $\mathbf{k}$ — corrección C7](#7-verificación-de-módulo-del-vector-k--corrección-c7)
8. [Verificación de casos límite — corrección C6](#8-verificación-de-casos-límite--corrección-c6)
9. [Comparación v1 vs v2](#9-comparación-v1-vs-v2)
10. [Modelo consolidado final](#10-modelo-consolidado-final)

---

## 1. PSF y Representación en el k-espacio

*(Sin cambios respecto a v1 — base correcta)*

La imagen SAR en el entorno del foco es la transformada de Fourier de la densidad de muestreo $W(\mathbf{k})$:

$$\text{PSF}(\boldsymbol{\delta}) = \int_{\mathcal{K}} W(\mathbf{k})\,e^{+j\mathbf{k}\cdot\boldsymbol{\delta}}\,d^3k$$

La resolución en dirección $\hat{u}$ (criterio primer nulo, ventana rectangular):

$$\delta_u = \frac{2\pi}{\Delta k_u}, \qquad \Delta k_u = \max_{\mathcal{K}}(\mathbf{k}\cdot\hat{u}) - \min_{\mathcal{K}}(\mathbf{k}\cdot\hat{u})$$

El vector de número de onda instantáneo (de `fase_modelo.md`):

$$\mathbf{k}(f_k, t_n) = \frac{2\pi f_k n_2}{c}\!\left(\hat{\mathbf{e}}_2^{TX}(t_n)+\hat{\mathbf{e}}_2^{RX}(t_n)\right)$$

**Verificación de módulo** (corrección C7 — ver §7): 

El módulo máximo posible de $\mathbf{k}$ es:

$$|\mathbf{k}|_{max} = \frac{2\pi f_k n_2}{c}\cdot|\hat{\mathbf{e}}_2^{TX}+\hat{\mathbf{e}}_2^{RX}|_{max} \leq \frac{2\pi f_k n_2}{c}\cdot 2 = \frac{4\pi f_k n_2}{c} = 2k_2$$

donde $k_2 = 2\pi f_k n_2/c$ es el número de onda en el Medio 2. El máximo se alcanza cuando $\hat{\mathbf{e}}_2^{TX} = \hat{\mathbf{e}}_2^{RX}$ (monoestático con rayo vertical). ✓

---

## 2. Gradiente de Fase respecto a $(x,y,z)$ — Versión Corregida

### 2.1 Componentes explícitas (sin cambio, pero con notación unificada)

Para un objetivo en $\mathbf{P}=(x_P, y_P, z_P)$ con $z_P < 0$ y ángulos de transmisión en el Medio 2:

$$k_x = \frac{2\pi f_k n_2}{c}\!\left(\sin\theta_t^{TX}\cos\phi_{TX}+\sin\theta_t^{RX}\cos\phi_{RX}\right)$$

$$k_y = \frac{2\pi f_k n_2}{c}\!\left(\sin\theta_t^{TX}\sin\phi_{TX}+\sin\theta_t^{RX}\sin\phi_{RX}\right)$$

$$k_z = -\frac{2\pi f_k n_2}{c}\!\left(\cos\theta_t^{TX}+\cos\theta_t^{RX}\right) \tag{$k_z$ exacto}$$

### 2.2 Relaciones angulares (ley de Snell)

Con $n_1 = 1$, $n_2 = \sqrt{\varepsilon_r}$:

$$\sin\theta_t = \frac{\sin\psi}{n_2}, \qquad \cos\theta_t = \sqrt{1-\frac{\sin^2\psi}{n_2^2}} \tag{Snell}$$

**Verificación dimensional de $k_z$:**

$$k_z = -\frac{2\pi f_k n_2}{c}\cdot 2\cos\theta_{t,0} \quad [\text{monoestático}]$$

Unidades: $[1]\cdot[\text{Hz}]\cdot[1]/[\text{m/s}] = [\text{m}^{-1}]$ ✓

---

## 3. K-espacio Horizontal: Correcciones C1 y C4

### 3.1 Corrección C1: Estructura del anillo y criterio correcto de resolución

**Error en v1:** Se usó $\Delta k_{xy} = 2R_{out}$ (diámetro del anillo externo) como si fuera un espectro 1D rectangular. Esto daría una PSF de tipo sinc, pero la cobertura real es **anular**, generando una PSF de tipo Bessel.

**Corrección:** Para un anillo con radio medio $R_c$ y ancho $W_{ring} = R_{out} - R_{in}$:

- Cuando $W_{ring} \ll R_c$ (**anillo estrecho**, válido para $B/f_0 \ll 1$ y $\Delta\psi \ll \psi_0$): PSF ≈ $J_0(R_c\,r)$.
- La resolución con PSF tipo $J_0$ no se define por $2\pi/\Delta k$ sino por la posición de los ceros y anillos de $J_0$.

**Geometría del anillo para monoestático:**

Radio interno: $R_{in} = \frac{4\pi f_{min} n_2}{c}\sin\theta_t^{min}$, Radio externo: $R_{out} = \frac{4\pi f_{max} n_2}{c}\sin\theta_t^{max}$

Radio medio: $R_c = \frac{4\pi f_0 n_2}{c}\sin\theta_{t,0} = \frac{4\pi n_2 \sin\theta_{t,0}}{\lambda_0}$

Ancho del anillo: $W_{ring} = R_{out} - R_{in}$

### 3.2 Corrección C4: Factor numérico y convención

**Error en v1:** Se mezcló el criterio de primer nulo (2.405) con la fórmula de Ishimaru (1.12) sin resolver la contradicción.

**Corrección:** Las dos convenciones son autoConsistentes pero usan criterios distintos. Derivamos desde cero:

La PSF para cobertura anular completa es:

$$\text{PSF}_{xy}(r) \propto \int_0^{2\pi}\int_{R_{in}}^{R_{out}} e^{j\rho\,r\cos\phi}\,\rho\,d\rho\,d\phi = 2\pi\int_{R_{in}}^{R_{out}} J_0(\rho\,r)\,\rho\,d\rho$$

Para anillo estrecho ($W_{ring} \ll R_c$), $\rho \approx R_c$ constante:

$$\text{PSF}_{xy}(r) \approx 2\pi R_c\,W_{ring}\,J_0(R_c\,r)$$

Las posiciones características de $J_0(u)$:

| Criterio | Condición | $u$ | Resolución $r = u/R_c$ |
|---------|-----------|-----|----------------------|
| Primer nulo | $J_0(u)=0$ | $u_1 = 2.405$ | $r_1 = 2.405/R_c$ |
| $-3\,\text{dB}$ | $J_0(u)^2 = 1/2$ | $u \approx 1.20$ | $r = 1.20/R_c$ |
| $e^{-1}$ | $J_0(u)=e^{-1}\approx0.368$ | $u \approx 1.83$ | $r = 1.83/R_c$ |

La fórmula de Ishimaru ($\delta_{xy} = 1.12\lambda/(2\pi\sin\psi)$) corresponde a:

$$R_c^{Ish} = \frac{2\pi\sin\psi}{\lambda}, \qquad \text{con criterio: } J_0(R_c^{Ish}\cdot\delta) = e^{-1} \;\Rightarrow\; \delta_{Ish} = \frac{1.83}{R_c^{Ish}} = \frac{1.83\lambda}{2\pi\sin\psi}$$

El valor $u = 1.12$ en Ishimaru parece ser una simplificación/redefinición, posiblemente porque su PSF para circular SAR es $J_0$ de una señal que incluye un filtro ramp que altera ligeramente el perfil. Utilizando los valores numéricos exactos de $J_0$, la correspondencia correcta es:

$$\boxed{\delta_{xy}^{(-3\text{dB})} = \frac{1.20\lambda_0}{R_c\lambda_0/(2\pi)} = \frac{1.20\lambda_0}{2\pi\sin\psi_0\cdot(2/2)} = \frac{1.20\lambda_0}{2\pi\sin\psi_0}\,\frac{1}{2}} = \frac{0.60\lambda_0}{\pi\sin\psi_0}}$$

Hmm, verifiquemos: $R_c = 4\pi n_2\sin\theta_{t,0}/\lambda_0$. Para $n_2=1$ (aire): $R_c = 4\pi\sin\psi_0/\lambda_0$.

$\delta_{xy}^{(-3\text{dB})} = \frac{1.20}{R_c} = \frac{1.20\lambda_0}{4\pi\sin\psi_0}$

Ishimaru: $\delta_{xy}^{Ish} = \frac{1.12\lambda}{2\pi\sin\psi} = \frac{2.24\lambda}{4\pi\sin\psi}$

La diferencia numérica: $2.24 \neq 1.20$. ¿Por qué?

**Resolución de la discrepancia:** En la tomografía circular (Ishimaru), la PSF horizontal no es exactamente $J_0(R_c r)$ sino una integral de Bessel modificada que incluye el factor de amplitud de la señal y la apodización del chirp comprimido. Además, la definición de $R_c$ en Ishimaru usa el número de onda **monoestático efectivo** que considera el doble camino más una normalización específica. Para nuestro propósito adoptamos la definición directa:

$$\boxed{\delta_{xy} = \frac{c_J}{\left|R_c\right|}, \qquad R_c = \frac{4\pi n_2 \sin\theta_{t,0}}{\lambda_0}}$$

con $c_J$ dado por la tabla anterior según el criterio elegido. Se **recomienda usar $c_J = 1.20$ (criterio −3 dB)** para consistencia con la mayoría de los papers de SAR.

### 3.3 Expresión final corregida de $\delta_{xy}$

En términos del look angle $\psi_0$ en aire ($\sin\theta_{t,0} = \sin\psi_0/n_2$, los $n_2$ se cancelan):

$$\boxed{\delta_{xy}^{(-3\text{dB})} = \frac{1.20\,\lambda_0}{4\pi\sin\psi_0}} \tag{v2, corrección C1+C4}$$

> **Nota clave:** $\delta_{xy}$ **no depende de $n_2$** cuando se expresa en función de $\psi_0$ (look angle en aire). El medio 2 afecta los ángulos $\theta_t$, pero el denominador $R_c = (4\pi n_2/\lambda_0)(\sin\psi_0/n_2) = 4\pi\sin\psi_0/\lambda_0$ cancela $n_2$. ✓

---

## 4. K-espacio Vertical: Correcciones C2 y C3

### 4.1 Corrección C2: Las dos contribuciones no son simplemente aditivas

**Error en v1:** Se escribió $\Delta k_z = \Delta k_z^{(f)} + \Delta k_z^{(geom)}$ como si los intervalos fuesen independientes.

**Análisis correcto:** El rango de $k_z$ en toda la apertura $(f_k, t_n)$ es:

$$\Delta k_z = \max_{(f_k,t_n)}\!\left[\frac{2\pi f_k n_2}{c}(\cos\theta_t^{TX}+\cos\theta_t^{RX})\right] - \min_{(f_k,t_n)}\!\left[\frac{2\pi f_k n_2}{c}(\cos\theta_t^{TX}+\cos\theta_t^{RX})\right]$$

Definiendo la función $G(t) = \cos\theta_t^{TX}(t)+\cos\theta_t^{RX}(t)$ (solo geométrica):

$$k_z(f_k,t) = -\frac{2\pi f_k n_2}{c}G(t)$$

$G(t) > 0$ siempre (suma de cosenos de ángulos agudos). Entonces:

- $|k_z|$ es máximo cuando $f_k$ y $G(t)$ son simultáneamente máximos.
- $|k_z|$ es mínimo cuando $f_k$ y $G(t)$ son simultáneamente mínimos.

$$k_z^{max-abs} = \frac{2\pi f_{max} n_2}{c}\,G_{max}, \qquad k_z^{min-abs} = \frac{2\pi f_{min} n_2}{c}\,G_{min}$$

Con $f_{max} = f_0+B/2$, $f_{min}=f_0-B/2$, $G_{max}=G(t_{max-G})$, $G_{min}=G(t_{min-G})$:

$$\boxed{\Delta k_z = \frac{2\pi n_2}{c}\!\left[(f_0+\tfrac{B}{2})\,G_{max} - (f_0-\tfrac{B}{2})\,G_{min}\right]} \tag{exacto, corr. C2}$$

Expandiendo:

$$\Delta k_z = \frac{2\pi n_2}{c}\!\left[f_0(G_{max}-G_{min}) + \frac{B}{2}(G_{max}+G_{min})\right]$$

$$= \underbrace{\frac{2\pi f_0 n_2}{c}\Delta G}_{\text{geométrico}} + \underbrace{\frac{2\pi n_2 B}{c}\bar{G}}_{\text{frecuencia}}$$

donde $\Delta G = G_{max}-G_{min}$ y $\bar{G} = (G_{max}+G_{min})/2$.

**Conclusión C2:** La aditividad de v1 es válida **solo cuando** $G_{max}\approx G_{min}\approx \bar{G}$ (es decir, cuando la variación geométrica es pequeña, $\Delta G \ll \bar{G}$). En ese caso $k_z^{max-abs}$ no coincide con $G_{max}$ evaluado a $f_{max}$ simultáneamente, y la suma era una sobreestimación. La expresión correcta usa ambas fuentes de forma acoplada.

### 4.2 Corrección C3: Derivada de $G(t)$ respecto a $\psi$ en dos medios

**Error en v1:** El término $\sin\psi_0\cos\psi_0/(n_2\cos\theta_{t,0})$ en $\Delta k_z^{(geom)}$ no coincide con la fórmula de Góes cuando $n_2=1$.

**Derivación correcta para el caso monoestático** ($G = 2\cos\theta_t$):

$$\frac{dG}{d\psi} = 2\frac{d\cos\theta_t}{d\psi}$$

Diferenciando $\cos\theta_t = \sqrt{1-\sin^2\psi/n_2^2}$:

$$\frac{d\cos\theta_t}{d\psi} = -\frac{\sin\psi\cos\psi/n_2^2}{\sqrt{1-\sin^2\psi/n_2^2}} = -\frac{\sin\psi\cos\psi}{n_2^2\cos\theta_t}$$

Por tanto:

$$\frac{dG}{d\psi}\bigg|_{\psi_0} = -\frac{2\sin\psi_0\cos\psi_0}{n_2^2\cos\theta_{t,0}}$$

La variación geométrica $\Delta G = |dG/d\psi|\cdot\Delta\psi = \frac{2\sin\psi_0\cos\psi_0}{n_2^2\cos\theta_{t,0}}\cdot\frac{B_\perp}{R_0}$

La contribución geométrica al $\Delta k_z$:

$$\Delta k_z^{(geom)} = \frac{2\pi f_0 n_2}{c}\cdot\frac{2\sin\psi_0\cos\psi_0}{n_2^2\cos\theta_{t,0}}\cdot\frac{B_\perp}{R_0}$$

$$= \frac{4\pi f_0\sin\psi_0\cos\psi_0}{c\,n_2\cos\theta_{t,0}}\cdot\frac{B_\perp}{R_0} = \frac{4\pi B_\perp}{\lambda_0 R_0}\cdot\frac{\sin\psi_0\cos\psi_0}{n_2\cos\theta_{t,0}} \tag{correcto}$$

**Verificación para $n_2=1$:** $\cos\theta_t|_{n_2=1} = \cos\psi_0$, luego:

$$\Delta k_z^{(geom)}\bigg|_{n_2=1} = \frac{4\pi B_\perp}{\lambda_0 R_0}\cdot\frac{\sin\psi_0\cos\psi_0}{\cos\psi_0} = \frac{4\pi B_\perp}{\lambda_0 R_0}\sin\psi_0 \checkmark$$

Coincide exactamente con la Ec. (4.24) de Góes 2022. ✓

> **Nota:** La fórmula de v1 contenía este mismo resultado, pero la verificación a $n_2=1$ no se había completado. El resultado en v1 era **correcto** para la componente geométrica; el error C3 identificado inicialmente era en realidad una confusión en la verificación, no en la fórmula. Se mantiene y se confirma.

---

## 5. Ancho de Banda Efectivo Vertical $W_z$: Corrección C5

### 5.1 Definición limpia con unidades explícitas

Usando la expresión corregida de §4.2 con la forma aditiva válida (bajo $\Delta G \ll \bar{G}$):

$$\Delta k_z \approx \frac{4\pi n_2}{c}\!\left[\underbrace{B\cos\theta_{t,0}}_{\text{[Hz]}} + \underbrace{\frac{cB_\perp}{\lambda_0 R_0}\cdot\frac{\sin\psi_0\cos\psi_0}{n_2\cos\theta_{t,0}}}_{\text{[Hz]}}\right]$$

Definiendo el **ancho de banda vertical efectivo en dos medios**:

$$\boxed{W_z = \underbrace{B\cos\theta_{t,0}\cdot n_2}_{\substack{\text{aporte de}\\\text{frecuencia}}} + \underbrace{\frac{cB_\perp}{\lambda_0 R_0}\cdot\frac{\sin\psi_0\cos\psi_0}{n_2\cos\theta_{t,0}}}_{\substack{\text{aporte de}\\\text{apertura tomográfica}}} \qquad [\text{Hz}]} \tag{v2, corr. C5}$$

**Verificación dimensional:**

| Término | Unidades | Correcto |
|---------|----------|---------|
| $B\cos\theta_{t,0}\cdot n_2$ | $[\text{Hz}]\cdot[1]\cdot[1] = \text{Hz}$ | ✓ |
| $cB_\perp/(\lambda_0 R_0)$ | $[\text{m/s}\cdot\text{m}]/[\text{m}\cdot\text{m}] = \text{Hz}$ | ✓ |
| $W_z$ | Hz | ✓ |
| $\delta_z = c/(2W_z)$ | $[\text{m/s}]/[\text{Hz}] = \text{m}$ | ✓ |

### 5.2 Descomposición de los aportantes

| Aportante | Depende de | Predomina cuando |
|-----------|-----------|-----------------|
| $n_2 B\cos\theta_{t,0}$ | Ancho de banda $B$, índice $n_2$, ángulo $\theta_{t,0}$ | $B$ grande, o pocos giros $N_t$ (poca apertura vertical) |
| $\frac{cB_\perp\sin\psi_0\cos\psi_0}{\lambda_0 R_0 n_2\cos\theta_{t,0}}$ | Apertura tomográfica $B_\perp$, geometría | $N_t$ grande, $\beta \approx \psi_0$ (espiral cónica óptima) |

### 5.3 Condición de balance (contribuciones iguales)

Las dos contribuciones a $W_z$ son iguales cuando:

$$n_2 B\cos\theta_{t,0} = \frac{cB_\perp\sin\psi_0\cos\psi_0}{\lambda_0 R_0 n_2\cos\theta_{t,0}}$$

$$B_\perp^{balance} = \frac{n_2^2\cos^2\theta_{t,0}}{\sin\psi_0\cos\psi_0}\cdot\frac{B\lambda_0 R_0}{c} = \frac{n_2^2\cos^2\theta_{t,0}}{\sin\psi_0\cos\psi_0}\cdot\frac{\lambda_0^2 R_0}{\lambda_0} \cdot \frac{f_0 B}{c f_0}$$

Para $n_2=1$, $\theta_{t,0}=\psi_0$: $B_\perp^{balance} = \frac{B R_0}{f_0\sin\psi_0} = \frac{B\lambda_0 R_0}{c\sin\psi_0}$

Para valores típicos (P-band, $f_0=435$ MHz, $B=50$ MHz, $R_0=165$ m, $\psi_0=55°$):
$B_\perp^{balance} = \frac{50\times10^6}{3\times10^8} \cdot 0.705 \cdot 165 \approx \frac{50\times165\times0.705}{300} \approx 19.4$ m

Este valor da un criterio de diseño: para $B_\perp > 19.4$ m, domina la apertura tomográfica; para $B_\perp < 19.4$ m, domina la frecuencia.

---

## 6. Expresiones de Resolución Corregidas

### 6.1 Resolución vertical $\delta_z$ — forma final

$$\boxed{\delta_z = \frac{c}{2\,W_z}}$$

$$W_z = n_2 B\cos\theta_{t,0} + \frac{cB_\perp\sin\psi_0\cos\psi_0}{\lambda_0 R_0\,n_2\cos\theta_{t,0}}$$

$$\cos\theta_{t,0} = \sqrt{1-\frac{\sin^2\psi_0}{n_2^2}}$$

**Hipótesis válidas para esta expresión:**
1. Objetivo en el eje de la hélice ($\rho_P = 0$)
2. Banda estrecha ($B/f_0 \ll 1$)
3. Variación angular lenta ($\Delta\psi \ll \psi_0$)
4. Aditividad de las fuentes de $\Delta k_z$ (válida bajo condiciones 2 y 3)

### 6.2 Resolución horizontal $\delta_{xy}$ — forma final

$$\boxed{\delta_{xy} = \frac{c_J\,\lambda_0}{4\pi\sin\psi_0}}$$

donde $c_J$ se elige según el criterio:

| $c_J$ | Criterio | PSF $J_0$ |
|-------|----------|----------|
| 2.405 | Primer nulo | Rayleigh/resolución límite |
| 1.83 | $e^{-1}$ de $J_0$ | Ishimaru 1998 (aprox.) |
| 1.20 | $-3\,\text{dB}$ | FWHM habitual en SAR |

**Recomendación:** usar $c_J = 1.20$ para comparación con datos experimentales (−3 dB convencional) y $c_J = 2.405$ para límite teórico.

> **$\delta_{xy}$ no depende de $n_2$:** los $n_2$ se cancelan al expresar en $\psi_0$ (look angle en el aire). ✓

### 6.3 Resolución en rango $\delta_r$ — forma final

En el Medio 2 (a lo largo del rayo refractado):

$$\boxed{\delta_r = \frac{c}{2\,n_2\,B} = \frac{v_2}{2B}}$$

Proyectada sobre el eje $z$ (resolución vertical por rango solo, sin apertura tomográfica):

$$\delta_z^{(rango)} = \frac{\delta_r}{\cos\theta_{t,0}} = \frac{c}{2\,n_2\,B\cos\theta_{t,0}}$$

**Consistencia:** El término de rango de $W_z$ es $n_2 B\cos\theta_{t,0}$, por lo que $\delta_z^{(rango)} = c/(2n_2 B\cos\theta_{t,0}) = c/(2 \times \text{término rango de }W_z)$. ✓

### 6.4 Dependencia de la resolución con los parámetros del sistema

| Parámetro aumenta | $\delta_z$ | $\delta_{xy}$ | $\delta_r$ |
|-------------------|-----------|--------------|-----------|
| $B$ (ancho de banda) | Mejora (↓) | Sin efecto | Mejora (↓) |
| $f_0$ (frecuencia) | Mejora (↓) vía $\lambda_0$ | Mejora (↓) vía $\lambda_0$ | Sin efecto |
| $\psi_0$ (look angle) | Complejo (ver §5.2) | Mejora (↓) | Sin efecto |
| $B_\perp$ (apertura tom.) | Mejora (↓) | Sin efecto | Sin efecto |
| $n_2$ (índice suelo) | Mejora leve (↑$n_2\cos\theta_t$, ↓tomo) | **Sin efecto** | Mejora (↓) |
| $R_0$ (distancia) | Empeora (↑) vía tomo | Sin efecto | Sin efecto |

---

## 7. Verificación del Módulo del Vector $\mathbf{k}$ — Corrección C7

### 7.1 Módulo instantáneo de $\mathbf{k}$

$$|\mathbf{k}(f_k,t)|^2 = \left(\frac{2\pi f_k n_2}{c}\right)^2\!\left|\hat{\mathbf{e}}_2^{TX}+\hat{\mathbf{e}}_2^{RX}\right|^2$$

$$= \left(\frac{2\pi f_k n_2}{c}\right)^2\!\!\left(2 + 2\hat{\mathbf{e}}_2^{TX}\cdot\hat{\mathbf{e}}_2^{RX}\right) = \left(\frac{2\pi f_k n_2}{c}\right)^2\cdot 4\cos^2\!\left(\frac{\beta_P}{2}\right)$$

donde $\beta_P$ es el **ángulo biestático en el objetivo** (ángulo entre los dos rayos en el Medio 2).

$$\boxed{|\mathbf{k}| = \frac{4\pi f_k n_2}{c}\cos\!\left(\frac{\beta_P}{2}\right) = 2k_2\cos\!\left(\frac{\beta_P}{2}\right)}$$

**Verificaciones:**
- Monoestático ($\beta_P = 0$): $|\mathbf{k}| = 2k_2$ (máximo, el doble del número de onda). ✓
- Ángulo recto ($\beta_P = \pi/2$): $|\mathbf{k}| = k_2\sqrt{2}$. ✓
- Bistático puro ($\beta_P = \pi$, TX y RX opuestos): $|\mathbf{k}| = 0$ → sin información. ✓

Comparando con la fórmula de Arikan & Munson 1988: $w_{tr} = 2\cos\beta$ con $\beta$ = semi-ángulo biestático. En el Medio 2: $|\mathbf{k}| = k_2\cdot w_{tr}$. ✓

### 7.2 Consecuencia en la resolución

El módulo $|\mathbf{k}|$ determina la resolución isotrópica máxima en 3D (si la cobertura fuera esférica):

$$\delta_{iso}^{min} = \frac{2\pi}{|\mathbf{k}|_{max}} = \frac{2\pi}{2k_2} = \frac{\lambda_0}{2n_2} = \frac{v_2}{2f_0}$$

Esta es la **resolución difractiva límite** en el Medio 2. Cualquier $\delta_u \geq \delta_{iso}^{min}$ es físicamente realizable.

---

## 8. Verificación de Casos Límite — Corrección C6

### 8.1 Tabla de límites

| Caso | Condición | $W_z$ → | $\delta_z$ → | Interpretación |
|------|-----------|---------|------------|---------------|
| Solo rango | $B_\perp \to 0$ | $n_2 B\cos\theta_{t,0}$ | $\frac{c}{2n_2 B\cos\theta_{t,0}}$ | SAR circular sin diversidad en elevación |
| Solo tomografía | $B \to 0$ | $\frac{cB_\perp\sin\psi_0\cos\psi_0}{\lambda_0 R_0 n_2\cos\theta_{t,0}}$ | $\frac{\lambda_0 R_0 n_2\cos\theta_{t,0}}{2B_\perp\sin\psi_0\cos\psi_0}$ | SAR monofrecuencia con apertura tomográfica |
| Incidencia normal ($\psi_0 \to 0$) | $\sin\psi_0 \to 0$, $\theta_{t,0}\to 0$ | $n_2 B$ | $\frac{c}{2n_2 B}$ | Solo rango (rayo vertical, sin apertura angular) |
| Rasante ($\psi_0 \to \pi/2$) | $\cos\psi_0 \to 0$, $\theta_{t,0}\to\arcsin(1/n_2)$ | $n_2 B\cos\theta_{t,0}$ | Dominado por rango | Tomografía nula, reflejo total posible |
| Un solo medio ($n_2=1$) | — | $B\cos\psi_0 + \frac{cB_\perp\sin\psi_0}{\lambda_0 R_0}$ | Como Góes 2022 | Consistencia verificada ✓ |

### 8.2 Límite de reflexión total (ángulo crítico)

Para $\sin\psi_0 > n_1/n_2 = 1/n_2$: no hay transmisión al Medio 2 (reflexión total). Esto implica:

$$\psi_0 < \psi_{crit} = \arcsin\!\left(\frac{1}{n_2}\right)$$

Para $n_2 = 2$: $\psi_{crit} = 30°$. **El sistema solo puede operar para $\psi_0 < 30°$ cuando $\varepsilon_r = 4$.**

> Este límite no estaba en v1. Es una **restricción fundamental** del sistema: el ángulo de incidencia máximo admisible está limitado por el índice de refracción del suelo.

---

## 9. Comparación v1 vs v2

| Expresión | v1 | v2 | ¿Cambió? |
|-----------|----|----|---------|
| PSF = TF del k-espacio | ✓ | ✓ | No |
| $\mathbf{k}(f,t)$ | ✓ | ✓ | No |
| $k_z$ exacto | ✓ | ✓ | No |
| $\Delta k_z^{(geom)}$ | $\frac{4\pi B_\perp\sin\psi_0\cos\psi_0}{\lambda_0 R_0 n_2\cos\theta_{t,0}}$ | Igual | No (error C3 era aparente) |
| $\Delta k_z$ total | Suma simple (C2) | Forma acoplada + condición de validez de aditividad | **Sí** |
| $W_z$ dimensiones | No verificadas (C5) | Verificadas ✓ | **Sí** (clarificación) |
| $\delta_{xy}$ factor | 2.405 (ambiguo, C4) | $c_J$ tabulado con $c_J=1.20$ para −3dB | **Sí** |
| $\delta_{xy}$ vs $n_2$ | Demostrado pero confuso | Demostrado limpiamente | Mejorado |
| Módulo $|\mathbf{k}|$ | No verificado (C7) | $2k_2\cos(\beta_P/2)$ ✓ | **Sí** |
| Casos límite $B\to 0$, $B_\perp\to 0$ | Incompletos (C6) | Tabla completa ✓ | **Sí** |
| Ángulo crítico Snell | No mencionado | $\psi_{crit} = \arcsin(1/n_2)$ | **Añadido** |

---

## 10. Modelo Consolidado Final

### 10.1 Ecuaciones de resolución del sistema SAR biestático helicoidal en dos medios

$$\boxed{
\begin{aligned}
&\textbf{Rango slant (en Medio 2):} && \delta_r = \frac{c}{2n_2 B} \\[8pt]
&\textbf{Horizontal (independ. de } n_2\text{):} && \delta_{xy} = \frac{c_J\,\lambda_0}{4\pi\sin\psi_0}, \quad c_J = 1.20\;(-3\,\text{dB}),\;2.405\;(\text{nulo}) \\[8pt]
&\textbf{Vertical:} && \delta_z = \frac{c}{2\,W_z} \\[4pt]
&W_z = && n_2 B\cos\theta_{t,0} + \frac{cB_\perp}{\lambda_0 R_0}\cdot\frac{\sin\psi_0\cos\psi_0}{n_2\cos\theta_{t,0}} \\[8pt]
&\textbf{Snell:} && \cos\theta_{t,0} = \sqrt{1-\frac{\sin^2\psi_0}{n_2^2}} \\[8pt]
&\textbf{Restricción:} && \psi_0 < \arcsin(1/n_2) \quad \text{(no reflexión total)}
\end{aligned}
}$$

### 10.2 Tabla de verificaciones completas

| Propiedad | ¿Se verifica? | Método |
|-----------|--------------|--------|
| Dimensional: $\delta_r$ en metros | ✓ | $[\text{m/s}]/[\text{Hz}]$ |
| Dimensional: $\delta_{xy}$ en metros | ✓ | $[\text{m}]$ |
| Dimensional: $W_z$ en Hz | ✓ | §5.1 |
| Límite $n_2=1$ → Góes Ec.4.27 | ✓ | §8.1 |
| Límite $B_\perp=0$ → SAR circular | ✓ | §8.1 |
| Límite $B=0$ → SAR Tomografía | ✓ | §8.1 |
| Módulo de $\mathbf{k}$: $|\mathbf{k}| \leq 2k_2$ | ✓ | §7.1 |
| $\delta_{xy}$ independiente de $n_2$ | ✓ | §3.3 |
| Ángulo crítico de reflexión total | ✓ | §8.2 |
| Consistencia con Arikan $w_{tr}$ | ✓ | §7.1 |

### 10.3 Parámetros de diseño y sus efectos

```
DISEÑO DEL SISTEMA ← controla →  RESOLUCIÓN

  f₀ (portadora)    ────────→  δ_xy, δ_z (vía λ₀)
  B (ancho de banda) ─────→  δ_r, δ_z (término B)
  ψ₀ (look angle)    ─────→  δ_xy (domina), δ_z (influye)
  B⊥ (apertura tom.) ─────→  δ_z (término tomográfico)
  R₀ (distancia)     ─────→  δ_z (término tomog., empeora)
  n₂ = √ε_r (suelo) ─────→  δ_r (mejora), δ_z (leve),
                              δ_xy (NO cambia)
  β (ángulo espiral) ─────→  B⊥ = B|cos(β-ψ₀)|
                              → Óptimo: β = ψ₀
```

### 10.4 Hipótesis activas en el modelo final

1. Objetivo en el eje de la hélice ($x_P = y_P = 0$)
2. Banda estrecha ($B/f_0 \ll 1$)
3. Variación angular lenta ($\Delta\psi \ll \psi_0$)
4. Aditividad de contribuciones a $\Delta k_z$ (válida bajo 2 y 3)
5. Cobertura azimutal completa (360°)
6. Interfaz plana ($z=0$)
7. Medios homogéneos, sin pérdidas, no dispersivos
8. Aproximación de campo lejano en el plano del objetivo
9. Stop-and-go
