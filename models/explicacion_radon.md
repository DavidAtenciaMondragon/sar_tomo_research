# La Transformada de Radon — Explicación Detallada

**Objetivo:** entender, paso a paso y con ecuaciones explícitas, qué es la transformada de Radon, cómo se invierte (retroproyección filtrada) y por qué es la base matemática del algoritmo de **backprojection** usado en este proyecto para reconstruir imágenes 3D a partir del historial de fase SAR.

---

## Índice

1. [Motivación: ¿por qué nos importa la transformada de Radon?](#1-motivación-por-qué-nos-importa-la-transformada-de-radon)
2. [Definición formal](#2-definición-formal)
3. [Geometría de un rayo de proyección, paso a paso](#3-geometría-de-un-rayo-de-proyección-paso-a-paso)
4. [Ejemplo más simple: la transformada de Radon de un punto](#4-ejemplo-más-simple-la-transformada-de-radon-de-un-punto)
5. [Propiedades fundamentales](#5-propiedades-fundamentales)
6. [El teorema de la rebanada de Fourier (Fourier Slice Theorem)](#6-el-teorema-de-la-rebanada-de-fourier-fourier-slice-theorem)
7. [La transformada inversa: retroproyección filtrada (FBP)](#7-la-transformada-inversa-retroproyección-filtrada-fbp)
8. [Retroproyección simple vs. filtrada — visualización](#8-retroproyección-simple-vs-filtrada--visualización)
9. [Ejemplo completo: el phantom de Shepp-Logan](#9-ejemplo-completo-el-phantom-de-shepp-logan)
10. [Conexión con el backprojection SAR de este proyecto](#10-conexión-con-el-backprojection-sar-de-este-proyecto)
11. [Resumen de ecuaciones clave](#11-resumen-de-ecuaciones-clave)

---

## 1. Motivación: ¿por qué nos importa la transformada de Radon?

En tomografía (médica, sísmica, o **SAR**), nunca medimos directamente la "imagen" $f(x,y)$ que queremos reconocer (por ejemplo, la reflectividad del subsuelo). En su lugar medimos **proyecciones**: integrales de $f$ a lo largo de líneas (rayos de un escáner CT, rayos sísmicos, o — en nuestro caso — el camino óptico de la señal de radar entre la antena y cada punto del terreno).

La **transformada de Radon** $\mathcal{R}$ es la operación matemática que convierte la imagen $f(x,y)$ en el conjunto de todas esas proyecciones. El problema de "reconstruir la imagen" es entonces el problema de **invertir** la transformada de Radon: dado el conjunto de proyecciones, recuperar $f(x,y)$.

El algoritmo de **backprojection** (retroproyección) que se usa en `simulations/proc/` para formar la imagen 3D del subsuelo es, en esencia, una versión generalizada (con refracción y geometría biestática) de la inversión de la transformada de Radon que se explica aquí. Entender Radon en su forma clásica 2D es el primer paso para entender por qué backprojection funciona y por qué necesita un **filtro** para dar una imagen nítida.

---

## 2. Definición formal

### 2.1 La idea: "sumar a lo largo de una línea"

Sea $f(x,y)$ una función 2D (la imagen, p. ej. densidad o reflectividad). Tomamos una **línea recta** $L$ en el plano y calculamos la integral de $f$ a lo largo de esa línea:

$$
p_L = \int_{L} f(x,y)\, d\ell
$$

Esto es una **única proyección**: un número que resume cuánta "masa" de $f$ atraviesa esa línea.

### 2.2 Parametrización de la línea: ángulo $\theta$ y distancia $s$

Para describir *todas* las líneas posibles del plano usamos dos parámetros:

- $\theta \in [0,\pi)$: el ángulo de la **normal** a la línea respecto al eje $x$.
- $s \in \mathbb{R}$: la distancia con signo desde el origen hasta la línea, medida a lo largo de esa normal.

La ecuación de la línea $L_{\theta,s}$ es entonces:

$$
\boxed{x\cos\theta + y\sin\theta = s}
$$

Es decir: **todos los puntos $(x,y)$ cuya proyección sobre la dirección $(\cos\theta,\sin\theta)$ vale exactamente $s$**.

La figura siguiente ilustra esta definición directamente:

![Definición de la recta L_theta_s](figures_radon/fig1b_definicion_recta.png)

**Lectura de la figura:**

- El vector $\hat{n}=(\cos\theta,\sin\theta)$ (línea azul punteada) es la **normal** a la recta, y forma un ángulo $\theta$ con el eje $x$.
- El vector $s\hat{n}$ (flecha azul sólida) marca el punto de la normal que está a distancia $s$ del origen — éste es el **pie de la perpendicular** desde el origen hasta la recta.
- La **recta roja $L_{\theta,s}$** pasa por ese punto y es **perpendicular** a $\hat{n}$ (es decir, paralela a la dirección $(-\sin\theta,\cos\theta)$).
- Los dos puntos grises $(x_1,y_1)$ y $(x_2,y_2)$ están sobre $L_{\theta,s}$, en posiciones distintas a lo largo de la recta, pero **ambos satisfacen** $x\cos\theta+y\sin\theta=s$ — la ecuación no fija un punto, sino *toda la recta*.

En otras palabras: $\theta$ fija la **orientación** de la recta (a través de su normal $\hat{n}$), y $s$ fija **a qué distancia del origen** está esa recta. Variando $(\theta,s)$ recorremos *todas* las rectas posibles del plano.

### 2.3 La integral de línea con la delta de Dirac

Para escribir la integral de línea de forma compacta (integrando sobre todo el plano $(x,y)$), usamos la función delta de Dirac $\delta(\cdot)$, que vale $0$ en todas partes salvo cuando su argumento es cero:

$$
\boxed{\,p(\theta,s) \;=\; \mathcal{R}\{f\}(\theta,s) \;=\; \int_{-\infty}^{\infty}\int_{-\infty}^{\infty} f(x,y)\,\delta(x\cos\theta+y\sin\theta - s)\,dx\,dy\,}
$$

**Interpretación de cada paso:**

1. La doble integral recorre **todo** el plano $(x,y)$.
2. La delta $\delta(x\cos\theta+y\sin\theta-s)$ actúa como un "filtro perfecto": vale $1$ (en sentido distribucional) sólo sobre los puntos que satisfacen $x\cos\theta+y\sin\theta=s$, es decir, sólo sobre la línea $L_{\theta,s}$.
3. El resultado $p(\theta,s)$ es entonces la integral de $f$ restringida a esa línea — exactamente lo que queríamos.

A la función $p(\theta,s)$, vista como imagen 2D en los ejes $(\theta, s)$, se le llama **sinograma** (porque la transformada de Radon de un punto es una curva senoidal, como veremos en la Sección 4).

---

## 3. Geometría de un rayo de proyección, paso a paso

Para visualizar la definición anterior, construimos un **sistema de coordenadas rotado** $(s,u)$, donde:

- $s$ es el eje **a lo largo** del rayo de integración.
- $u$ es el eje **perpendicular**, en la dirección en que se "barren" las distintas líneas paralelas para formar una proyección completa.

La relación entre $(x,y)$ y $(s,u)$ es una simple rotación por el ángulo $\theta$:

$$
\begin{pmatrix}u\\ s\end{pmatrix}
=
\begin{pmatrix}\cos\theta & \sin\theta\\ -\sin\theta & \cos\theta\end{pmatrix}
\begin{pmatrix}x\\ y\end{pmatrix}
\quad\Longleftrightarrow\quad
\begin{pmatrix}x\\ y\end{pmatrix}
=
\begin{pmatrix}\cos\theta & -\sin\theta\\ \sin\theta & \cos\theta\end{pmatrix}
\begin{pmatrix}u\\ s\end{pmatrix}
$$

En estas coordenadas, **una proyección completa** $p(\theta,\cdot)$ se obtiene integrando $f$ a lo largo de $s$ (manteniendo $u$ fijo) para cada valor de $u$:

$$
p(\theta,u) = \int_{-\infty}^{\infty} f(u\cos\theta - s\sin\theta,\; u\sin\theta + s\cos\theta)\,ds
$$

(Aquí renombramos la variable de barrido $u\to s$ en la notación final $p(\theta,s)$ de la Sección 2 — es la misma cantidad, sólo que la integral de la definición con la delta ya hace ese cambio de variable automáticamente.)

La figura siguiente resume toda la geometría: el objeto $f(x,y)$, la línea de integración $L_{\theta,s}$, la distancia $s_0$ desde el origen, y los ejes rotados $(s,u)$.

![Geometría de un rayo de proyección](figures_radon/fig1_geometria.png)

**Lectura de la figura:**

- El **eje azul $s$** es paralelo a la línea roja (el rayo): integramos $f$ moviéndonos a lo largo de esta dirección.
- El **eje verde $u$** es perpendicular: cada valor de $u$ (aquí, $s_0$ en la notación de la Sección 2) define una línea distinta, y al recorrer todos los valores de $u$ obtenemos **una proyección completa** $p(\theta,\cdot)$ para ese ángulo $\theta$.
- Repitiendo esto para todo $\theta\in[0,\pi)$ obtenemos el sinograma completo $p(\theta,s)$.

---

## 4. Ejemplo más simple: la transformada de Radon de un punto

El ejemplo más instructivo es calcular $\mathcal{R}\{f\}$ cuando $f$ es una **fuente puntual** desplazada del origen:

$$
f(x,y) = \delta(x-x_0)\,\delta(y-y_0)
$$

### Paso 1 — Sustituir en la definición

$$
p(\theta,s) = \int\int \delta(x-x_0)\delta(y-y_0)\,\delta(x\cos\theta+y\sin\theta-s)\,dx\,dy
$$

### Paso 2 — Las dos primeras deltas "evalúan" la integral en $(x,y)=(x_0,y_0)$

Por la propiedad de muestreo de la delta de Dirac, $\int g(x)\delta(x-x_0)dx = g(x_0)$. Aplicándola dos veces:

$$
p(\theta,s) = \delta\big(x_0\cos\theta + y_0\sin\theta - s\big)
$$

### Paso 3 — Interpretación

La proyección de un punto $(x_0,y_0)$ no es cero **únicamente** cuando

$$
\boxed{s = x_0\cos\theta + y_0\sin\theta}
$$

Es decir: en el plano $(\theta,s)$, la transformada de Radon de un punto es una **curva sinusoidal** (suma de un seno y un coseno con la misma frecuencia angular = una sinusoide pura en $\theta$, de amplitud $\sqrt{x_0^2+y_0^2}$ y fase $\arctan(y_0/x_0)$). Esto explica el nombre **sinograma**.

La figura siguiente muestra esto numéricamente: a la izquierda, un punto (suavizado a un disco pequeño) fuera del centro; a la derecha, su sinograma, con la curva analítica $s=x_0\cos\theta+y_0\sin\theta$ superpuesta en cian.

![Punto y su sinograma](figures_radon/fig2_punto_sinograma.png)

**Por qué esto importa para SAR:** un blanco puntual enterrado genera, en el dominio de la señal (frecuencia-tiempo lento), una curva característica análoga a esta sinusoide — es la "firma" que el backprojection debe reconocer y enfocar de vuelta a un punto.

---

## 5. Propiedades fundamentales

Estas propiedades se obtienen directamente de la definición y son útiles para entender cómo se comporta $\mathcal{R}$:

### 5.1 Linealidad

$$
\mathcal{R}\{a f_1 + b f_2\} = a\,\mathcal{R}\{f_1\} + b\,\mathcal{R}\{f_2\}
$$

**Demostración:** la integral es un operador lineal, y la delta de Dirac no depende de $f$, así que la linealidad de la integral se transmite directamente.

### 5.2 Desplazamiento (shift)

Si $g(x,y) = f(x-x_0,y-y_0)$, entonces

$$
\mathcal{R}\{g\}(\theta,s) = \mathcal{R}\{f\}(\theta,\, s - x_0\cos\theta - y_0\sin\theta)
$$

**Idea de la demostración:** cambio de variable $x'=x-x_0$, $y'=y-y_0$ dentro de la integral; el argumento de la delta se convierte en $x'\cos\theta+y'\sin\theta - (s - x_0\cos\theta-y_0\sin\theta)$, que es justo $\mathcal{R}\{f\}$ evaluado en $s' = s-x_0\cos\theta-y_0\sin\theta$.

Esto generaliza el resultado de la Sección 4: un punto desplazado $(x_0,y_0)$ desplaza el sinograma del punto en el origen por $x_0\cos\theta+y_0\sin\theta$.

### 5.3 Periodicidad / simetría

$$
p(\theta+\pi, s) = p(\theta, -s)
$$

**Por qué:** $\cos(\theta+\pi)=-\cos\theta$ y $\sin(\theta+\pi)=-\sin\theta$, así que la condición $x\cos(\theta+\pi)+y\sin(\theta+\pi)=s$ es equivalente a $x\cos\theta+y\sin\theta=-s$. Físicamente: ver la línea "desde el otro lado" sólo invierte el signo de $s$.

### 5.4 Derivada / escalado (mencionadas para referencia)

- **Escalado:** si $g(x,y)=f(ax,ay)$, entonces $\mathcal{R}\{g\}(\theta,s)=\frac{1}{|a|}\mathcal{R}\{f\}(\theta,as)$.
- **Rotación:** si $g$ es $f$ rotado un ángulo $\alpha$, entonces $\mathcal{R}\{g\}(\theta,s)=\mathcal{R}\{f\}(\theta-\alpha,s)$ — rotar la imagen simplemente desplaza el sinograma en $\theta$.

---

## 6. El teorema de la rebanada de Fourier (Fourier Slice Theorem)

Este es el resultado **más importante** porque conecta la transformada de Radon con la transformada de Fourier 2D, y es la llave para invertir $\mathcal{R}$.

### 6.1 Enunciado

> La transformada de Fourier 1D de una proyección $p(\theta,s)$ (respecto a $s$) es igual a una **rebanada radial** de la transformada de Fourier 2D de $f(x,y)$, tomada a lo largo de la dirección $\theta$.

Formalmente, definiendo:

$$
F(u,v) = \mathcal{F}_{2D}\{f\}(u,v) = \int\int f(x,y)\,e^{-j2\pi(ux+vy)}\,dx\,dy
$$

$$
P(\theta,\omega) = \mathcal{F}_{1D}\{p(\theta,\cdot)\}(\omega) = \int p(\theta,s)\,e^{-j2\pi\omega s}\,ds
$$

el teorema afirma:

$$
\boxed{\,P(\theta,\omega) = F(\omega\cos\theta,\;\omega\sin\theta)\,}
$$

Es decir, $P(\theta,\omega)$ es exactamente $F(u,v)$ evaluado a lo largo de la línea radial $(u,v)=(\omega\cos\theta,\omega\sin\theta)$.

### 6.2 Demostración paso a paso

**Paso 1 — Partimos de $P(\theta,\omega)$ y sustituimos la definición de $p(\theta,s)$:**

$$
P(\theta,\omega) = \int_{-\infty}^{\infty}\left[\int\int f(x,y)\,\delta(x\cos\theta+y\sin\theta-s)\,dx\,dy\right]e^{-j2\pi\omega s}\,ds
$$

**Paso 2 — Intercambiamos el orden de integración** (Fubini, asumiendo $f$ suficientemente regular) para integrar primero en $s$:

$$
P(\theta,\omega) = \int\int f(x,y)\left[\int_{-\infty}^{\infty}\delta(x\cos\theta+y\sin\theta-s)\,e^{-j2\pi\omega s}\,ds\right]dx\,dy
$$

**Paso 3 — Resolvemos la integral interna usando la propiedad de muestreo de la delta** ($\int g(s)\delta(a-s)ds=g(a)$, con $a=x\cos\theta+y\sin\theta$):

$$
\int_{-\infty}^{\infty}\delta(x\cos\theta+y\sin\theta-s)\,e^{-j2\pi\omega s}\,ds = e^{-j2\pi\omega(x\cos\theta+y\sin\theta)}
$$

**Paso 4 — Sustituimos de vuelta:**

$$
P(\theta,\omega) = \int\int f(x,y)\,e^{-j2\pi\omega(x\cos\theta+y\sin\theta)}\,dx\,dy
$$

**Paso 5 — Reconocemos la transformada de Fourier 2D.** Comparando con la definición de $F(u,v)$, vemos que el exponente $-j2\pi(ux+vy)$ coincide si identificamos:

$$
u = \omega\cos\theta, \qquad v=\omega\sin\theta
$$

Por lo tanto:

$$
P(\theta,\omega) = F(\omega\cos\theta,\;\omega\sin\theta) \qquad \blacksquare
$$

### 6.3 Visualización

La figura siguiente muestra los tres objetos del teorema para el phantom de Shepp-Logan, a un ángulo $\theta_0=35^\circ$:

1. **Izquierda:** la proyección $p(\theta_0,s)$ (un perfil 1D).
2. **Centro:** su transformada de Fourier 1D, $|P(\theta_0,\omega)|$.
3. **Derecha:** la transformada de Fourier 2D de la imagen completa, $|F(u,v)|$ (en escala log), con la línea radial a $\theta_0$ marcada en rojo — el teorema dice que el centro coincide exactamente con los valores de $|F(u,v)|$ sobre esa línea.

![Teorema de la rebanada de Fourier](figures_radon/fig4_fourier_slice.png)

**Consecuencia clave:** si tenemos proyecciones para **todos** los ángulos $\theta\in[0,\pi)$, entonces — vía sus FFT 1D — tenemos $F(u,v)$ en **todo el plano de frecuencias** (en coordenadas polares). Recuperar $f(x,y)$ es entonces, en principio, tan simple como aplicar una FFT 2D inversa a $F(u,v)$. El problema práctico es que las muestras de $F(u,v)$ así obtenidas están en una **rejilla polar**, no cartesiana — de ahí surge la necesidad del filtro que veremos a continuación.

---

## 7. La transformada inversa: retroproyección filtrada (FBP)

### 7.1 Punto de partida: Fourier inversa 2D en coordenadas polares

La transformada de Fourier 2D inversa es:

$$
f(x,y) = \int_{-\infty}^{\infty}\int_{-\infty}^{\infty} F(u,v)\,e^{j2\pi(ux+vy)}\,du\,dv
$$

Cambiamos a coordenadas polares en el plano de frecuencias: $u=\omega\cos\theta$, $v=\omega\sin\theta$, con Jacobiano $du\,dv = |\omega|\,d\omega\,d\theta$, y $\theta\in[0,\pi)$, $\omega\in(-\infty,\infty)$ (permitir $\omega$ negativo cubre el círculo completo con $\theta\in[0,\pi)$ en vez de $[0,2\pi)$):

$$
f(x,y) = \int_0^{\pi}\int_{-\infty}^{\infty} F(\omega\cos\theta,\omega\sin\theta)\,e^{j2\pi\omega(x\cos\theta+y\sin\theta)}\,|\omega|\,d\omega\,d\theta
$$

### 7.2 Aplicar el teorema de la rebanada de Fourier

Por el resultado de la Sección 6, $F(\omega\cos\theta,\omega\sin\theta)=P(\theta,\omega)$. Sustituyendo:

$$
f(x,y) = \int_0^{\pi}\left[\int_{-\infty}^{\infty} P(\theta,\omega)\,|\omega|\,e^{j2\pi\omega(x\cos\theta+y\sin\theta)}\,d\omega\right]d\theta
$$

### 7.3 Identificar el filtro y la retroproyección

Definimos $s = x\cos\theta+y\sin\theta$ (la misma $s$ de la Sección 2 — ¡es la coordenada del punto $(x,y)$ proyectada sobre la dirección $\theta$!). La integral interna es entonces:

$$
q(\theta,s) \;=\; \int_{-\infty}^{\infty} \underbrace{P(\theta,\omega)\,|\omega|}_{\text{proyección filtrada en frecuencia}}\,e^{j2\pi\omega s}\,d\omega \;=\; \mathcal{F}^{-1}_{1D}\big\{P(\theta,\omega)|\omega|\big\}(s)
$$

Es decir, $q(\theta,s)$ es la **proyección filtrada**: tomamos cada proyección $p(\theta,s)$, la pasamos al dominio de Fourier, la multiplicamos por $|\omega|$ (el **filtro rampa**), y volvemos al dominio espacial.

La fórmula de inversión completa queda:

$$
\boxed{\,f(x,y) = \int_0^{\pi} q\big(\theta,\;x\cos\theta+y\sin\theta\big)\,d\theta\,}
$$

### 7.4 Lectura del algoritmo (FBP en dos pasos)

1. **Filtrar:** para cada ángulo $\theta$, calcular $q(\theta,s) = \mathcal{F}^{-1}_{1D}\{|\omega|\cdot \mathcal{F}_{1D}\{p(\theta,s)\}\}$. Esto realza altas frecuencias (bordes) y atenúa el "borroneo" de bajas frecuencias.
2. **Retroproyectar:** para cada píxel $(x,y)$ de la imagen de salida, y para cada ángulo $\theta$, calcular $s=x\cos\theta+y\sin\theta$, leer (interpolando) el valor $q(\theta,s)$, y **sumarlo** a $f(x,y)$. Repetir para todos los $\theta$ e integrar (sumar).

Este es el famoso algoritmo de **retroproyección filtrada (Filtered Back-Projection, FBP)**, el estándar en CT médico y la base conceptual del backprojection SAR.

### 7.5 El filtro rampa $|\omega|$

El factor $|\omega|$ proviene directamente del Jacobiano del cambio a coordenadas polares (paso 7.1) — **no es un ajuste arbitrario**, es una consecuencia geométrica de que las muestras de frecuencia obtenidas del teorema de la rebanada están más densas cerca del origen ($\omega$ pequeño) que lejos. El filtro $|\omega|$ compensa esa sobre-densidad.

En la práctica, $|\omega|$ crece sin límite, así que se trunca a una frecuencia de corte y se suaviza con una ventana (Hann, Hamming, etc.) para reducir ruido de alta frecuencia:

![Filtro rampa](figures_radon/fig5_filtro_rampa.png)

---

## 8. Retroproyección simple vs. filtrada — visualización

### 8.1 ¿Qué pasa si **no** filtramos? (sólo el paso 2 de la Sección 7.4)

Si simplemente sumamos las proyecciones sin pasar por el filtro $|\omega|$ (es decir, hacemos $q(\theta,s)=p(\theta,s)$), obtenemos la **retroproyección simple (laminograma)**:

$$
f_{\text{simple}}(x,y) = \int_0^{\pi} p\big(\theta,\;x\cos\theta+y\sin\theta\big)\,d\theta
$$

Esto **no** recupera $f(x,y)$ exactamente. Se puede demostrar que:

$$
f_{\text{simple}}(x,y) = f(x,y) * h(x,y), \qquad h(x,y)=\frac{1}{\sqrt{x^2+y^2}}
$$

es decir, la imagen verdadera **convolucionada con un núcleo $1/r$** que la difumina (efecto de "halo" o borrosidad radial). Esto es exactamente lo que corrige el filtro $|\omega|$: en frecuencia, $\mathcal{F}\{1/r\} \propto 1/|\omega|$, así que multiplicar por $|\omega|$ **cancela** esa convolución.

### 8.2 Comparación numérica

Para el phantom de Shepp-Logan, comparamos la imagen original, la retroproyección simple (sin filtro) y la retroproyección filtrada:

![Reconstrucción: simple vs filtrada](figures_radon/fig6_reconstruccion.png)

La versión sin filtrar es una mancha borrosa (cada punto se "esparce" como $1/r$); la versión filtrada recupera correctamente los bordes y estructuras del phantom.

### 8.3 Cómo se "construye" la imagen: superposición de "estrellas"

Otra forma de ver la retroproyección simple: cada proyección, al ser "untada" de vuelta sobre la imagen a lo largo de su línea original, genera una franja. Sumando franjas de muchos ángulos sobre un punto, se forma un patrón tipo **estrella** que se concentra cada vez más en el punto verdadero a medida que aumentan los ángulos — pero **nunca** converge a una delta perfecta sin el filtro (queda el halo $1/r$):

![Construcción por superposición de proyecciones](figures_radon/fig7_retroproyeccion_acumulada.png)

Con 1 proyección se ve una franja; con 4, una estrella de 4 puntas; con 16, ya casi un círculo borroso; con 180, el halo $1/r$ es evidente alrededor del punto central.

---

## 9. Ejemplo completo: el phantom de Shepp-Logan

El **phantom de Shepp-Logan** es la imagen de prueba estándar en tomografía (simula un corte de cráneo con elipses de distinta densidad). Aplicando $\mathcal{R}$ para $\theta\in[0^\circ,180^\circ)$ obtenemos el sinograma:

![Phantom y su sinograma](figures_radon/fig3_phantom_sinograma.png)

Cada **fila horizontal** del sinograma (a un $\theta$ fijo) es una proyección 1D de la imagen vista desde ese ángulo. Cada **columna** (a un $s$ fijo) muestra cómo varía, con el ángulo, la integral a lo largo de las líneas que pasan por esa distancia $s$ del centro. Las curvas sinusoidales visibles son la superposición de las "firmas" (Sección 4) de cada una de las elipses que componen el phantom — cada elipse, vista como una colección de puntos, contribuye con su propio conjunto de sinusoides desplazadas.

---

## 10. Conexión con el backprojection SAR de este proyecto

El algoritmo de imagen 3D usado en `simulations/proc/` (ver `models/modelo_fisico_base.md` y `models/fase_modelo.md`) es conceptualmente el **mismo principio de retroproyección**, generalizado de tres maneras:

| Radon clásico (CT) | Backprojection SAR (este proyecto) |
|---|---|
| Líneas rectas $L_{\theta,s}$ | Caminos ópticos refractados $R_{OP}(\mathbf{P})$ entre TX, interfaz, blanco $\mathbf{P}$ y RX (con ley de Snell) |
| Proyección $p(\theta,s)$ = integral de $f$ sobre la línea | Historial de fase $S(f,t)=A\,e^{-j2\pi f R_{OP}/c}$, medido en frecuencia $f$ y tiempo lento $t$ |
| Parámetro de barrido: ángulo $\theta\in[0,\pi)$ | Trayectoria helicoidal del sensor (apertura sintética 3D) |
| Retroproyección: sumar $p(\theta, x\cos\theta+y\sin\theta)$ sobre $\theta$ | Para cada vóxel $\mathbf{P}$ de la grilla 3D, sumar (coherentemente) $S(f,t)$ evaluado en la fase correspondiente a $R_{OP}(\mathbf{P})$, sobre todo $(f,t)$ |
| Filtro rampa $\lvert\omega\rvert$ corrige el Jacobiano polar→cartesiano | La densidad no uniforme de cobertura en el "k-espacio" (ver `derivacion_modelo_resolucion.md`, Sección 6-7) juega un papel análogo — determina la forma de la PSF |
| Resultado: imagen 2D $f(x,y)$ | Resultado: imagen 3D de reflectividad del subsuelo |

La idea central que **se traslada directamente** es:

> **Cada medición (proyección, o muestra de fase) "vota" a lo largo de un lugar geométrico (línea, o superficie isócrona) en el dominio de la imagen. Sumar coherentemente todos los votos para cada punto de la imagen reconstruye $f$ — pero el resultado está "borroneado" por la cobertura no uniforme del dominio dual (frecuencia/k-espacio), lo cual determina la forma de la PSF (Función de Dispersión de Punto) y, por extensión, la resolución del sistema.**

Esto es exactamente el enfoque de `derivacion_modelo_resolucion.md`: en vez de derivar un filtro $|\omega|$ explícito como en FBP, ese documento analiza directamente la **cobertura del k-espacio** $\mathbf{k}(f,t)$ generada por la trayectoria helicoidal y la refracción, y de ahí obtiene la PSF como transformada de Fourier de esa cobertura — el análogo 3D y biestático del teorema de la rebanada de Fourier presentado aquí.

---

## 11. Resumen de ecuaciones clave

| # | Ecuación | Significado |
|---|---|---|
| 1 | $p(\theta,s)=\displaystyle\int\int f(x,y)\,\delta(x\cos\theta+y\sin\theta-s)\,dx\,dy$ | Definición de la transformada de Radon (proyección a lo largo de $L_{\theta,s}$) |
| 2 | $x\cos\theta+y\sin\theta=s$ | Ecuación de la línea de proyección |
| 3 | $\mathcal{R}\{\delta(x-x_0,y-y_0)\}=\delta(s-x_0\cos\theta-y_0\sin\theta)$ | Sinograma de un punto: curva sinusoidal |
| 4 | $P(\theta,\omega)=F(\omega\cos\theta,\omega\sin\theta)$ | Teorema de la rebanada de Fourier |
| 5 | $f(x,y)=\displaystyle\int_0^{\pi} q(\theta,x\cos\theta+y\sin\theta)\,d\theta$ | Fórmula de inversión (FBP) |
| 6 | $q(\theta,s)=\mathcal{F}^{-1}_{1D}\{\,\lvert\omega\rvert\,\mathcal{F}_{1D}\{p(\theta,s)\}\,\}$ | Proyección filtrada (filtro rampa $\lvert\omega\rvert$) |
| 7 | $f_{\text{simple}}=f * \dfrac{1}{\sqrt{x^2+y^2}}$ | Retroproyección sin filtrar = imagen verdadera difuminada por núcleo $1/r$ |

---

## Reproducibilidad

Todas las figuras de este documento se generan con:

```bash
python utils/generate_figures_radon.py
```

que guarda los archivos PNG en `models/figures_radon/`, usando `numpy`, `matplotlib` y `scikit-image` (`skimage.transform.radon`/`iradon`, `skimage.data.shepp_logan_phantom`).
