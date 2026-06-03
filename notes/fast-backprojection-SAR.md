# Fast Backprojection Algorithm for Synthetic Aperture Radar

**Autores:** Ali F. Yegulalp  
**Año:** 1999  
**Fuente/Publicación:** Proceedings of the IEEE Radar Conference 1999, pp. 60–65. MIT Lincoln Laboratory.

---

## 1. Geometría del Sistema

Geometría monoestática general con **trayectoria de vuelo arbitraria**. El algoritmo es aplicable a cualquier forma de trayectoria (lineal, circular, helicoidal, no lineal). La antena se desplaza a lo largo de la trayectoria $\vec{q}(s)$, donde $s$ es la longitud de arco. Se procesa en modo **spotlight** (apertura completa para cada píxel), aunque también aplica a stripmap. La imagen se forma sobre un plano horizontal o sobre la superficie del terreno 3D.

El algoritmo descompone la apertura completa en **subapturas** y procesa cada subapertura sobre una rejilla polar local en lugar de una rejilla Cartesiana, explotando el bajo ancho de banda en la dirección angular de cada subapertura.

---

## 2. Ecuaciones de Resolución SAR

### Señal de datos (Phase History)

$$F_0(s, t) = \frac{1}{R^2} e^{2\pi i \nu (t - 2R/c)} \, \text{sinc}(B(t - 2R/c))$$

**Variables:**
- `s` — posición en la trayectoria (longitud de arco)
- `t` — variable temporal (fast time)
- `R = |\vec{q}(s) - \vec{d}|` — distancia antena-objetivo
- `ν` — frecuencia portadora central, $\nu = \frac{1}{2}(\nu_{min} + \nu_{max})$
- `B = ν_{max} - ν_{min}` — ancho de banda del pulso
- `c` — velocidad de la luz

### Filtro rampa (Ramp Filter)

$$\tilde{F}(s, \nu) = \tilde{F}_0(s, \nu) |\nu|$$

**Variables:**
- `\tilde{F}_0(s,\nu)` — TF de los datos respecto a $t$
- `|\nu|` — factor de ponderación (filtro rampa en frecuencia)

### Integral de Backprojection (continua)

$$I(\vec{p}) = \int F\!\left(s, \frac{2}{c}|\vec{p} - \vec{q}(s)|\right) ds$$

**Variables:**
- `\vec{p}` — posición del píxel en la imagen
- `\vec{q}(s)` — posición de la antena en el instante $s$
- `F(s,t)` — datos filtrados (tras filtro rampa)

### Integral de Backprojection (discreta)

$$I(\vec{p}) = \sum_k F\!\left(s_k, \frac{2}{c}|\vec{p} - \vec{q}(s_k)|\right) \Delta s_k$$

**Variables:**
- `s_k` — posiciones discretas en la apertura (no necesariamente uniformes)
- `\Delta s_k` — incremento de arco correspondiente

### Apertura de integración (límite por ángulo)

$$\cos(\theta_2) \leq \frac{(\vec{p} - \vec{q}(s)) \cdot \vec{q}'(s)}{|\vec{p} - \vec{q}(s)|} \leq \cos(\theta_1)$$

**Variables:**
- `\theta_1, \theta_2` — ángulos de integración (desde la trayectoria hasta el píxel)
- `\vec{q}'(s)` — vector tangente unitario a la trayectoria

### Complejidad computacional: Backprojection estándar

$$N_{op} \approx N^2 N_{pulse}$$

**Variables:**
- `N` — dimensión de la imagen ($N \times N$ píxeles)
- `N_{pulse}` — número de pulsos en la apertura

### Imagen de subapertura (Fast BP)

$$I_n(\vec{p}) = \int_{-l/2}^{l/2} F\!\left(s_n + \xi, \frac{2}{c}|\vec{p} - \vec{q}(s_n + \xi)|\right) d\xi$$

**Variables:**
- `l` — longitud de la subapertura
- `s_n = (n - \frac{1}{2})l` — centro de la n-ésima subapertura
- `\xi` — desplazamiento local dentro de la subapertura

### Subapertura en coordenadas polares

$$I_n(r, \alpha) = \int_{-l/2}^{l/2} d\xi \int_{\nu_{min}}^{\nu_{max}} d\nu \, \tilde{F}(s_n+\xi, \nu) \, e^{\frac{4\pi i\nu}{c}\sqrt{r^2 + \xi^2 - 2r\xi\alpha}}$$

**Variables:**
- `r = |\vec{p} - \vec{q}(s_n)|` — distancia radial desde el centro de la subapertura al píxel
- `\alpha = \frac{1}{r}(\vec{p}-\vec{q}(s_n))\cdot\vec{q}'(s_n)` — coseno del ángulo (coordenada angular)

### Criterio de muestreo Nyquist (polar)

$$\Delta\alpha \leq \frac{c}{2\nu_{max} \, l}$$

$$\Delta r \leq \frac{c}{2(\nu_{max} - \nu_{min})}$$

**Variables:**
- `\Delta\alpha` — paso de muestreo en la coordenada angular
- `\Delta r` — paso de muestreo radial
- `\nu_{max}` — frecuencia máxima del pulso

### Complejidad: Fast Backprojection

$$N_{op} \approx N^{3/2} N_{pulse}$$

Obtenida al elegir el tamaño óptimo de subapertura $m = \sqrt{N}$.

**Variables:**
- `N` — dimensión de imagen ($N\times N$)
- `N_{pulse}` — número de pulsos totales

---

## 3. Suposiciones del Modelo

1. El radar no se mueve significativamente durante la transmisión/recepción de un pulso (stop-and-go assumption).
2. Los datos tienen magnitud plana y fase lineal en el dominio espectral (respuesta de blanco ideal).
3. La trayectoria puede ser arbitraria, pero cada subapertura individual es suficientemente corta como para aproximarse a una línea recta a efectos del cálculo del Nyquist.
4. La longitud de la subapertura $l \ll r$ (subapertura mucho menor que la distancia al objetivo), lo que hace la relación de dispersión transparente y simplifica el análisis.
5. La imagen es 2D (plano horizontal); la extensión a 3D es directa.
6. No se modela dispersión (propagación no dispersiva, velocidad de fase constante $c$).
7. El factor de pérdida de propagación $1/R^2$ se compensa antes del procesado.

---

## 4. Notas Adicionales

- El algoritmo produce imágenes **pixel-for-pixel** idénticas al BP estándar, sin aproximaciones en la formación de la imagen. La diferencia es únicamente computacional.
- El factor de aceleración es $\approx \sqrt{N}$: para una imagen de $2000 \times 2000$ píxeles, la mejora es ~44×.
- En la práctica, con factor de sobrenmuestreo 8, la velocidad relativa al BP estándar es ~6.8× (tabla del paper).
- El máximo residuo en la imagen respecto al BP estándar es −37.5 dB, equivalente a un error de 1.3% en el valor complejo del píxel.
- Desarrollado para el programa FOPEN (Foliage Penetration) de DARPA usando el SAR CARABAS II VHF.
- La rejilla polar usada en Fast BP es en el **dominio imagen**, no en dominio frecuencia (a diferencia del Polar Format Algorithm, PFA).
