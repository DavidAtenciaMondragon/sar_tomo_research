# A Tomographic Formulation of Bistatic Synthetic Aperture Radar

**Autores:** Orhan Arikan, David C. Munson, Jr.  
**Año:** ~1988–1989 (conferencia)  
**Fuente/Publicación:** University of Illinois at Urbana-Champaign. Actas de conferencia (ICASSP o similar). Trabajo financiado por SDIO/IST, U.S. Army Research Office, contrato DAAL03-86-K-0111. Páginas 289–302 del documento.

---

## 1. Geometría del Sistema

El sistema es un SAR en **modo spotlight biestático** (BSSAR — Bistatic Spotlight-Mode SAR): el transmisor y el receptor están ubicados en **plataformas separadas**, y la antena se orienta para iluminar continuamente el mismo parche de terreno.

- El sistema de coordenadas tiene su origen en el centro del parche de terreno (ground patch).
- Posición del **transmisor:** $(x_t, y_t)$ con distancia $R_t = \sqrt{x_t^2 + y_t^2}$
- Posición del **receptor:** $(x_r, y_r)$ con distancia $R_r = \sqrt{x_r^2 + y_r^2}$
- El parche de terreno tiene radio circular $L$ ($R_t, R_r \gg L$).
- El eje $u$ es la normal al eje $v$ (dirección del frente de onda "efectivo"); los ejes $u$ y $v$ están definidos por la geometría biestática.
- **Las iso-rangos biestáticas son elipses** con focos en el transmisor y el receptor. Los puntos del parche que se encuentran sobre una misma elipse están a la misma **distancia total de propagación** (TX a punto a RX).
- El ángulo $\beta$ es el **ángulo biestático**: la mitad del ángulo entre las líneas focales desde el transmisor y desde el receptor hasta el origen. Se cumple: $\beta = 0$ para el caso monoestático.
- Para el caso 3D (no planar): el eje $u$ y los parámetros geométricos se redefinen usando las posiciones 3D $(x_t, y_t, z_t)$ y $(x_r, y_r, z_r)$; las elipses se reemplazan por **elipsoides** con focos en el transmisor y el receptor.
- La altitud del radar se toma como cero para simplicidad; el efecto de altitud no nula puede incorporarse como en el caso monoestático.

---

## 2. Ecuaciones de Resolución SAR

### Señal transmitida (chirp lineal FM)

$$s(t) = \begin{cases} \exp\!\left(j(\omega_0 t + \alpha t^2)\right), & |t| \leq T/2 \\ 0, & \text{otherwise} \end{cases}$$

**Variables:**
- `\omega_0` — frecuencia portadora (rad/s)
- `2\alpha` — tasa FM del chirp (rad/s²)
- `T` — duración del pulso (s)

---

### Señal de retorno del parche completo (integral de proyección biestática)

$$r_\theta(t) = A \int_{-L}^{L}\!\!\int_{-L}^{L} g(u\cos\theta - v\sin\theta,\; u\sin\theta + v\cos\theta)\, s\!\left(t - \tau(u,v)\right)\, du\, dv$$

**Variables:**
- `g(x,y)` — densidad de reflectividad compleja del parche de terreno (independiente del ángulo de vista y de la frecuencia)
- `\tau(u,v)` — retraso total de propagación desde el transmisor al punto $(u,v)$ al receptor (s)
- `A` — factor de atenuación de espacio libre (constante para todos los puntos del parche)
- `\theta` — ángulo de proyección (orientación del frente de onda efectivo)

---

### Retraso total de propagación biestático (aproximación)

$$\tau(u) \approx \frac{R_t + R_r + u \cdot w_{tr}}{c}$$

donde:

$$w_{tr} = 2|\cos(\beta)| = \left(2 + 2\cos(2\beta)\right)^{1/2} = \left(2 + 2\hat{R}_t \cdot \hat{R}_r\right)^{1/2}$$

y en términos de las coordenadas del transmisor y receptor:

$$w_{tr} = \left[2\!\left(1 + \frac{x_t}{R_t}\frac{x_r}{R_r} + \frac{y_t}{R_t}\frac{y_r}{R_r}\right)\right]^{1/2}$$

**Variables:**
- `R_t, R_r` — distancias desde el transmisor y el receptor al origen (m)
- `\beta` — semi-ángulo biestático (ángulo entre las líneas focales y la bisectriz)
- `\hat{R}_t, \hat{R}_r` — vectores unitarios desde el origen al transmisor y al receptor respectivamente
- `w_{tr}` — factor de escala biestático del número de onda: $w_{tr} = 2$ para el caso monoestático ($\beta=0$), y $w_{tr} < 2$ para $\beta > 0$
- Para el caso 3D: $w_{tr} = \left[2(1 + x_t x_r/(R_t R_r) + y_t y_r/(R_t R_r) + z_t z_r/(R_t R_r))\right]^{1/2}$

---

### Señal BSSAR demodulada (baseband)

Multiplicando $r_\theta(t)$ por $\mathrm{Re}(s(t))$ y $\mathrm{Im}(s(t))$, cada uno retrasado por $(R_t + R_r)/c$, y aplicando filtro paso-bajo:

$$C_\theta(t) \approx \frac{A}{2} P_\theta\!\left[\frac{w_{tr}}{c}\!\left(\omega_0 + 2\alpha\!\left(t - \frac{R_t + R_r}{c}\right)\right)\right]$$

**Variables:**
- `C_\theta(t)` — señal baseband biestática procesada (compleja)
- `P_\theta(X)` — transformada de Fourier 1D de la proyección $p_\theta(u)$ evaluada en $X$
- La señal BSSAR procesada es la **transformada de Fourier de la proyección** de $g(x,y)$, al igual que en el caso monoestático, pero escalada por el factor $w_{tr}/c$ (en lugar de $2/c$)

---

### Teorema de la proyección-rebanada para BSSAR

Por el teorema de proyección-rebanada:

$$C_\theta(t) \approx \frac{A}{2}\, G\!\left(\frac{w_{tr}}{c}\!\left[\omega_0 + 2\alpha\!\left(t-\frac{R_t+R_r}{c}\right)\right]\cos\theta,\; \frac{w_{tr}}{c}\!\left[\omega_0 + 2\alpha\!\left(t-\frac{R_t+R_r}{c}\right)\right]\sin\theta\right)$$

**Variables:**
- `G(X,Y)` — transformada de Fourier 2D de la reflectividad $g(x,y)$

---

### Rango de frecuencias espaciales adquiridas (BSSAR)

$$X_1 = \frac{w_{tr}}{c}\!\left[\omega_0 + 2\alpha\!\left(-\frac{T}{2} + \frac{Lw_{tr}}{c}\right)\right] \approx \frac{w_{tr}}{c}\left[\omega_0 - \alpha T\right]$$

$$X_2 = \frac{w_{tr}}{c}\!\left[\omega_0 + 2\alpha\!\left(\frac{T}{2} - \frac{Lw_{tr}}{c}\right)\right] \approx \frac{w_{tr}}{c}\left[\omega_0 + \alpha T\right]$$

**Variables:**
- `X_1, X_2` — frecuencias espaciales mínima y máxima adquiridas en cada rebanada radial
- Para el caso monoestático ($\beta=0$, $w_{tr}=2$): $X_1 \approx \frac{2}{c}[\omega_0 - \alpha T]$, $X_2 \approx \frac{2}{c}[\omega_0 + \alpha T]$ (coincide con Munson 1983)
- Para BSSAR con $\beta > 0$: $w_{tr} < 2$, por lo que $X_2 - X_1$ es menor, i.e., **menos información de Fourier** se adquiere en cada rebanada radial

---

### Factor de escala biestático y relación con el ángulo biestático

$$w_{tr} = 2\cos(\beta)$$

Para el caso monoestático: $\beta = 0 \Rightarrow w_{tr} = 2$. Para el caso biestático: $0 < \beta < 90° \Rightarrow 0 < w_{tr} < 2$.

---

### Error de fase por curvatura del frente de onda biestático

$$|\Delta\phi| \leq \alpha \frac{T}{c} L^2 \cdot \frac{a^5 b^5}{\left[a^4 y_0^2 + b^4 x_0^2\right]^{3/2}\!\left[a^4 - (a^2-b^2)x_0^2\right]^{1/2}}$$

donde $a$ y $b$ son los semiejes mayor y menor de la elipse de iso-rango biestática centrada en el origen con focos en el transmisor y el receptor.

**Para el caso monoestático** (con $a = (R_t + R_r)/2$):

$$|\Delta\phi| \leq \frac{\alpha T L^2}{c \cdot a}$$

La curvatura en BSSAR es menor que en MSSAR cuando el parche se encuentra en la región donde:

$$|\gamma| \leq \tan\!\left[\left(\frac{b}{a}\right)^{-3/2}\right]$$

siendo $\gamma$ el ángulo definido en el sistema de coordenadas de la elipse.

---

## 3. Suposiciones del Modelo

1. El transmisor y el receptor se encuentran en el **mismo plano** que el parche de terreno (caso 2D); el caso 3D modifica la expresión de $w_{tr}$ usando posiciones 3D y reemplaza elipses por elipsoides.
2. Las elipses de iso-rango biestáticas tienen intersecciones de segmento de línea con el parche de terreno, lo que es una aproximación válida cuando $R_t, R_r \gg L$. El error de fase introducido debe ser menor que $\pi/10$ para que la aproximación sea válida.
3. El retraso $\tau(u,v)$ es prácticamente independiente de $v$ para puntos en el parche (perpendicular al eje $u$), lo que permite reducir la integral 2D a una integral 1D (convolución con $s(t)$).
4. La **reflectividad** $g(x,y)$ es independiente del ángulo de vista y de la frecuencia (dentro del rango de operación del radar).
5. El factor de atenuación $A$ se asume constante para todos los puntos del parche (sin corrección por pérdida de propagación variable).
6. **Superposición válida:** no hay dispersión múltiple entre los objetos del parche (Born approximation).
7. El término de fase cuadrático en la señal demodulada biestática se desprecia (condición análoga a la del caso monoestático).
8. Para el análisis de curvatura del frente de onda, la elipse de iso-rango se aproxima localmente por su círculo osculador con radio de curvatura $\rho$ en el punto de interés.

---

## 4. Notas Adicionales

- Este paper extiende la formulación tomográfica del SAR monoestático (Munson, O'Brien y Jenkins, 1983 — Paper 5 de esta colección) al **caso biestático**, demostrando que el BSSAR también puede interpretarse bajo el teorema de proyección-rebanada de la CT.
- La diferencia fundamental con el caso monoestático es que el factor de escala del número de onda es $w_{tr}/c$ en lugar de $2/c$, donde $w_{tr} = 2\cos(\beta) \leq 2$.
- El ángulo biestático $\beta$ tiene un efecto directo en la **cantidad de información de Fourier** adquirida en cada rebanada radial: a mayor ángulo biestático, menor $w_{tr}$, menor rango de frecuencias espaciales, y por tanto **menor resolución efectiva** comparado con el caso monoestático con el mismo sistema.
- Se analizan **tres configuraciones de adquisición** distintas para BSSAR:
  1. **Transmisor fijo, receptor en movimiento con muestreo uniforme:** rejilla de Fourier con arcos circulares (centros en el eje X); requiere interpolación polar-a-cartesiana de alto orden.
  2. **Ambos TX y RX rotan con la misma velocidad angular:** rejilla polar uniforme; permite usar el mismo algoritmo de reconstrucción MSSAR (interpolación polar-a-cartesiana + IFFT 2D).
  3. **Rejilla trapezoidal:** reducción en 1D de la carga de interpolación; puede obtenerse con velocidades angulares desiguales o muestreo A/D variable; permite aplicar 1D IFFT + interpolación vertical.
- La curvatura del frente de onda en BSSAR puede ser **mayor o menor** que en MSSAR dependiendo de si el parche está cerca del eje mayor o menor de la elipse de iso-rango biestática.
- El paper concluye que, aunque el BSSAR tiene ventajas estratégicas operacionales (receptor encubierto), el problema de reconstrucción de imagen es igualmente difícil que en el caso monoestático; no se encontraron ventajas significativas en la forma de la rejilla de Fourier para simplificar la reconstrucción.
- Las posiciones de las muestras de Fourier adquiridas se sitúan sobre líneas que pasan por el origen del plano Fourier (condición fundamental de la CT), lo que implica que ni BSSAR ni MSSAR pueden adquirir directamente datos de Fourier cartesianos (solo radiales).
