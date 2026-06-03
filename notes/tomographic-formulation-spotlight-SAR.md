# A Tomographic Formulation of Spotlight-Mode Synthetic Aperture Radar

**Autores:** David C. Munson, Jr., James Dennis O'Brien, W. Kenneth Jenkins  
**Año:** 1983  
**Fuente/Publicación:** Proceedings of the IEEE, Vol. 71, No. 8, August 1983, pp. 917–925. DOI: 0018-9219/83/0800-0917$01.00

---

## 1. Geometría del Sistema

El sistema es un SAR en **modo spotlight monoestático** (MSSAR): el haz de la antena radar se orienta continuamente hacia un parche de terreno fijo mientras la aeronave se desplaza a lo largo de la trayectoria de vuelo.

- El sistema de coordenadas tiene su origen en el centro del parche de terreno, con eje $x$ (acimut) y eje $y$ (rango). La coordenada $u$ del eje de proyección es perpendicular al eje $v$ (dirección de integración) y forma el ángulo $\theta$ con el eje $x$.
- El radar transmite desde diferentes ángulos $\theta$ a medida que la aeronave se mueve, registrando diferentes vistas del parche de terreno.
- La distancia desde el radar hasta el centro del parche es $R$ (generalmente $R \gg L$, la mitad del tamaño del parche).
- Para el caso de **depresión cero** (ground plane geometry): los puntos equidistantes del radar se sitúan en un arco, aproximado por una línea recta cuando $R \gg L$.
- Para geometría de **plano slant** (altura $h \neq 0$): se introduce una corrección de estiramiento en la dimensión de rango por el factor $R_0/R$ (donde $R_0$ es la distancia slant real y $R$ la proyección al suelo).

**Conexión con tomografía:** El SAR spotlight es interpretado como un sistema tomográfico: cada posición de la aeronave produce una **proyección** del parche de terreno bajo el **teorema de la proyección-rebanada** (projection-slice theorem). La imagen SAR se reconstruye con algoritmos análogos a los de CT.

---

## 2. Ecuaciones de Resolución SAR

### Señal transmitida (chirp lineal FM)

$$s(t) = \begin{cases} e^{j(\omega_0 t + \alpha t^2)}, & |t| \leq T/2 \\ 0, & \text{otherwise} \end{cases}$$

**Variables:**
- `\omega_0` — frecuencia portadora RF (rad/s)
- `2\alpha` — tasa FM (FM rate, rad/s²); el ancho de banda del chirp es $2\alpha T$
- `T` — duración del pulso (s)

---

### Señal de retorno de un diferencial de reflectividad en $(x_0, y_0)$ a distancia $R_0$

$$r_0(t) = A \cdot g(x_0, y_0) \cdot s\!\left(t - \frac{2R_0}{c}\right)$$

**Variables:**
- `A` — factor de atenuación de propagación (considerado constante para el parche)
- `g(x_0, y_0)` — densidad de reflectividad compleja del parche de terreno; amplitud $|g| = |g(x,y)|$ escalada por $|\text{sinusoide reflejado}|$, fase $\angle g(x_0, y_0)$
- `R_0` — distancia desde la antena al punto $(x_0, y_0)$
- `c` — velocidad de la luz (m/s)

---

### Señal de retorno integrada del parche completo (proyección)

$$r_\theta(t) = A \cdot \mathrm{Re}\!\left\{\int_{-L}^{L} p_\theta(u)\, s\!\left(t - \frac{2(R + u)}{c}\right) du\right\}$$

**Variables:**
- `p_\theta(u)` — proyección de la reflectividad $g(x,y)$ al ángulo $\theta$: integral de línea en la dirección $v$ (perpendicular a $u$)
- `R` — distancia al centro del parche
- `u` — coordenada en la dirección de proyección (a lo largo del frente de onda)
- `L` — semi-longitud del parche (en la dirección $u$)
- La proyección es válida cuando $R \gg L$ (la curvatura del frente de onda es despreciable)

---

### Proyección de la reflectividad

$$p_\theta(u) = \int_{-\infty}^{\infty} g(u\cos\theta - v\sin\theta,\; u\sin\theta + v\cos\theta)\, dv$$

**Variables:**
- `\theta` — ángulo de proyección (ángulo de observación del radar)
- `u` — coordenada a lo largo del eje de proyección (perpendicular a la dirección de propagación)
- `v` — coordenada a lo largo de la dirección de integración (paralela a la propagación)

---

### Teorema de la proyección-rebanada (Projection-Slice Theorem)

$$P_\theta(U) = G(U\cos\theta,\; U\sin\theta)$$

**Variables:**
- `P_\theta(U)` — transformada de Fourier 1D de la proyección $p_\theta(u)$: $P_\theta(U) = \int_{-\infty}^{\infty} p_\theta(u) e^{-jUu}\,du$
- `G(X,Y)` — transformada de Fourier 2D de la reflectividad $g(x,y)$: $G(X,Y) = \int\!\!\int g(x,y) e^{-j(xX+yY)}\,dx\,dy$
- La transformada de la proyección a ángulo $\theta$ es una **rebanada** (slice) de la transformada 2D $G$ a lo largo de la línea a ángulo $\theta$

---

### Señal de retorno procesada (tras demodulación con chirp de referencia)

Multiplicando $r_\theta(t)$ por el chirp de referencia $\cos[\omega_0(t-\tau_0) + \alpha(t-\tau_0)^2]$ y aplicando filtro paso-bajo:

$$\bar{C}_\theta(t) \approx \frac{A}{2} P_\theta\!\left[\frac{2}{c}\!\left(\omega_0 + 2\alpha(t - \tau_0)\right)\right]$$

**Variables:**
- `\tau_0 = 2R/c` — retraso de ida y vuelta al centro del parche
- La señal procesada $\bar{C}_\theta(t)$ es **la transformada de Fourier de la proyección** evaluada en la frecuencia espacial $X = \frac{2}{c}(\omega_0 + 2\alpha(t-\tau_0))$

---

### Frecuencias espaciales accesibles (dominio de Fourier SAR)

Las muestras del dominio de Fourier $G(X,Y)$ adquiridas por el SAR se sitúan en un **segmento anular** (annulus segment) definido por:

$$X_1 = \frac{2}{c}(\omega_0 - \alpha T), \qquad X_2 = \frac{2}{c}(\omega_0 + \alpha T)$$

para $|\theta| \leq \theta_M$ (ángulo máximo de apertura).

**Variables:**
- `X_1, X_2` — frecuencias espaciales mínima y máxima (proporcionales a las frecuencias mínima y máxima del chirp)
- `\theta_M` — semi-ángulo de apertura máximo (rad)
- Las muestras de Fourier disponibles ocupan solo el segmento $[X_1, X_2]$ de cada rebanada radial

---

### Resolución SAR (tomográfica)

$$\delta_x = \frac{2\pi}{\Delta X} = \frac{\pi c}{\alpha T} = \frac{\pi c}{2\alpha T}$$

(resolución en rango, función solo del ancho de banda del chirp)

$$\delta_y = \frac{2\pi}{\Delta Y} \approx \frac{2\pi}{2\omega_0 \sin\theta_M} = \frac{\pi c}{\omega_0 \sin\theta_M}$$

(resolución en acimut, función de la frecuencia central y del ángulo de apertura)

**Variables:**
- `\delta_x` — resolución en la dirección $x$ (acimut / cross-range)
- `\delta_y` — resolución en la dirección $y$ (rango)
- `\alpha T` — producto frecuencia-tiempo del chirp: $\alpha T = B/2$ donde $B$ es el ancho de banda
- `\theta_M` — semi-ángulo de apertura efectivo de la observación SAR
- `\omega_0` — frecuencia portadora (rad/s)

En notación simplificada:

$$\delta_x \approx \frac{\pi c}{2\alpha T}, \qquad \delta_y \approx \frac{\pi c}{2\omega_0 \sin\theta_M}$$

---

### Condición para despreciar el término de fase cuadrática (RVPE)

$$\frac{4\alpha L^2}{c^2} \ll \frac{\pi}{2}$$

equivalentemente:

$$N^2 \ll 2\,\text{TBW}$$

donde $\text{TBW} = T \cdot \frac{2\alpha T}{2\pi} = \frac{\alpha T^2}{\pi}$ es el producto tiempo-ancho de banda del chirp transmitido.

---

### Corrección por efecto Doppler (velocidad de la plataforma)

La velocidad radial $v_r = v\sin\theta$ produce un factor de escala $a \approx 1 + 2v_r/c$ en la tasa de demodulación; el waveform de referencia corregido es:

$$\cos\!\left[\omega_0(at - \tau_0) + \alpha(at - \tau_0)^2\right]$$

Condición para alta resolución en acimut (velocidad de plataforma limitada):

$$\frac{v}{c} < \frac{\lambda_0}{8L}$$

---

### Curvatura del frente de onda — condición de validez de la proyección rectilínea

$$\frac{L^2}{2R} < \delta_x$$

y para coherencia entre distintos ángulos de vista:

$$\frac{L^2 \sin 2\theta_M}{4R} \ll \frac{\lambda}{8}$$

---

## 3. Suposiciones del Modelo

1. El radar opera en **modo spotlight**: la antena se orienta continuamente hacia el mismo parche de terreno durante la colección de datos; la posición del radar cambia pero el parche permanece iluminado.
2. **Imagen tomográfica, no Doppler:** El mecanismo de formación de imagen es de tipo tomográfico (diferentes ángulos de vista), no basado en desplazamientos Doppler; no se requiere movimiento relativo entre plataforma y escena durante la transmisión/recepción individual.
3. **Ángulo de depresión cero** (caso base): el radar vuela a altura $h=0$; para $h \neq 0$ se aplica la corrección de geometría slant-plane mediante estiramiento en la dimensión de rango.
4. **$R \gg L$:** La distancia al parche es mucho mayor que el tamaño del parche, por lo que el frente de onda puede aproximarse como una línea recta (la proyección $p_\theta(u)$ es una integral de línea genuina).
5. **Señal narrow-band:** El SAR es una versión de banda estrecha de la CT: los datos del dominio de transformada se restringen al segmento anular $[X_1, X_2] \times [-\theta_M, \theta_M]$ en el plano polar de Fourier.
6. **El término de fase cuadrática** en la señal demodulada puede despreciarse si el TBW es suficientemente grande respecto al número de celdas de resolución ($N^2 \ll 2\,\text{TBW}$).
7. La **reflectividad** $g(x,y)$ es compleja (amplitud y fase): la amplitud está escalada por $|g|$ y la fase es $\angle g$; la función es independiente de la frecuencia y del ángulo de vista dentro del rango de operación del radar.
8. **Superposición válida** (linealidad): la respuesta del parche completo es la integral de las respuestas individuales de cada elemento diferencial (no hay dispersión múltiple).

---

## 4. Notas Adicionales

- Este paper es un trabajo **seminal** (1983) que estableció la conexión formal entre SAR spotlight y tomografía computarizada (CAT), demostrando que el "polar format Doppler processing" usado en SAR es equivalente a la reconstrucción tomográfica mediante el teorema de proyección-rebanada.
- La principal contribución es mostrar que la señal demodulada en cada ángulo de observación $\theta$ es la **transformada de Fourier de la proyección** de la reflectividad del terreno, lo que permite aplicar algoritmos de CT al SAR.
- El algoritmo de **convolution backprojection** (modificado para SAR) es introducido como alternativa al método directo de Fourier, con mejor calidad de imagen.
- Las diferencias entre SAR y CT clásica son: (a) en SAR la integral de línea es perpendicular a la dirección de propagación de las ondas de radio; (b) SAR es una versión narrow-band de CT; (c) los datos SAR se distribuyen en un segmento anular del plano de Fourier en lugar del disco completo.
- La coherencia (speckle) en SAR se explica naturalmente en este marco: como el SAR es una versión narrow-band de CT, la imagen es la suma coherente de muchos dispersores dentro de cada celda de resolución, con amplitud exponencial y fase uniforme.
- El efecto Doppler y el rango variable durante la apertura imponen limitaciones en el tamaño del parche de terreno que pueden procesarse sin degradación: $v/c < \lambda_0/(8L)$ para el efecto Doppler, y $L^2/(2R) < \delta_x$ para la curvatura del frente de onda.
- Este paper es el antecedente directo del Paper 6 de esta colección (Arikan y Munson, formulación tomográfica del SAR biestático), que extiende la misma metodología al caso biestático.
