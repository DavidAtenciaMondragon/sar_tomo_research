# A Fast Back-Projection Algorithm for Bistatic SAR Imaging

**Autores:** Yu Ding, David C. Munson Jr.  
**Año:** 2002  
**Fuente/Publicación:** IEEE International Conference on Image Processing (ICIP 2002), Vol. II, pp. 449–452. University of Illinois at Urbana-Champaign.

---

## 1. Geometría del Sistema

**Geometría biestática** en modo spotlight: el transmisor y el receptor se desplazan a lo largo de trayectorias **rectilíneas distintas y paralelas**. La dirección de rango cruzado $x$ es paralela a la trayectoria del receptor; la dirección de rango $y$ es perpendicular a $x$. Ambas antenas apuntan al mismo parche de terreno durante toda la operación (modo spotlight).

- Transmisor en posición $(x_t, y_t)$ → trayectoria rectilínea
- Receptor en posición $(x_r, y_c)$ → trayectoria rectilínea a distancia $y_c$ del centro de la escena
- Ángulo bistático $\beta$: ángulo entre la línea que une transmisor-receptor y la bisectriz
- La distancia más corta del receptor al centro de la escena es $y_c$

---

## 2. Ecuaciones de Resolución SAR

### Distancia de ida y vuelta (round-trip range)

$$R(x_t, y_t, x_r; x', y') = \sqrt{(x_t - x')^2 + (y_t - y')^2} + \sqrt{(x_r - x')^2 + (y_c - y')^2}$$

**Variables:**
- `(x', y')` — posición del reflector en la escena
- `(x_t, y_t)` — posición del transmisor
- `(x_r, y_c)` — posición del receptor
- `R_0` — distancia del transmisor al centro de la escena y de vuelta al receptor

### Pulso transmitido (chirp LFM)

$$s(\tau) = p(\tau) \cdot e^{j(\omega_0 \tau - \pi\alpha\tau^2)}$$

con $p(\tau) = 1$ si $|\tau| \leq T/2$, $0$ en caso contrario.

**Variables:**
- `\omega_0` — frecuencia angular central
- `\alpha` — tasa de barrido del chirp (Hz/s)
- `T` — duración del pulso

### Señal recibida (retorno de un reflector puntual)

$$h(x_t, y_t, x_r; \tau; x', y') = \cos\!\left(\omega_0(\tau - \Delta_\tau) - \pi\alpha(\tau - \Delta_\tau)^2\right)$$

donde $\Delta_\tau = R(x_t, y_t, x_r; x', y')/c$.

### Señal SAR bistática demodulada

$$d(x_r, \omega) = \iint g(x', y') \exp\!\left(-j\frac{\omega}{c} R(x_r; x', y')\right) dx'\, dy'$$

donde

$$\omega = \omega_0 - 2\pi\alpha\left(\tau - \frac{R_0(x_r; x', y')}{c}\right)$$

**Variables:**
- `g(x',y')` — reflectividad de la escena en $(x', y')$
- `\omega` — frecuencia angular instantánea (tras demodulación)
- `R_0(x_r; x', y')` — distancia de referencia (transmisor→centro→receptor)

### Señal antes del backprojection

$$b(x_r, l) = \int_\omega d(x_r, \omega) \exp(-j\omega l) \, d\omega$$

**Variables:**
- `l` — variable de tiempo transformada (tiempo de viaje)

### Backprojection bistático (estándar)

$$\hat{g}(x, y) = \sum_{i=1}^{P} b\!\left(x_{ri},\, \frac{\sqrt{(x_{ti}-x)^2+(y_{ti}-y)^2} + \sqrt{(x_{ri}-x)^2+(y_c-y)^2}}{c}\right) \Delta x_{ri}$$

**Variables:**
- `P` — número total de proyecciones (posiciones del receptor)
- `(x_{ti}, y_{ti})` — posición del transmisor en la proyección $i$
- `(x_{ri})` — posición del receptor en la proyección $i$
- `\Delta x_{ri}` — espaciado entre posiciones del receptor

### Complejidad: BP bistático estándar

$$\mathcal{O}(N^3)$$

para una imagen de $N \times N$ con $N$ proyecciones.

### Complejidad: BP bistático rápido

$$\mathcal{O}(N^2 \log_2 N)$$

Reducida mediante descomposición recursiva de la imagen en 4 subimágenes, explotando la propiedad angular band-limited del radar bistático.

### Factor de escala bistático

La señal SAR bistática es la TF-2D de la reflectividad $g(x',y')$, distribuida a lo largo de elipses en el dominio de Fourier. El ancho de banda bistático efectivo en la dirección perpendicular a la bisectriz se escala por:

$$\frac{2}{c}\cos\!\left(\frac{\beta}{2}\right)$$

donde $\beta$ es el ángulo bistático, lo que limita la resolución respecto al caso monoestático.

---

## 3. Suposiciones del Modelo

1. Ambas antenas se desplazan en modo spotlight, apuntando siempre al mismo parche de terreno.
2. Las trayectorias del transmisor y del receptor son **rectilíneas** (vuelo en línea recta).
3. Aproximación far-field bistática: la variación del ángulo bistático durante la colección es pequeña (pocos grados en la simulación, β ≈ 90°).
4. El transmisor envía chirps LFM con ancho de banda plano.
5. La reflectividad de la escena es estacionaria durante la colección de datos.
6. Modelo de propagación en espacio libre (sin efectos atmosféricos ni de terreno).
7. Las trayectorias del transmisor y receptor son paralelas a la dirección de rango cruzado $x$.
8. La propiedad angular band-limited (usada en el Fast BP) es válida cuando la imagen es menor que la mitad del tamaño de las proyecciones.

---

## 4. Notas Adicionales

- El algoritmo Fast BP bistático es una extensión directa del Fast BP para tomografía (Basu & Bresler, 2000) al caso SAR bistático.
- Diferencias clave con el caso monoestático: (1) el BP se hace sobre trayectorias elípticas en lugar de circulares; (2) los datos SAR están en una región de cobertura angular limitada (~pocos grados); (3) las interpolaciones son en banda de paso, no paso bajo.
- Resultados de simulación: imagen de $256 \times 256$ con 4 blancos puntuales. Tiempo CPU BP estándar: 402 s; Fast BP: 59 s. Speedup ~7×.
- Para imágenes mayores de $256 \times 256$ el speedup aumenta según $N/\log N$.
- Parámetros de simulación: portadora 10 GHz, ancho de banda 500 MHz, ángulo de visión del receptor 3°, ángulo bistático ≈ 90°.
