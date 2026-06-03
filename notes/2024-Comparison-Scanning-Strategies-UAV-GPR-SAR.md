# Comparison of Scanning Strategies in UAV-Mounted Multichannel GPR-SAR Systems Using Antenna Arrays

**Autores:** María García-Fernández, Guillermo Álvarez-Narciandi, Fernando Las Heras, Yuri Álvarez-López  
**Año:** 2024  
**Fuente/Publicación:** IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, vol. 17, pp. 3571–3586. DOI: 10.1109/JSTARS.2024.3351602

---

## 1. Geometría del Sistema

El sistema es un GPR-SAR multicanal montado en UAV, con un array de antenas compuesto por 3 transmisores (TX) y 4 receptores (RX) → 12 canales TX-RX de adquisición. El UAV vuela autónomamente sobre el área inspeccionada a altura nominal $h \approx 1.5$ m, realizando barridos paralelos (along-track, eje y) separados una distancia across-track (eje x).

- **Array de antenas:** 2 subarrays Vivaldi (TX: 3 antenas + switch; RX: 4 antenas + switch). Espaciado físico entre elementos del array: 13.33 cm. Banda de trabajo: 1–6 GHz (procesado: 1–3 GHz). Adquisición completa de 12 canales en 280 ms.
- **Posicionamiento:** RTK multibanda + laser rangefinder (precisión cm-level).
- **Dominio de observación:** área de $4.5 \times 12$ m en la validación en vuelo (campo de tiro militar "El Palancar", norte de Madrid). 15 objetivos enterrados de distintos materiales, tamaños y profundidades (hasta 14 cm).
- **Permitividad del suelo estimada:** $\varepsilon_r \approx 4$.
- **Dos estrategias de escaneo comparadas:**
  - Esquema uniforme $U$-$\Delta x$: barridos along-track separados uniformemente una distancia $\Delta x$ en across-track.
  - Esquema no uniforme $3X$-$\Delta x$: el array realiza tres barridos separados $\lambda_{\min}/2$ entre sí, luego se desplaza $\Delta x - \lambda_{\min}$ en across-track y repite (esquema en zigzag).
- **Cuatro configuraciones evaluadas:** $U$-40 cm, $U$-20 cm, $3X$-80 cm, $3X$-40 cm.

---

## 2. Ecuaciones de Resolución SAR/GPR

### Resolución espacial SAR (along-track y across-track)

$$\delta_u \approx \lambda_c \frac{\sqrt{L^2/4 + h^2}}{2L}$$

**Variables:**
- `delta_u` — resolución espacial en along-track (y) o across-track (x) [m]
- `lambda_c` — longitud de onda a la frecuencia central $f_c$ [m]
- `L` — longitud de la apertura SAR (longitud de la máscara) [m]
- `h` — altura de vuelo del UAV [m]

Para los parámetros de este sistema ($h \approx 1.5$ m, $L = 2$ m, $\lambda_c = 0.15$ m a $f_c = 2$ GHz):
$$\delta_{\text{along-track}} \approx 0.45\lambda_c \approx 6.8 \text{ cm}, \qquad \delta_{\text{across-track}} \approx 0.79\lambda_c \approx 11.9 \text{ cm}$$

### Resolución en profundidad (eje z)

$$\delta_z \approx \frac{v}{2\,\text{BW}} = \frac{c}{2\sqrt{\varepsilon_r}\,\text{BW}} \approx \frac{7.5}{\sqrt{\varepsilon_r}} \text{ cm}$$

**Variables:**
- `v` — velocidad de propagación en el suelo [m/s]
- `BW` — ancho de banda del sistema [Hz] (aquí BW = 2 GHz, de 1 a 3 GHz)
- `c` — velocidad de la luz (0.3 m/ns)
- `epsilon_r` — permitividad relativa del suelo

### Criterio de muestreo espacial (intervalo de Nyquist en apertura)

$$\Delta_u \approx \lambda_{\min} \frac{\sqrt{L^2/4 + h^2}}{2L}$$

**Variables:**
- `Delta_u` — intervalo de muestreo espacial requerido en along-track y across-track [m]
- `lambda_min` — longitud de onda mínima (máxima frecuencia de procesado) [m]

Con $\lambda_{\min} = 0.10$ m ($f_{\max} = 3$ GHz), $h = 1.5$ m, $L = 2$ m:
$$\Delta_{\text{along-track}} \approx 0.45\lambda_{\min} \approx 4.5 \text{ cm}, \qquad \Delta_{\text{across-track}} \approx 0.79\lambda_{\min} \approx 7.9 \text{ cm}$$

En la práctica el criterio se fija en $\lambda_{\min}/2 = 5$ cm.

### Peak Signal-to-Clutter Ratio (PSCR)

$$\text{PSCR [dB]} = 10\log_{10}\left(\frac{\max_{x,y \in A_t} |\rho(x, y, z_t)|^2}{\frac{1}{N_c}\sum_{x,y \in A_c} |\rho(x, y, z_t)|^2}\right)$$

**Variables:**
- `A_t` — región del objetivo (rectángulo centrado en la posición del target que encierra su extensión física)
- `A_c` — región de clutter (área de 1 m × 1 m centrada en la posición del objetivo, excluyendo $A_t$)
- `rho(x, y, z_t)` — reflectividad SAR en el corte horizontal a la profundidad $z_t$ del objetivo
- `N_c` — número de píxeles en la región de clutter $A_c$

### Masked SAR — algoritmo DAS con máscara espacial

El algoritmo Masked SAR restringe las contribuciones al cálculo de cada píxel a las medidas tomadas dentro de una región llamada máscara, centrada en el punto de reconstrucción:

- Tamaño de la máscara: $L_{\text{along-track}} = 2$ m (equivalente al ancho de haz proyectado de la antena) × $L_{\text{across-track}} = 1$ m.
- Esto reduce el clutter de reflexiones fuera del lóbulo principal de la antena y disminuye el coste computacional.

### Número de celdas del dominio de imagen (particionado, $3X$ scheme)

Para el esquema $3X$, el dominio de imagen se reconstruye en celdas cuyo tamaño en across-track es $\Delta x$ (el mismo que el paso de muestreo de adquisición) y cuyo tamaño en along-track es el intervalo de muestreo along-track. La cobertura se evalúa como la fracción de celdas que contienen al menos una medida.

---

## 3. Suposiciones del Modelo

1. El UAV vuela a altura aproximadamente constante $h \approx 1.5$ m; las desviaciones se corrigen con el sistema de posicionamiento RTK + laser rangefinder.
2. El medio se trata como homogéneo con una permitividad relativa uniforme $\varepsilon_r$ (estimada de las medidas). Para los objetivos plásticos se asume $\varepsilon_r \approx 3$; para el objetivo de madera (PP, tarima) $\varepsilon_r \approx 2$.
3. El eje de simetría del array está orientado en la dirección across-track (eje x); los barridos along-track se realizan en la dirección y.
4. El espaciado físico del array (13.33 cm) no cumple el criterio de Nyquist para $f_{\max} = 3$ GHz ($\lambda_{\min}/2 = 5$ cm), lo que genera lóbulos de difracción (grating lobes). El esquema $3X$ está diseñado para minimizar este efecto añadiendo virtualmente elementos intermedios.
5. La técnica de co-registro (co-registration) consiste en registrar las imágenes de los barridos forward y backward del mismo canal según su intensidad para mejorar el enfoque; se aplica opcionalmente en el paso de postprocesado.
6. La detección automática de objetivos se realiza con un detector CFAR (Constant False Alarm Rate) aplicado a cada corte horizontal de la imagen 3D GPR-SAR.
7. El tiempo de adquisición completa de 12 canales es 280 ms; a velocidad de 50 cm/s, esto equivale a una medida cada 14 cm en along-track.

---

## 4. Notas Adicionales

- **Comparación de rendimiento (vuelo sobre área 4.5 × 12 m, Tabla IV):**
  - $U$-20 cm: PSCR medio = 16.7 dB, tiempo de vuelo = 10 min (1 vuelo, 24 barridos).
  - $3X$-40 cm: PSCR medio = 17.7 dB, tiempo de vuelo = 15 min (2 vuelos, 36 barridos).
  - $U$-40 cm (disperso): PSCR medio = 14.6 dB, tiempo = 5 min.
  - $3X$-80 cm (disperso): PSCR medio = 15.3 dB, tiempo = 7.5 min.
- **Ventaja del array multicanal respecto a un único canal (1TX-1RX):** para inspeccionar 4.5 × 12 m, el esquema $3X$-40 cm solo requiere el 39.6% de los barridos que necesitaría un sistema DLGPR convencional de un solo canal; el esquema $U$-20 cm requiere solo el 26.4%.
- **Cobertura en vuelo real:** $U$-20 cm → 83.1% (24 barridos); $3X$-40 cm → 92.7% (36 barridos). La cobertura del esquema $3X$ es mayor porque recopila más muestras por celda.
- **Targets detectados:** todos los objetivos excepto el target (xii) (pequeña mina AP de 8.5 cm) fueron detectados visualmente con ambos esquemas densos. El CFAR detectó todos excepto (x) (botella de agua pequeña) y (xii).
- **Esquemas dispersos** ($U$-40 y $3X$-80): adecuados para detección de objetivos medianos/grandes; la mitad del tiempo de vuelo. No detectan el objetivo AP más pequeño (vi).
- **Inspección de 60 m² en 10 min** con el esquema $U$-20 cm es el resultado clave de aplicabilidad operacional.
- La implementación del esquema $3X$ requiere software de control de misión que defina una trayectoria de vuelo no uniforme (más complejo que el esquema uniforme).
