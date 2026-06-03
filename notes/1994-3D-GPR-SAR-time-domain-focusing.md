# Three-Dimensional Ground-Penetrating Radar Imaging Using Synthetic Aperture Time-Domain Focusing

**Autores:** Erik M. Johansson, Jeffrey E. Mast  
**Año:** 1994  
**Fuente/Publicación:** SPIE Vol. 2275, Proceedings of the International Symposium on Optics, Imaging, and Instrumentation (Advanced Microwave and Millimeter-Wave Detectors), pp. 205–214. Lawrence Livermore National Laboratory.

---

## 1. Geometría del Sistema

El sistema es un GPR monoestático (o biestático/multiestático) operando en modo pulse-echo con radar de impulsos UWB (ultra-wideband). Una antena (o array de antenas) se mueve sobre una apertura sintética 2D plana (cuadrícula rectangular en los planos x e y), formando una apertura sintética que ilumina un volumen de objetos enterrados en el semieespacio inferior (eje z apuntando hacia abajo, positivo en profundidad).

- **Modo monoestático:** una sola antena opera como transmisor y receptor simultáneamente, desplazada a N posiciones $(x_n, y_n, z_n)$ en la apertura.
- **Modo biestático:** transmisor fijo en $(x_t, y_t, z_t)$, N receptores en $(x_{r_n}, y_{r_n}, z_{r_n})$.
- **Modo multiestático:** M transmisores y N receptores, obteniendo $M \times N$ formas de onda $R_{mn}(t)$.
- La apertura puede ser plana y uniformemente muestreada, pero el algoritmo es genérico para posiciones no planares y no uniformes.
- **Aplicación de validación:** losa de hormigón con barras de refuerzo metálicas (rebars). Apertura de 1.5 m × 1.5 m con paso de muestreo 1.27 cm. Antena de bocina con ancho de haz aproximado de 60°. Separación entre antena y superficie del hormigón: ~7.5 cm. Señal: pulso de baja potencia (~1 W pico), frecuencias de 1.2 a 3.5 GHz.
- **Medios multicapa:** la geometría considera hasta P capas planares con permitividades $\varepsilon_1, \varepsilon_2, \ldots, \varepsilon_P$ y límites a profundidades conocidas.

---

## 2. Ecuaciones de Resolución SAR/GPR

### SAFT monoestático — estimación de la distribución de objetos (dominio tiempo)

$$\hat{o}(x, y, z) = \frac{1}{N} \sum_{n} R_n\!\left(\frac{2r_n}{v}\right)$$

**Variables:**
- `o_hat(x,y,z)` — estimación de la distribución de dispersores en el punto $(x,y,z)$
- `N` — número de posiciones de antena (apertura sintética)
- `R_n(t)` — forma de onda recibida por la n-ésima antena (en el dominio del tiempo)
- `r_n` — distancia de la n-ésima antena al punto $(x,y,z)$ en el volumen de reconstrucción
- `v` — velocidad de propagación del pulso en el medio
- `2r_n/v` — tiempo de viaje de ida y vuelta (round-trip travel time)

### Distancia antena–punto de reconstrucción (monoestático)

$$r_n = \sqrt{(x - x_n)^2 + (y - y_n)^2 + (z - z_n)^2}$$

**Variables:**
- `(x_n, y_n, z_n)` — posición de la n-ésima antena en la apertura
- `(x, y, z)` — punto del volumen de reconstrucción

### SAFT biestático — estimación de la distribución de objetos

$$\hat{o}(x, y, z) = \frac{1}{N} \sum_{n} R_n\!\left(\frac{r_t + r_{r_n}}{v}\right)$$

**Variables:**
- `r_t` — distancia del transmisor (fijo) al punto de reconstrucción $(x,y,z)$
- `r_{r_n}` — distancia del n-ésimo receptor al punto de reconstrucción
- `R_n(t)` — forma de onda recibida por el n-ésimo receptor

### Distancia transmisor–punto (biestático)

$$r_t = \sqrt{(x - x_t)^2 + (y - y_t)^2 + (z - z_t)^2}$$

### Distancia receptor–punto (biestático)

$$r_{r_n} = \sqrt{(x - x_{r_n})^2 + (y - y_{r_n})^2 + (z - z_{r_n})^2}$$

### SAFT multiestático (M transmisores, N receptores)

$$\hat{o}(x, y, z) = \frac{1}{MN} \sum_{m} \sum_{n} R_{mn}\!\left(\frac{r_{t_m} + r_{r_n}}{v}\right)$$

**Variables:**
- `R_{mn}(t)` — forma de onda recibida usando el m-ésimo transmisor y el n-ésimo receptor
- `r_{t_m}` — distancia del m-ésimo transmisor al punto $(x,y,z)$

### Punto de inflexión en interfaz bicapa (método de aproximación)

$$x_b = x_2 + \sqrt{\frac{\varepsilon_{r1}}{\varepsilon_{r2}}}(x_1 - x_2)$$

**Variables:**
- `x_b` — punto de inflexión aproximado en la interfaz entre medios (donde el rayo se refracta)
- `x_1` — proyección lateral de la antena sobre la interfaz (punto de refracción si $\varepsilon_{r2} = \varepsilon_{r1}$, i.e., rayo recto)
- `x_2` — proyección lateral del punto de reconstrucción sobre la interfaz (punto de refracción si $\varepsilon_{r2} \gg \varepsilon_{r1}$, rayo casi vertical)
- `epsilon_{r1}` — permitividad relativa del primer medio (el que contiene la antena)
- `epsilon_{r2}` — permitividad relativa del segundo medio (el que contiene el objeto)
- La expresión usa $\sqrt{\varepsilon_{r1}/\varepsilon_{r2}}$ porque este término representa la razón de los índices de refracción

### Compensación de pérdidas por trayecto (path loss)

Para un sistema monoestático en medio homogéneo, cada punto del tiempo en la forma de onda se multiplica por $r^2$, donde $r = vt$ ($v$ velocidad, $t$ índice temporal):

$$R_n^{comp}(t) = R_n(t) \cdot r^2 = R_n(t) \cdot (vt)^2$$

Para geometría biestática se multiplica por $r_t \cdot r_r$ (producto de las distancias transmisor–objeto y objeto–receptor). En medios multicapa, $r_t$ y $r_r$ se calculan sumando las distancias recorridas dentro de cada medio.

---

## 3. Suposiciones del Modelo

1. El medio de propagación es (inicialmente) homogéneo; la extensión a múltiples capas requiere conocer las profundidades de las interfaces y las permitividades de cada capa.
2. El objetivo se modela como un conjunto de dispersores puntuales (point scatterers); el campo retrodispersado de cada punto es una versión escalada y retrasada del pulso transmitido.
3. El modo de adquisición puede ser tiempo-dominio directo (ADC rápido) o stepped-frequency (con IFFT para sintetizar el dominio temporal).
4. La deconvolución del pulso transmitido de la señal recibida mejora la relación señal-ruido y la resolución en rango; el pulso transmitido ideal sería un impulso de Dirac.
5. El "coupling pulse" (acoplamiento directo Tx-Rx) se elimina por sustracción de fondo (background subtraction, promedio sobre todas las posiciones) en el caso monoestático; en el biestático requiere técnicas alternativas (enmascaramiento).
6. La corrección del diagrama de antena (beam pattern correction) escala cada punto de la forma de onda por la inversa de la amplitud del diagrama en la dirección del punto de reconstrucción; también permite enmascarar datos fuera del lóbulo principal para suprimir lóbulos de difracción (grating lobes que aparecen como arcos o "cat whiskers").
7. Para medios multicapa, el tiempo de viaje correcto requiere encontrar el punto de inflexión en cada interfaz (problema no cerrado en forma analítica); se resuelve por el método de aproximación (Ec. 7) o por minimización del tiempo de viaje (método steepest descent).
8. La señal analítica (señal compleja via transformada de Hilbert) se usa en el procesado para tener un único pico por celda de resolución y poder analizar fase y envolvente por separado.

---

## 4. Notas Adicionales

- **Nombre del algoritmo:** SAFT (Synthetic Aperture Focusing Technique); equivalente al algoritmo Delay-and-Sum (DAS) en dominio tiempo con corrección de trayecto y diagrama de antena.
- **Comparación con tomografía por difracción (frecuencia):** el algoritmo tiempo-dominio es simple y eficiente para problemas pequeños, adaptable a geometrías no estándar (no planares, biestáticas, multiestáticas), pero ineficiente para problemas grandes y difícil de adaptar a parámetros dependientes de la frecuencia (dispersión, atenuación). El algoritmo frecuencia-dominio es muy eficiente para geometrías monoestáticas con cuadrícula plana y uniforme, y fácilmente adaptable a múltiples capas, pero requiere grilla plana e igual espaciado.
- **Resultado experimental:** imagen volumétrica 3D de losa de hormigón de 1.5 × 1.5 m con apertura de 1.27 cm de paso. Dos barras de refuerzo claramente visibles en la imagen; parte de una tercera barra ortogonal y una delaminación también detectables en cortes planares.
- **Velocidad de propagación en hormigón:** implícita en la corrección de la distancia $r_n$; para hormigón típico $\varepsilon_r \approx 4-9$, lo que da $v = c/\sqrt{\varepsilon_r} \approx 0.1$–0.15 m/ns.
- **Señal analítica:** se obtiene tomando la FFT de la señal real, enmascarando las componentes de frecuencia negativa y cero, y aplicando la IFFT inversa. La envolvente (módulo) tiene un único pico por celda de resolución, facilitando la interpretación.
- Este paper es considerado una de las referencias fundacionales del SAFT 3D para GPR y es citado como referencia [22] en el paper de Garcia-Fernandez 2019.
