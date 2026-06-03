# Autonomous Airborne 3D SAR Imaging System for Subsurface Sensing: UWB-GPR on Board a UAV for Landmine and IED Detection

**Autores:** Maria Garcia-Fernandez, Yuri Alvarez-Lopez, Fernando Las Heras  
**Año:** 2019  
**Fuente/Publicación:** Remote Sensing, MDPI, vol. 11, no. 20, art. 2357. DOI: 10.3390/rs11202357

---

## 1. Geometría del Sistema

El sistema consiste en un UAV multirrotor que vuela autónomamente sobre una región de interés a una altura constante (nominalmente 2.3 m sobre el suelo), realizando pasadas paralelas (along-track en el eje y, across-track en el eje x) sobre una cuadrícula rectangular plana. La geometría es de tipo Down-Looking GPR (DLGPR): las antenas apuntan verticalmente hacia abajo.

- **Antenas:** dos antenas UWB Vivaldi (una transmisora Tx, una receptora Rx), operando de 600 MHz a 3 GHz (procesadas), montadas bajo el UAV.
- **Radar:** módulo UWB M-sequence, rango de trabajo 100 MHz a 6 GHz.
- **Posicionamiento:** sistema RTK dual-banda de alta precisión (exactitud esperada: 0.5 cm horizontal, 1 cm vertical) + laser rangefinder para medir distancia al suelo.
- **Sistema de coordenadas:** se convierte de geodésico (lat/lon/alt) a un sistema local ENU (East–North–Up): ejes $(x_e, y_n, z_u)$.
- El camino de vuelo se rota para alinearlo con los ejes x–y (rotación según el rumbo principal sobre el terreno $\widehat{c_{og}}$).
- El dominio de investigación (donde se calcula la imagen SAR 3D) $\mathbf{r'} = (x', y', z')$ se define dentro del plano de observación, muestreado con pasos $\delta_t = \delta_{ct} = \lambda_{\min}/4 = 0.025$ m en xy y $\delta_z = 0.02$ m en profundidad.
- Geometría de validación experimental: dominio $x' \in [-0.4, 0.6]$ m, $y' \in [1, 5]$ m, $z' \in [-0.6, 0.4]$ m.

---

## 2. Ecuaciones de Resolución SAR/GPR

### Reflectividad SAR (Delay-and-Sum, dominio frecuencia)

$$\rho(\mathbf{r'_p}) = \sum_{n=1}^{N} \rho(\mathbf{r'_p}, n) = \sum_{n=1}^{N} \sum_{m=1}^{M} E_{scatt}(f_n, \mathbf{r}^m_t, \mathbf{r}^m_r) \exp(+j k_{0,n} R_{p,m})$$

**Variables:**
- `r'_p` — punto del dominio de investigación donde se calcula la reflectividad
- `f_n` — n-ésima frecuencia de interés (de N frecuencias totales)
- `r^m_t`, `r^m_r` — posiciones de la antena transmisora y receptora en el m-ésimo punto de medida
- `E_scatt(f_n, r^m_t, r^m_r)` — campo dispersado medido (en frecuencia) en el m-ésimo punto de adquisición
- `k_{0,n}` — número de onda en espacio libre para la n-ésima frecuencia
- `R_{p,m}` — longitud total del camino de propagación entre antena transmisora, píxel, y antena receptora
- `M` — número total de medidas (puntos de adquisición)
- `N` — número de frecuencias

### Longitud de camino en espacio libre

$$R_{p,m} = \|\mathbf{r}^m_t - \mathbf{r'_p}\| + \|\mathbf{r}^m_r - \mathbf{r'_p}\|$$

**Variables:**
- `r^m_t` — posición de la antena transmisora en la m-ésima medida
- `r^m_r` — posición de la antena receptora en la m-ésima medida
- `r'_p` — posición del píxel de reconstrucción

### Longitud de camino considerando permitividad del suelo (refracción)

Para píxeles bajo la superficie del suelo ($z'_p < 0$), se aplica la corrección de Snell con $n_s = \sqrt{\varepsilon_r - 1} - \sqrt{\varepsilon_r}$:

$$R_{p,m} = 2d\sqrt{\varepsilon_r - 1} + \frac{d_t(d_t - d_n \cos(2\phi_t))}{d_t + d_n \sin(2\phi_t)^2} + \frac{d_r(d_r - d_n \cos(2\phi_r))}{d_r + d_n \sin(2\phi_r)^2}$$

**Variables:**
- `d` — profundidad del punto bajo la superficie del suelo
- `epsilon_r` — permitividad relativa del suelo
- `d_t`, `d_r` — distancias horizontales desde la antena transmisora/receptora hasta el punto de refracción en la interfaz
- `phi_t`, `phi_r` — ángulos de incidencia en la interfaz aire-suelo del rayo transmisor/receptor
- `d_n` — parámetro geométrico auxiliar

### SAR con Ecualización de la respuesta en frecuencia

$$\rho(\mathbf{r'_p}) = \sum_{n=1}^{N} \frac{\rho(\mathbf{r'_p}, n)}{\max\{|\rho(\mathbf{r'_p}, n)|\}}$$

**Variables:**
- `rho(r'_p, n)` — imagen SAR calculada para la n-ésima frecuencia individual
- El denominador normaliza cada imagen por su valor máximo en valor absoluto antes de sumar, igualando el peso de todas las frecuencias

### Estimación del rumbo principal sobre el terreno

$$c_{og}[rad] = \text{mod}\left(\text{atan}\left(\frac{v_n}{v_e}\right), \pi\right) \approx \text{mod}\left(\text{atan}\left(\frac{\partial y_n}{\partial x_e}\right), \pi\right)$$

**Variables:**
- `c_og` — course over ground (rumbo sobre el terreno), en radianes
- `v_n`, `v_e` — velocidades del UAV en dirección norte y este, respectivamente
- `x_e`, `y_n` — coordenadas en el sistema ENU (este y norte)

### Tiempo-cero y corrección de rango

$$t_0 = \frac{2(d_{eg} - d_{rg})}{c}, \quad d_{eg} = \frac{c \cdot t_g}{2}, \quad r_{ng} = \frac{c \cdot t}{2}$$

**Variables:**
- `t_0` — corrección de tiempo-cero
- `d_eg` — distancia estimada por radar al suelo
- `d_rg` — distancia real al suelo (del laser rangefinder)
- `c` — velocidad de la luz
- `t_g` — instante de tiempo donde se detecta el suelo
- `r_ng` — distancia de rango en el eje corregido

### Detección aparente de objeto enterrado (efecto dieléctrico)

Si un objeto está enterrado a profundidad real $d_{obj}$ y la permitividad del suelo es $\varepsilon_r$, la imagen SAR (sin corregir) lo detecta a:

$$d_{aparente} = \sqrt{\varepsilon_r} \cdot d_{obj}$$

**Variables:**
- `d_obj` — profundidad real del objeto bajo la superficie
- `epsilon_r` — permitividad relativa del suelo
- `d_aparente` — profundidad donde aparece el objeto en la imagen SAR sin corrección dieléctrica

---

## 3. Suposiciones del Modelo

1. El UAV vuela a altura aproximadamente constante sobre una superficie plana (vuelo autónomo a 2.3 m).
2. Las medidas se geo-referencian con precisión sub-centimétrica mediante RTK dual-banda; la precisión requerida es $< \lambda_{\min}/4 = 2.5$ cm en el plano horizontal y $< \lambda_{\min}/8 = 1.25$ cm en la dirección vertical (rango).
3. Se asume propagación en espacio libre para el procesado básico; el efecto del suelo se introduce como corrección opcional mediante la permitividad relativa $\varepsilon_r$.
4. El radar transmite una secuencia binaria de máxima longitud (MLBS); la respuesta impulsional se obtiene por correlación cruzada con la secuencia ideal.
5. El fondo (clutter) se estima como el promedio de todas las medidas y se resta (background subtraction).
6. La corrección de altura consiste en desplazar cada medida en rango por $(z - \bar{z})$ para simular adquisición a altura constante, mejorando el enfoque SAR.
7. La banda de procesado usada en validación experimental es 600 MHz a 3 GHz (se descarta la parte de alta frecuencia por excesiva atenuación en suelo).
8. Se usa un único par Tx-Rx (configuración monoestática/cuasi-monoestática); las posiciones de Tx y Rx se calculan a partir de las coordenadas RTK más ángulos de actitud (roll, pitch, yaw) y distancias geométricas internas del sistema.

---

## 4. Notas Adicionales

- **Validación experimental:** dos medidas en campo abierto (aeródromo de la Universidad de Oviedo): (1) caja metálica cilíndrica de 9.5 cm de radio descubierta en un hoyo de 8 cm de profundidad; (2) la misma caja cubierta con tierra (enterrada a 8 cm). Permitividad del suelo estimada en $\varepsilon_r = 3$, lo que desplaza la detección de $d_{box} = 0.08$ m a $d_{aparente} \approx 0.14$ m en la imagen sin corrección.
- **Mejoras propuestas:** (a) corrección de altura (height correction) antes del procesado SAR, lo que reduce la reflexión de la interfaz aire-suelo y mejora la relación señal-clutter; (b) ecualización de la respuesta en frecuencia para aprovechar todo el ancho de banda UWB y mejorar la resolución en rango (eje z).
- **Resolución teórica en profundidad:** determinada por el ancho de banda: $\Delta R = c/(2\Delta BW)$; con BW = 2.4 GHz y $\varepsilon_r = 1$ se obtiene $\Delta R \approx 6.25$ cm en aire.
- **Selección de datos:** se descartan medidas con variaciones bruscas de actitud, velocidad nula (sobrelapso en los giros), y desviaciones del rumbo principal superiores a $th_c = 20°$. Umbral de mínimo desplazamiento: $th_\Delta = \lambda_{\min}/6 \approx 0.017$ m.
- **Patente relacionada:** WO/2017/125627, PCT/ES2017/000006, Universidad de Oviedo / Universidad de Vigo.
- La imagen SAR final es un volumen 3D $\rho(x', y', z')$ que puede visualizarse en cortes YZ (along-track vs. profundidad), YX (vista superior) y XZ (across-track vs. profundidad).
