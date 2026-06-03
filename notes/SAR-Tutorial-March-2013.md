# A Tutorial on Synthetic Aperture Radar

**Autores:** Alberto Moreira, Pau Prats-Iraola, Marwan Younis, Gerhard Krieger, Irena Hajnsek, Konstantinos P. Papathanassiou  
**Año:** 2013  
**Fuente/Publicación:** IEEE Geoscience and Remote Sensing Magazine, Vol. 1, No. 1, March 2013. DOI: 10.1109/MGRS.2013.2248301

---

## 1. Geometría del Sistema

El sistema SAR es un radar de apertura sintética montado sobre una plataforma en movimiento (aeronave o satélite). La geometría es **monoestática de visión lateral** (side-looking), donde el radar transmite pulsos electromagnéticos y recibe los ecos de retrodispersión desde la superficie terrestre.

- **Dimensión en rango (slant range / fast time):** dirección perpendicular a la trayectoria de vuelo, medida como tiempo de ida y vuelta de la señal.
- **Dimensión en acimut (azimuth / slow time):** dirección a lo largo de la trayectoria de vuelo.
- La plataforma se mueve a velocidad $v$ con altura $H$. La distancia mínima al objetivo es $r_0 = \sqrt{(H - \Delta h)^2 + x_0^2}$.

**Modos de operación cubiertos:**
- **Stripmap:** el haz de antena apunta en dirección fija, iluminando una franja continua. Resolución en acimut = $d_a/2$.
- **ScanSAR:** el haz se dirige sucesivamente a diferentes sub-franjas (sub-swaths) para aumentar el ancho de faja, sacrificando resolución en acimut.
- **Spotlight:** el haz se orienta continuamente hacia un parche fijo del terreno, aumentando el tiempo de iluminación y por tanto la resolución en acimut a costa de no generar una imagen continua.

También se describen técnicas avanzadas: **polarimetría**, **interferometría**, **interferometría diferencial (DInSAR)**, **tomografía SAR**, **beamforming digital**, **MIMO** y configuraciones **bi/multi-estáticas** (p. ej. TanDEM-X).

---

## 2. Ecuaciones de Resolución SAR

### Resolución en rango (slant range)

$$\delta_r = \frac{c_0}{2 B_r}$$

**Variables:**
- `c_0` — velocidad de la luz (m/s)
- `B_r` — ancho de banda del chirp transmitido (Hz); $B_r = k_r \cdot \tau$ donde $k_r$ es la tasa de chirp y $\tau$ es la duración del pulso

---

### Resolución en acimut (SAR enfocado)

$$\delta_a = \frac{d_a}{2}$$

**Variables:**
- `d_a` — longitud física de la antena en la dirección de acimut (m)
- El factor 2 proviene del procesado coherente de la apertura sintética (longitud $L_{sa} = \Theta_a \cdot r_0 = \frac{\lambda}{d_a} r_0$)

---

### Resolución en acimut del sistema SLAR (sin apertura sintética)

$$\delta_a^{SLAR} = \Theta_a \cdot r_0 = \frac{\lambda}{d_a} r_0$$

**Variables:**
- `\Theta_a` — ancho del haz en acimut (rad) $= \lambda/d_a$
- `r_0` — distancia de aproximación mínima (m)
- `\lambda` — longitud de onda (m)

---

### Distancia instantánea radar–objetivo (historia de rango)

$$r(t) = \sqrt{r_0^2 + (vt)^2} \approx r_0 + \frac{(vt)^2}{2r_0} \quad \text{para } vt/r_0 \ll 1$$

**Variables:**
- `t` — tiempo en acimut (slow time), con $t=0$ en el punto de máximo acercamiento
- `v` — velocidad de la plataforma (m/s)
- `r_0` — distancia mínima al objetivo (m)

---

### Señal de acimut (modelo del eco)

$$s_a(t) = A\sqrt{\sigma_0} \exp(i\varphi^{scat}) \exp\!\left(-i\frac{4\pi}{\lambda} r(t)\right)$$

**Variables:**
- `A` — factor de amplitud (incluye potencia transmitida, pérdidas, patrón de antena)
- `\sigma_0` — sección transversal de retrodispersión normalizada por área
- `\varphi^{scat}` — fase de dispersión del objetivo
- `\lambda` — longitud de onda (m)
- `r(t)` — distancia instantánea radar–objetivo (m)

---

### Frecuencia Doppler instantánea

$$f_D = -\frac{1}{2\pi}\frac{\partial}{\partial t}\frac{4\pi r(t)}{\lambda} = -\frac{2v^2 t}{\lambda r_0}$$

**Variables:**
- `f_D` — frecuencia Doppler instantánea (Hz)
- `v` — velocidad de la plataforma (m/s)

---

### Migración de celda de rango (Range Cell Migration, RCM)

$$RCM(t) = \sqrt{r_0^2 + (vt)^2} - r_0 \approx \frac{(vt)^2}{2r_0}$$

**Variables:**
- Describe el desplazamiento del objetivo a través de celdas de rango durante la iluminación

---

### Tasa de muestreo mínima en acimut (PRF)

$$PRF \geq B_D = \frac{2v}{d_a}$$

**Variables:**
- `PRF` — frecuencia de repetición de pulsos (Hz)
- `B_D` — ancho de banda Doppler de la señal de acimut (Hz)
- `v` — velocidad de la plataforma (m/s)
- `d_a` — longitud de antena en acimut (m)

---

### Diferencia de rango interferométrica

$$\Delta r \simeq \frac{B_\perp}{r_0 \sin(\theta_i)} \cdot \Delta h$$

**Variables:**
- `B_\perp` — línea de base perpendicular a la línea de visión (m)
- `\theta_i` — ángulo de incidencia local
- `\Delta h` — diferencia de altura del terreno (m)

---

### Fase interferométrica

$$\Delta\varphi = m\frac{2\pi}{\lambda}\Delta r$$

**Variables:**
- `m` — factor de modo: 1 para interferómetro single-pass con un TX y dos RX; 2 para repeat-pass
- `\lambda` — longitud de onda (m)

---

## 3. Suposiciones del Modelo

1. Plataforma en movimiento rectilíneo uniforme a velocidad $v$ y altura $H$ constantes (modelo simplificado).
2. Aproximación de Taylor de segundo orden para $r(t)$ (válida cuando $vt \ll r_0$); el término cuadrático domina la historia de fase.
3. Modelo de señal de acimut como chirp lineal FM ("narrow-band approximation"): la frecuencia Doppler varía linealmente con el tiempo lento.
4. Modelo "stop-and-go": la plataforma se supone estacionaria durante la transmisión y recepción de cada pulso; válido cuando $v \ll c_0$.
5. La resolución en rango depende únicamente del ancho de banda $B_r$ del chirp transmitido, independientemente del rango.
6. La resolución en acimut de un SAR enfocado es $d_a/2$, independiente del rango (ventaja fundamental sobre SLAR).
7. Existe una limitación fundamental en sistemas SAR de canal único: mejorar la resolución en acimut (mayor PRF) reduce el ancho de faja en rango (swath width), y viceversa.
8. La reflectividad de la escena $\sigma_0$ es compleja (amplitud y fase); la fase contiene información de la topografía (InSAR) y los desplazamientos (DInSAR).
9. Para polarimetría monoestática, la matriz de dispersión $[S]$ es simétrica ($S_{HV} = S_{VH}$), reduciendo los parámetros independientes a cinco.

---

## 4. Notas Adicionales

- El paper es un tutorial exhaustivo publicado como artículo de revista de IEEE GRSS. Cubre los principios básicos de SAR y también técnicas avanzadas: polarimetría (descomposición de Freeman-Durden, entropía H/A/alpha), interferometría (InSAR, DInSAR, PSInSAR), interferometría SAR polarimétrica (Pol-InSAR), tomografía SAR, y tecnologías futuras (beamforming digital, MIMO, configuraciones bi/multi-estáticas).
- La resolución SAR fundamental se expresa en el dominio del número de onda como: área de celda de resolución $\Delta A_{SAR} = \frac{\lambda_c}{2(\vartheta_2 - \vartheta_1)} \cdot \frac{c}{2B}$, donde $\vartheta_2 - \vartheta_1$ es el ángulo de apertura y $B$ es el ancho de banda.
- En el contexto de **tomografía SAR**, el paper menciona que múltiples pasadas con diferentes líneas de base permiten reconstruir el perfil vertical de retrodispersión (perfil de altura), superando la ambigüedad de layover en zonas urbanas y vegetadas. La tomografía holográfica permite obtener una vista de 360°.
- Sistemas de referencia SAR mencionados: ERS-1/2, JERS-1, Radarsat-1, SRTM, ENVISAT/ASAR, ALOS/PalSAR, TerraSAR-X/TanDEM-X, Radarsat-2, COSMO-SkyMed, Sentinel-1.
- Bandas de frecuencia comúnmente usadas en SAR: P (0.5–0.25 GHz), L (2–1 GHz), S (3.75–2 GHz), C (7.5–3.75 GHz), X (12–7.5 GHz), Ku (17.6–12 GHz), Ka (40–25 GHz).
