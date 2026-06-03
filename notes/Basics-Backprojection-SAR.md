# Basics of Backprojection Algorithm for Processing Synthetic Aperture Radar Images

**Autores:** Armin W. Doerry (Sandia National Laboratories), Edward E. Bishop, John A. Miller (General Atomics Aeronautical Systems, Inc.)  
**Año:** 2016 (Febrero)  
**Fuente/Publicación:** Sandia Report SAND2016-1682, Unlimited Release, Sandia National Laboratories, Albuquerque, NM. ResearchGate: https://www.researchgate.net/publication/295905425

---

## 1. Geometría del Sistema

El sistema es un SAR **monoestático**, cubriendo tanto modos **spotlight** como **VideoSAR** (Circular SAR / Motion Imagery SAR, MISAR).

- **Sistema de coordenadas:** Marco $(x, y, z)$ centrado en el **Scene Reference Point (SRP)**, también llamado Motion Compensation Point (MCP). El plano $xy$ es localmente horizontal respecto al SRP.
- **Vector de posición del radar:** $\mathbf{r}_{c,n}$ = vector del SRP al radar en el pulso $n$.
- **Vector de posición del objetivo:** $\mathbf{s}$ = vector del SRP al dispersor objetivo.
- **Coordenadas del radar:** expresadas mediante ángulo de apertura $\alpha_n$ (ángulo azimutal en la apertura) y ángulo de grazing $\psi_n$ (ángulo de elevación respecto a la horizontal):

$$x_c = |\mathbf{r}_{c,n}|\cos\psi_n \sin\alpha_n, \quad y_c = -|\mathbf{r}_{c,n}|\cos\psi_n \cos\alpha_n, \quad z_c = |\mathbf{r}_{c,n}|\sin\psi_n$$

**Modos:**
- **Spotlight mode:** Cada pulso se asocia a un único número de índice conocido antes de la colección; permite máxima compensación de movimiento via modulación de parámetros del waveform.
- **VideoSAR / Circular SAR (MISAR):** Un pulso puede usarse en múltiples imágenes SAR, con diferentes índices según el subconjunto de pulsos seleccionado para cada parche de imagen.

---

## 2. Ecuaciones de Resolución SAR

### Señal transmitida (TX) — pulso LFM chirp

$$x_T(t, n) = A_T \,\mathrm{rect}\!\left(\frac{t - t_n}{T_{TX}}\right)\cos\!\left(\Phi_T(t - t_n, n)\right)$$

**Variables:**
- `n` — índice de muestra azimutal, $-N/2 \leq n < N/2$
- `N` — número total de muestras azimutales
- `t_n` — tiempo de referencia para el $n$-ésimo pulso TX
- `A_T` — amplitud arbitraria del pulso TX
- `T_{TX}` — duración del pulso TX (s)
- `\Phi_T(t,n)` — función de fase del pulso TX

---

### Función de fase LFM chirp

$$\Phi_T(t, n) = \phi_{T,n} + \omega_{T,n}\,t + \frac{\gamma_{T,n}}{2}t^2$$

**Variables:**
- `\phi_{T,n}` — fase de referencia del $n$-ésimo pulso
- `\omega_{T,n}` — frecuencia angular central de referencia del $n$-ésimo pulso (rad/s)
- `\gamma_{T,n}` — tasa de chirp de referencia del $n$-ésimo pulso (rad/s²)

---

### Señal de eco recibida de un dispersor puntual

$$x_R(t, n) = A_R \,\mathrm{rect}\!\left(\frac{t - t_n - t_{s,n}}{T_{TX}}\right)\cos\!\left(\Phi_T(t - t_n - t_{s,n}, n)\right)$$

**Variables:**
- `A_R` — amplitud de la señal de eco recibida
- `t_{s,n}` — tiempo de retraso del eco para el dispersor puntual en el pulso $n$

---

### Señal de video (Phase History data) — modelo baseband demodulado

$$x_V(i, n) = A_R \exp\!\left\{j\!\left[\!-(\omega_{T,n} + \gamma_{T,n}T_s i)\,\kappa_n\,\tau_{sc,n} + \frac{\gamma_{T,n}\kappa_n}{2}\left(\tau_{sc,n}\right)^2\right]\right\}$$

donde el **Phase History (PH) data model** refinado es:

$$\Phi_V(i, n) = \left(-\frac{2}{c}(\omega_0 + \gamma_0 T_s i)\,\kappa_n\,s_{r,n} + \frac{2\gamma_0\kappa_n}{c^2}(s_{r,n})^2\right)$$

**Variables:**
- `i` — índice de muestra en rango (fast-time), $-I/2 \leq i < I/2$
- `T_s` — espaciado de muestras en fast-time (s)
- `\tau_{sc,n} = t_{s,n} - t_{c,n}` — retraso relativo entre dispersor y SRP para el pulso $n$
- `s_{r,n} = |\mathbf{r}_{s,n}| - |\mathbf{r}_{c,n}|` — desplazamiento de rango slant del objetivo respecto al SRP
- `\kappa_n` — factor de compensación de proyección de número de onda ("wavenumber projection compensation factor")
- `\omega_0` — frecuencia angular nominal de referencia de la apertura (rad/s)
- `\gamma_0` — tasa de chirp nominal de referencia de la apertura (rad/s²)
- `c` — velocidad de la luz (m/s)
- El segundo término es el **Residual Video Phase Error (RVPE)**, artefacto del procesado stretch

---

### Modelo PH completo para un dispersor en s

$$x_V(i, n) = A_R \exp j\!\left(\!\!-\frac{2}{c}(\omega_0 + \gamma_0 T_s i)\,\kappa_n\!\left(|\mathbf{r}_{c,n} - \mathbf{s}| - |\mathbf{r}_{c,n}|\right) + \frac{2\gamma_0\kappa_n}{c^2}\!\left(|\mathbf{r}_{c,n} - \mathbf{s}| - |\mathbf{r}_{c,n}|\right)^2\right)$$

---

### Retraso de eco (echo delay)

$$\tau_{sc,n} = (t_{s,n} - t_{c,n}) = \frac{2}{c}\,s_{r,n}$$

Expandiendo la geometría:

$$\tau_{sc,n} = \frac{2}{c}\left(\sqrt{(x_c - s_x)^2 + (y_c - s_y)^2 + (z_c - s_z)^2} - \sqrt{x_c^2 + y_c^2 + z_c^2}\right)$$

---

### Factor de modulación de parámetros del waveform

$$\kappa_n = \frac{\cos\psi_0}{\cos\psi_n \cos\alpha_n}$$

**Variables:**
- `\psi_0` — ángulo de grazing nominal (en el centro de la apertura)
- `\psi_n` — ángulo de grazing del $n$-ésimo pulso
- `\alpha_n` — ángulo de apertura del $n$-ésimo pulso respecto al centro de apertura

---

### Filtro adaptado (Matched Filter) — salida de imagen

$$y(\hat{\mathbf{s}}) = \sum_i \sum_n \left\{x_V(i, n)\, h(i, n, \hat{\mathbf{s}})^*\right\}$$

donde la función filtro es:

$$h(i, n, \hat{\mathbf{s}}) = e^{j\!\left(\!-\frac{2}{c}(\omega_0+\gamma_0 T_s i)\kappa_n\!\left(|\mathbf{r}_{c,n}-\hat{\mathbf{s}}|-|\mathbf{r}_{c,n}|\right) + \frac{2\gamma_0\kappa_n}{c^2}\!\left(|\mathbf{r}_{c,n}-\hat{\mathbf{s}}|-|\mathbf{r}_{c,n}|\right)^2\right)}$$

**Variables:**
- `\hat{\mathbf{s}}` — posición de prueba (test location / píxel de imagen)
- `*` — denota complejo conjugado
- La doble suma sobre $i$ (rango) y $n$ (acimut) implementa la correlación 2D del filtro adaptado

---

## 3. Suposiciones del Modelo

1. **Modelo "stop & go":** El tiempo de retraso del eco $t_{s,n}$ se asume constante durante la recepción del pulso (la plataforma no se mueve durante el pulso). Las correcciones para el efecto Doppler intra-pulso son posibles pero generalmente menores.
2. **Procesado stretch (de-chirp):** Se asume que la señal LO (Local Oscillator) es también un chirp con la misma fase de referencia, frecuencia y tasa de chirp que el pulso TX; el tiempo de referencia $t_m = t_{c,n}$ se ajusta al retraso de eco al SRP.
3. Todas las muestras se toman en el intervalo definido por los retornos superpuestos de todos los objetivos de interés; no hay muestras "vacías".
4. El factor $\kappa_n$ puede tomar cualquier valor (incluyendo $\kappa_n = 1$ para colección sin modulación de parámetros); el algoritmo FBP lo accomoda en todos los casos.
5. La imagen SAR se forma realizando la doble suma del filtro adaptado para cada posición de prueba $\hat{\mathbf{s}}$ en una rejilla 2D; esto es computacionalmente intensivo si se implementa directamente.
6. El procesado FBP permite reducir la carga computacional mediante descomposición en etapas y representación en coordenadas locales polares de las imágenes intermedias.
7. La corrección RVPE (Residual Video Phase Error) es necesaria en el preprocesado antes de la formación de imagen para compensar el segundo término de fase cuadrática en $s_{r,n}$.

---

## 4. Notas Adicionales

- El reporte es un **manual práctico de implementación** del algoritmo Filtered Backprojection (FBP) para formación de imágenes SAR en tiempo real. Cubre preprocesado de Phase History data (rampa de filtro, compresión en rango, compensación RVPE, corrección por pérdida de rango), formación de imagen en la rejilla de muestreo, y post-procesado (correcciones de amplitud, correcciones de espectro de imagen, autofocus).
- El algoritmo FBP tiene sus raíces en **tomografía computarizada (CT)**, siendo matemáticamente equivalente al filtro adaptado para imágenes SAR.
- El paper fundamental que estableció BP como algoritmo viable de formación de imágenes SAR fue Munson et al. (Paper 5 de esta colección), y un buen tutorial de implementación en MATLAB es Gorham y Moore.
- Para **VideoSAR / MISAR**, el FBP es especialmente ventajoso porque no requiere conocer el índice del pulso antes de la colección, a diferencia de los algoritmos range-Doppler tradicionales.
- El sistema de coordenadas 3D $(x, y, z)$ usado en este reporte es completamente general: no requiere ninguna aproximación de "apertura plana" ni de campo lejano para la formación de imagen. La aproximación de campo lejano ($s_{r,n}$ aproximado según Eq. 29) se menciona solo para propósitos informativos.
- La **modulación de parámetros del waveform** (ajuste de $\omega_{T,n}$ y $\gamma_{T,n}$ en función de $\kappa_n$) es una conveniencia para algoritmos range-Doppler tradicionales (PFA, OSA); el FBP no la requiere pero puede acomodarla.
