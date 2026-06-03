# Synthetic-Aperture Radar Processing Using Fast Factorized Back-Projection

**Autores:** Lars M. H. Ulander, Hans Hellsten, Gunnar Stenstrom  
**Año:** 2003  
**Fuente/Publicación:** IEEE Transactions on Aerospace and Electronic Systems, Vol. 39, No. 3, July 2003. IEEE Log No. T-AES/39/3/816615.

---

## 1. Geometría del Sistema

El sistema considerado es un SAR **monoestático de apertura lineal**, con la antena moviéndose a lo largo de una trayectoria nominalmente rectilínea. La geometría general permite trayectorias no lineales.

- **Sistema de coordenadas:** cilíndrico $(\rho, \theta, x)$ donde $x$ es el eje a lo largo de la trayectoria (along-track), $\rho$ es el rango slant, y $\theta$ el ángulo de elevación.
- Para la apertura lineal existe simetría cilíndrica alrededor del eje de pista, lo que permite simplificar el modelo de eco.
- La historia de rango para un objeto en $(\rho_0, \theta_0, x_0)$ tiene dependencia hiperbólica en $x$:

$$g(x, R) = \frac{p\!\left(R - \sqrt{(x-x_0)^2 + \rho_0^2}\right)}{\left(\sqrt{(x-x_0)^2 + \rho_0^2}\right)^2}$$

- Se cubre el caso general de apertura no lineal (trayectorias curvas de aviones), motivado por el sistema CARABAS SAR (VHF/UHF ultra-wideband de la agencia sueca FOI), que requiere apertura de haz ancha.

**Dominio de números de onda (k-space):**

La función de transferencia del sistema SAR, antes del procesado, forma un sector de anillo (annular sector) en el dominio bidimensional de números de onda $(k_R, k_x)$:
- Anillo definido por frecuencias mínima y máxima: $2\omega_0/c$ y $2\omega_1/c$.
- Sector angular definido por el intervalo de ángulo Doppler $\vartheta_2 - \vartheta_1$.

La resolución SAR fundamental es:

$$\Delta A_{SAR} = \frac{\lambda_c}{2(\vartheta_2 - \vartheta_1)} \cdot \frac{c}{2B}$$

---

## 2. Ecuaciones de Resolución SAR

### Resolución SAR fundamental (dominio de número de onda)

$$\Delta A_{SAR} = \frac{\lambda_c}{2(\vartheta_2 - \vartheta_1)} \cdot \frac{c}{2B}$$

**Variables:**
- `\lambda_c` — longitud de onda del centro de frecuencia: $\lambda_c = c/f_c$, donde $f_c = (f_{max} + f_{min})/2$
- `\vartheta_2 - \vartheta_1` — ángulo de apertura total (intervalo de ángulo Doppler, en radianes)
- `B` — ancho de banda del sistema: $B = f_{max} - f_{min}$ (Hz)
- `c` — velocidad de la luz (m/s)

---

### Modelo de eco radar para objeto puntual (caso monoestático)

$$g(\vec{r}, R) = \frac{p(R - |\vec{r} - \vec{r}_o|)}{|\vec{r} - \vec{r}_o|^2}$$

**Variables:**
- `p(R)` — pulso radar band-limited tras compresión en rango
- `R` — rango (tiempo de retraso escalado por $c/2$)
- `\vec{r}` — posición de la antena
- `\vec{r}_o` — posición del objeto puntual
- La dependencia $|\vec{r} - \vec{r}_o|^{-2}$ refleja la atenuación de espacio libre (one-way $R^{-2}$ para el campo, two-way $R^{-4}$ para la potencia; aquí se trabaja con amplitud del campo)

---

### Integral de backprojection directa

$$f(\vec{r}_o) = \int p\!\left(R - |\vec{r} - \vec{r}_o|\right) d\vec{r}$$

**Variables:**
- `f(\vec{r}_o)` — imagen SAR reconstruida (reflectividad estimada en posición $\vec{r}_o$)
- La integración se realiza sobre todas las posiciones $\vec{r}$ de la apertura (a lo largo de la trayectoria de vuelo)
- Para cada punto imagen, se suma la contribución del eco evaluado en el retardo correspondiente a la distancia $|\vec{r} - \vec{r}_o|$

---

### Número de operaciones: backprojection directa vs. FBP

**Backprojection directa:**
$$N_{ops}^{direct} \propto N^3$$

**Fast Factorized Back-Projection (FBP) multi-nivel:**
$$N_{ops}^{FBP} \propto N^2 \log_2(N)$$

**Variables:**
- `N` — número de píxeles en una dimensión (imagen de $N \times N$ píxeles con $N$ posiciones de apertura)
- La FBP logra una reducción de 2–3 órdenes de magnitud en tiempo de procesado para valores típicos de $N$

---

### Modelo de eco en coordenadas cilíndricas (apertura lineal)

$$g(x, R) = \frac{p\!\left(R - \sqrt{(x-x_0)^2 + \rho_0^2}\right)}{\left(\sqrt{(x-x_0)^2 + \rho_0^2}\right)^2}$$

**Variables:**
- `x` — posición a lo largo de la pista (along-track coordinate)
- `x_0` — posición en acimut del objeto
- `\rho_0` — rango slant mínimo del objeto (distancia en el plano perpendicular a la pista)
- `R` — variable de rango (tiempo de retraso × c/2)

---

## 3. Suposiciones del Modelo

1. Radar monoestático: antena transmisora y receptora colocadas en el mismo punto; las antenas TX y RX son colocadas ("collocated").
2. Medio de propagación lineal, isotrópico, homogéneo y no dispersivo (velocidad de onda $c$ constante en el espacio).
3. La escena es una colección de objetos de dispersión simple ("Born approximation"): no hay dispersión múltiple, la superposición es válida.
4. Aproximación "start-stop": la plataforma se supone estacionaria durante la transmisión y recepción de cada pulso; válida cuando $v \ll c$.
5. Las respuestas de antena y objeto (en función de frecuencia, polarización y dirección) se asumen iguales a 1 para simplicidad; en la práctica se evalúan mediante métodos de fase estacionaria.
6. Para la apertura lineal existe simetría cilíndrica: el ángulo de elevación $\theta_o$ del objeto no está determinado (ambigüedad izquierda-derecha), lo que es mitigado por la directividad de la antena en elevación.
7. El algoritmo FBP se basa en la factorización recursiva de la integral de backprojection, representando las imágenes intermedias en coordenadas polares locales para reducir la carga computacional.
8. La factorización es exacta cuando los subaperturas son suficientemente pequeñas; la acumulación de errores de fase se controla manteniendo el error por debajo de $\pi/4$ por etapa.

---

## 4. Notas Adicionales

- El paper introduce el algoritmo **Fast Factorized Back-Projection (FBP)** como alternativa computacionalmente eficiente al backprojection directo, logrando complejidad $O(N^2 \log N)$ frente a $O(N^3)$.
- El algoritmo FBP es una generalización unificadora que incluye como casos especiales el backprojection directo (single-stage), los algoritmos two-stage y multiple-stage, y los algoritmos de fast transform para apertura lineal.
- La representación de imágenes intermedias en **coordenadas polares locales** es el elemento clave que permite la factorización eficiente.
- El algoritmo es validado con datos reales del sistema **CARABAS SAR** (VHF ultra-wideband, Swedish Defence Research Agency FOI), demostrando reducción de 2–3 órdenes de magnitud en tiempo de procesado.
- A diferencia de los algoritmos de dominio frecuencial (Polar Format, omega-k), FBP maneja naturalmente geometrías de apertura no lineal y no requiere la aproximación de campo lejano para la corrección de migración de celda de rango.
- La conexión con la tomografía computarizada (CT): la backprojection en SAR es matemáticamente equivalente a la inversión de Fourier-Hankel para SAR de apertura lineal, análogo a la retro-proyección en CT médica.
- El sistema de transferencia SAR en el dominio de número de onda ocupa una región anular sectorial, cuyos límites radiales están dados por las frecuencias mínima y máxima del chirp y cuyos límites angulares por el ángulo de apertura.
