# Refraction Effect Correction in Back Projection Algorithm for Subsurface SAR Imaging

**Autores:** (Conferencia IMOC 2023 — autores no visibles en el texto extraído; paper relacionado con el grupo de Goes/UNICAMP)  
**Año:** 2023  
**Fuente/Publicación:** MOMAG/IMOC 2023 — Conferencia Internacional de Microondas y Optoelectrónica (Brasil)

---

## 1. Geometría del Sistema

El sistema modela la propagación de ondas radar a través de dos medios con diferentes velocidades de propagación:

- **Medio 1 (aire):** el radar se desplaza a altura $H$ sobre la superficie. La onda se propaga a velocidad $c$ (velocidad de la luz).
- **Medio 2 (subsuelo/hielo/suelo):** la onda penetra en el subsuelo con velocidad de propagación $v = c/\sqrt{\varepsilon_r}$, donde $\varepsilon_r$ es la permitividad relativa del medio.

La geometría considera una interfaz plana entre el aire y el subsuelo. Al penetrar la interfaz, la onda sufre **refracción** según la Ley de Snell, cambiando su dirección de propagación y su velocidad. El algoritmo Back Projection estándar, diseñado para propagación en aire (velocidad $c$), introduce errores de localización (mislocation) y desenfoque (defocusing) en targets subsuperficiales si no se corrige la refracción.

El escenario de referencia es tomografía SAR de hielo glacial (similar al de Banda et al. 2015), con:
- Plataforma aérea (avión o dron) volando sobre superficie de hielo
- Permitividad relativa del hielo: $\varepsilon_r \approx 3$
- Blancos subsuperficiales a profundidades de decenas de metros

---

## 2. Ecuaciones de Resolución y Modelo de Refracción

### 2.1 Velocidad de Propagación en el Subsuelo

$$v = \frac{c}{\sqrt{\varepsilon_r}}$$

**Variables:**
- `c` — velocidad de la luz en el vacío ($3 \times 10^8$ m/s)
- `\varepsilon_r` — permitividad relativa del subsuelo (hielo: $\varepsilon_r \approx 3$)
- `v` — velocidad de fase de la onda en el subsuelo

### 2.2 Ley de Snell en la Interfaz Aire–Subsuelo

$$c \sin\theta_s = v \sin\theta$$

equivalentemente:

$$\sin\theta_s = \frac{1}{\sqrt{\varepsilon_r}}\sin\theta$$

**Variables:**
- `\theta` — ángulo de incidencia en el aire (respecto a la normal)
- `\theta_s` — ángulo de refracción en el subsuelo (respecto a la normal)

### 2.3 Profundidad Aparente (Supuesto $c$ constante, interfaz plana)

Sin corrección de refracción, el BP ubica el blanco a una profundidad aparente:

$$z_{\mathrm{air}} = -\frac{c}{2}(t - t_0)\cos\theta$$

**Variables:**
- `t` — tiempo de llegada del eco
- `t_0` — tiempo de llegada del eco de la interfaz aire-subsuelo ($z = 0$)
- `\theta` — ángulo de incidencia en el radar

### 2.4 Profundidad Real en el Subsuelo (Interfaz Plana)

$$z = -\frac{v}{2}(t - t_0)\cos\theta_s$$

### 2.5 Corrección de Profundidad — Interfaz Plana

Relacionando $z$ con $z_{\mathrm{air}}$ mediante las ecuaciones anteriores y la Ley de Snell:

$$z = z_{\mathrm{air}} \cdot \frac{1}{\sqrt{\varepsilon_r}} \cdot \frac{\cos\theta_s}{\cos\theta}$$

### 2.6 Corrección de Profundidad — Interfaz Inclinada (pendiente $\alpha$)

Para una interfaz con pendiente local $\alpha$ en el plano $y$–$z$:

$$z = z_{\mathrm{air}} \cdot \frac{1}{\sqrt{\varepsilon_r}} \cdot \frac{\cos(\theta_s + \alpha)}{\cos\theta}$$

donde $\theta_s = \arcsin\!\left(\sin(\theta - \alpha)/\sqrt{\varepsilon_r}\right)$ (ley de Snell modificada para superficie inclinada).

### 2.7 Distancia de Propagación Corregida en BP Refractive

El algoritmo BP refractive calcula la distancia total recorrida por la onda como suma de dos tramos:

$$R_{\mathrm{total}} = R_{\mathrm{aire}} + \frac{R_{\mathrm{sub}}}{1} \cdot \sqrt{\varepsilon_r}$$

equivalentemente, el retardo de fase total para la compensación en el BP es:

$$\varphi = \frac{4\pi}{\lambda}\left(R_{\mathrm{aire}} + \sqrt{\varepsilon_r}\, R_{\mathrm{sub}}\right)$$

**Variables:**
- `R_{\mathrm{aire}}` — distancia recorrida en aire (antena a punto de entrada en la interfaz)
- `R_{\mathrm{sub}}` — distancia recorrida en el subsuelo (interfaz a target)
- `\lambda` — longitud de onda en el aire

### 2.8 Punto de Entrada en la Interfaz (Trazado de Rayos)

El punto de cruce de la onda con la interfaz se determina minimizando el tiempo de tránsito total (principio de Fermat):

$$\mathbf{p}^* = \arg\min_{\mathbf{p} \in \text{interfaz}} \left[ \frac{\|\mathbf{r}_k - \mathbf{p}\|}{c} + \frac{\|\mathbf{p} - \mathbf{h}_m\|}{v} \right]$$

**Variables:**
- `\mathbf{r}_k` — posición del radar en el pulso $k$
- `\mathbf{h}_m` — posición del píxel subsuperficial $m$
- `\mathbf{p}` — punto candidato sobre la interfaz

---

## 3. Suposiciones del Modelo

1. La interfaz aire-subsuelo es **localmente plana** (o con pendiente constante $\alpha$ localmente).
2. La permitividad relativa $\varepsilon_r$ es **homogénea** en el subsuelo (sin gradientes internos).
3. Se desprecia la atenuación del medio (solo se modela el cambio de velocidad de fase, no la absorción).
4. El modelo de dos capas (aire + subsuelo) es suficiente para describir la geometría de propagación.
5. El DEM de la superficie es conocido (o estimado) para determinar el punto de entrada de los rayos en la interfaz.

---

## 4. Notas Adicionales

- El algoritmo BP refractive es una extensión directa del BP estándar: se reemplaza la distancia euclidiana $R_{m,k} = \|\mathbf{h}_m - \mathbf{r}_k\|$ por la distancia de camino óptico equivalente a través de los dos medios.
- Sin corrección, los blancos subsuperficiales aparecen desplazados hacia arriba (más cerca de la superficie de lo real) y con menor nitidez (defocusing) a medida que aumenta la profundidad.
- Para hielo ($\varepsilon_r \approx 3$): la velocidad en el subsuelo es $v \approx 1.73 \times 10^8$ m/s, y los blancos a 40 m de profundidad real aparecen en $z_{\mathrm{air}} \approx 23$ m si no se corrige.
- El paper demuestra la corrección con datos simulados y/o experimentales, mostrando la mejora en localización y nitidez de los blancos tras aplicar el BP refractive.
- Este trabajo es directamente relevante para SAR-GPR embarcado en drones operando sobre suelo, hielo glacial o capas de nieve.
- El trazado de rayos para encontrar $\mathbf{p}^*$ puede resolverse analíticamente en 2D para superficie plana, o numéricamente para geometrías 3D con superficie irregular.
