# Single and Multipolarimetric P-Band SAR Tomography of Subsurface Ice Structure

**Autores:** Francesco Banda, Jørgen Dall (Member, IEEE), Stefano Tebaldini  
**Año:** 2015 (aceptado noviembre 2015)  
**Fuente/Publicación:** IEEE Transactions on Geoscience and Remote Sensing, DOI: 10.1109/TGRS.2015.2506399

---

## 1. Geometría del Sistema

El sistema opera en el marco de la **campaña ESA IceSAR 2012** sobre el glaciar Russell en el suroeste de Groenlandia (K-transecto, ~67°N):

- **Plataforma:** avión airborne POLARIS (radar P-band, 435 MHz) del DTU (Technical University of Denmark).
- **Modo:** SAR tomográfico (TomoSAR) — múltiples pasadas sobre la misma escena a lo largo de una dirección transversal, formando una apertura sintética 2D en el plano elevación–azimut.
- **Geometría:** sistema cartesiano $(x, y, z)$ donde $x$ = azimut, $y$ = ground range, $z$ = elevación. La elevación se expresa respecto al elipsoide WGS84.
- **Ángulos de incidencia:** 25°–45°.
- **Altitud nominal del sensor:** 4 km sobre la superficie del hielo.
- **Sitios investigados:**
  - **SHR** (~50°W, ~700 m altitud): zona de ablación, 8 pasadas NENW + 8 pasadas SWSE (mayo y junio 2012), polarización HH.
  - **S10** (~47°W, ~1850 m altitud): zona de percolación/acumulación, 4 pasadas, polarización completa HH/HV/VV/VH.

**Parámetros del sistema (SHR y S10):**

| Parámetro | SHR | S10 |
|---|---|---|
| Altitud sensor | 4 km | 4 km |
| Frecuencia portadora | 435 MHz | 435 MHz |
| Longitud de onda | 68 cm | 68 cm |
| Ancho de banda | 85 MHz | 85 MHz |
| Frecuencia muestreo | 125 MHz | 125 MHz |
| Res. slant range | 1.8 m | 1.8 m |
| Res. azimut | 5 m | 5 m |
| N° pasadas | 8 | 4 |
| Apertura baseline total | 20 m | 30 m |
| Baseline nominal | 0 m | 10 m |
| Polarización | HH | HH/HV/VV/VH |

---

## 2. Ecuaciones de Resolución SAR

### 2.1 Resolución en Slant Range

$$\delta_r = \frac{c}{2W}$$

**Variables:**
- `c` — velocidad de la luz
- `W` — ancho de banda del pulso (85 MHz) → $\delta_r \approx 1.8$ m

### 2.2 Resolución Tomográfica Vertical (Fourier Beamforming)

La resolución vertical en TomoSAR (apertura en elevación) es:

$$\delta_z = \frac{\lambda R \sin\theta}{2 B_{\perp,\mathrm{max}}}$$

**Variables:**
- `\lambda` — longitud de onda (68 cm)
- `R` — slant range
- `\theta` — ángulo de incidencia
- `B_{\perp,\mathrm{max}}` — apertura máxima de la baseline perpendicular

> La resolución vertical varía con $R$ y $\theta$, siendo coarsa en la mayor parte de la escena para esta campaña (apertura de vuelo ajustada para penetración de cientos de metros, no resolución).

### 2.3 Profundidad Aparente (sin corrección de refracción)

$$z_{\mathrm{air}} = -\frac{c}{2}(t - t_0)\cos\theta$$

**Variables:**
- `t` — tiempo de llegada del eco
- `t_0` — tiempo de llegada del eco de la interfaz hielo-aire
- `\theta` — ángulo de incidencia

### 2.4 Profundidad Real en el Hielo (Interfaz Plana)

$$z = -\frac{v}{2}(t - t_0)\cos\theta_s$$

$$z = z_{\mathrm{air}} \cdot \frac{1}{\sqrt{\varepsilon_r}} \cdot \frac{\cos\theta_s}{\cos\theta}$$

con $\varepsilon_r = 3$ para hielo → $v = c/\sqrt{3} \approx 1.73 \times 10^8$ m/s.

### 2.5 Ley de Snell en Interfaz Aire–Hielo

$$\sin\theta_s = \frac{\sin\theta}{\sqrt{\varepsilon_r}} = \frac{\sin\theta}{\sqrt{3}}$$

### 2.6 Corrección para Interfaz Inclinada (pendiente $\alpha$)

$$z = z_{\mathrm{air}} \cdot \frac{1}{\sqrt{\varepsilon_r}} \cdot \frac{\cos(\theta_s + \alpha)}{\cos\theta}$$

$$\theta_s = \arcsin\!\left(\frac{\sin(\theta - \alpha)}{\sqrt{\varepsilon_r}}\right)$$

**Variables:**
- `\alpha` — pendiente local de la interfaz hielo-aire en el plano $y$–$z$

### 2.7 Número de Onda Vertical (Height-to-Phase)

$$k_z = \frac{4\pi B_\perp}{\lambda R \sin\theta}$$

**Variables:**
- `B_\perp` — baseline perpendicular entre dos pasadas
- `k_z` — factor de conversión altura-fase (rad/m)

### 2.8 Modelo de Coherencia — Inversión de Profundidad de Dispersión

La matriz de coherencia del stack multibaseline se modela como suma de contribuciones:

$$\mathbf{R} = \sigma_s \mathbf{R}_s + \sigma_v \mathbf{R}_v + \sigma_n \mathbf{R}_n$$

**Variables:**
- `\mathbf{R}_s = \mathbf{1}_N` — matriz $N \times N$ de unos (superficie coherente)
- `\mathbf{R}_v^{mn} = 1/(1 + jk_z^{mn} d_p)` — matriz de volumen con perfil de extinción exponencial
- `\mathbf{R}_n = \mathbf{I}_N` — matriz identidad (ruido no correlado)
- `d_p` — profundidad de dispersión de dos vías (two-way scattering depth)
- `\sigma_s, \sigma_v, \sigma_n` — potencias de superficie, volumen y ruido

### 2.9 Descomposición Polarimétrica (Algebraic Synthesis — S10)

La matriz de covarianza polarimétrica multibaseline se modela como:

$$\mathbf{W} = \mathbb{E}[\mathbf{d}\mathbf{d}^H] \cong \mathbf{W}_2 = \sum_{k=1}^{2} \mathbf{C}_k(a,b) \otimes \mathbf{R}_k(a,b)$$

**Variables:**
- `\mathbf{d}` — vector $3N \times 1$ de píxeles complejos para $N$ adquisiciones y 3 canales de polarización
- `\mathbf{C}_k` — matrices $3 \times 3$ de firma polarimétrica de los mecanismos de dispersión (SMs)
- `\mathbf{R}_k` — matrices $N \times N$ de coherencias interferométricas del SM $k$-ésimo
- `\otimes` — producto de Kronecker

---

## 3. Suposiciones del Modelo

1. La interfaz hielo-aire es **localmente plana** (primera aproximación; se extiende a superficies inclinadas).
2. Permitividad relativa del hielo: $\varepsilon_r = 3$ (valor promedio estándar para hielo seco de glaciar).
3. El perfil de extinción del hielo sigue un modelo **exponencial** (análogo al usado en forests).
4. **Reciprocidad:** HV = VH (válido para medios naturales).
5. **Invarianza estructural respecto a polarización** para cada mecanismo de dispersión (SM).
6. **Estacionariedad** de los datos entre pasadas (cada SM puede describirse independientemente por canal polarimétrico).
7. La resolución tomográfica vertical limitada por la apertura de la baseline es la restricción principal para profundidades de penetración < resolución.

---

## 4. Notas Adicionales

- El paper demuestra la **viabilidad del TomoSAR P-band para hielo subsuperficial** sobre dos sitios con condiciones físicas muy diferentes:
  - **SHR (ablación):** dispersión principalmente superficial; penetración limitada (~20–40 m verdaderos, ~0–40 m aparentes); morfología subsuperficial clara (estructuras tipo wishbone y flower).
  - **S10 (acumulación/percolación):** decorrelación volumétrica hasta ~100 m; separación volumen/superficie con Algebraic Synthesis exitosa.
- El **estimador de Capon** supera el beamforming de Fourier en resolución y supresión de lóbulos laterales.
- El estimador **MUSIC** se usa para cuantificar la profundidad de dispersión (boundaries del tomograma subsuperficial).
- Los efectos de refracción deben corregirse para obtener profundidades reales: blancos a ~23 m aparente equivalen a ~40 m real en hielo.
- Sin corrección de refracción y con velocidad $c$ constante, el impacto en la **resolución vertical 3D es mínimo** (confirmado numéricamente), pero el **mislocation es significativo**.
- Futuras misiones BIOMASS (ESA, P-band) y SAOCOM-CS (L-band) se beneficiarán de estas técnicas para criósfera y biomasa forestal.
- Los datos IceSAR 2012 se adquirieron en mayo y junio 2012; la diferencia entre meses refleja cambios en la tasa de fusión superficial y su efecto en la penetración radar.
