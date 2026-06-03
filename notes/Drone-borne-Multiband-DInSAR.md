# A Drone-Borne Multiband DInSAR: Results and Applications

**Autores:** (Paper de conferencia — grupo de investigación de SAR con UAV, probablemente DLR o FGAN, Alemania; autores exactos no visibles en el extracto)  
**Año:** ~2020–2022 (paper de conferencia reciente)  
**Fuente/Publicación:** Conferencia IEEE IGARSS o EUSAR — Multiband Drone SAR Interferometry

---

## 1. Geometría del Sistema

El sistema consiste en un **radar SAR multibanda embarcado en un dron (UAV)** que opera en modo **DInSAR (Differential Interferometric SAR)**:

- La plataforma UAV realiza múltiples pasadas sobre la misma escena en diferentes instantes (pases repetidos), generando pares interferométricos.
- El sistema opera en múltiples bandas de frecuencia simultáneamente o en pasadas separadas, permitiendo análisis multiescala de deformación y penetración.
- La geometría de adquisición es **side-looking** con ángulo de incidencia variable entre near-range y far-range.
- La **línea de base interferométrica** (baseline) entre dos pasadas puede ser espacial (diferencia de posición transversal) o temporal (diferencia en tiempo entre adquisiciones).
- En DInSAR diferencial, se elimina la contribución topográfica (usando un DEM externo) para aislar la deformación del terreno entre pasadas.

Parámetros típicos del sistema dron SAR multibanda:
- Altitud de vuelo: decenas a pocos cientos de metros
- Bandas: múltiples (p.ej. X, C, L o combinaciones)
- Resolución: centimétrica a métrica según banda y apertura
- Línea de base espacial: controlada por precisión de vuelo repetido (RTK/PPK GPS)

---

## 2. Ecuaciones de Resolución SAR

### 2.1 Resolución en Rango

$$\delta_r = \frac{c}{2W}$$

**Variables:**
- `c` — velocidad de la luz
- `W` — ancho de banda del chirp transmitido

### 2.2 Resolución en Azimut

$$\delta_a = \frac{d_a}{2}$$

**Variables:**
- `d_a` — dimensión física de la antena en azimut

### 2.3 Fase Interferométrica

La fase medida en un interferograma SAR es:

$$\phi_{\mathrm{int}} = \phi_{\mathrm{topo}} + \phi_{\mathrm{def}} + \phi_{\mathrm{atm}} + \phi_{\mathrm{ruido}}$$

**Variables:**
- `\phi_{\mathrm{topo}}` — contribución topográfica (relacionada con el DEM)
- `\phi_{\mathrm{def}}` — fase de deformación del terreno (señal de interés en DInSAR)
- `\phi_{\mathrm{atm}}` — fase atmosférica
- `\phi_{\mathrm{ruido}}` — ruido de fase

### 2.4 Fase Topográfica

$$\phi_{\mathrm{topo}} = \frac{4\pi}{\lambda} \cdot \frac{B_\perp \cdot h}{R \sin\theta}$$

**Variables:**
- `B_\perp` — línea de base perpendicular al slant range
- `h` — elevación del terreno (del DEM)
- `R` — distancia slant range
- `\theta` — ángulo de incidencia
- `\lambda` — longitud de onda

### 2.5 Fase de Deformación (DInSAR)

Eliminando la componente topográfica con un DEM de referencia:

$$\phi_{\mathrm{def}} = \phi_{\mathrm{int}} - \phi_{\mathrm{topo,DEM}} = \frac{4\pi}{\lambda} \cdot d_{\mathrm{LOS}}$$

**Variables:**
- `d_{\mathrm{LOS}}` — desplazamiento del terreno proyectado en la línea de visión (LOS) del radar
- `\phi_{\mathrm{topo,DEM}}` — fase topográfica simulada con el DEM externo

### 2.6 Altura de Ambigüedad

La altura de ambigüedad define la diferencia de elevación que produce un ciclo completo ($2\pi$) de fase topográfica:

$$h_{\mathrm{amb}} = \frac{\lambda R \sin\theta}{2 B_\perp}$$

**Variables:**
- `h_{\mathrm{amb}}` — altura de ambigüedad (m/ciclo)

### 2.7 Coherencia Interferométrica

$$\gamma = \frac{|\langle s_1 s_2^* \rangle|}{\sqrt{\langle|s_1|^2\rangle\langle|s_2|^2\rangle}}$$

**Variables:**
- `s_1, s_2` — señales SAR complejas de las dos pasadas
- `\langle\cdot\rangle` — promedio espacial (multilook)
- `\gamma \in [0,1]` — coherencia (1 = perfecta, 0 = incoherente)

### 2.8 Ventaja Multiband: Desambiguación de Fase

Usando dos bandas con longitudes de onda $\lambda_1$ y $\lambda_2$, se forma una longitud de onda sintética equivalente:

$$\Lambda = \frac{\lambda_1 \lambda_2}{|\lambda_1 - \lambda_2|}$$

$$\phi_{\Lambda} = \phi_1 - \phi_2 = \frac{4\pi}{\Lambda} d_{\mathrm{LOS}}$$

**Variables:**
- `\Lambda` — longitud de onda sintética (mucho mayor que $\lambda_1$ o $\lambda_2$)
- `\phi_\Lambda` — fase diferencial entre bandas (ambigüedad extendida)

---

## 3. Suposiciones del Modelo

1. La superficie del terreno no cambia significativamente en amplitud entre pasadas (coherencia temporal mantenida).
2. La trayectoria del UAV entre pasadas se controla con precisión sub-centimétrica (RTK/PPK GPS) para garantizar líneas de base pequeñas y repetibilidad.
3. El DEM externo es suficientemente preciso para eliminar la contribución topográfica en DInSAR.
4. Las contribuciones atmosféricas son despreciables o se corrigen con modelos meteorológicos (el efecto es menor en UAV a baja altitud que en satélite).
5. La deformación del terreno es **lenta** respecto al intervalo entre pasadas (compatible con la hipótesis de pequeñas deformaciones en la LOS).
6. La resolución en azimut puede degradarse si la trayectoria del UAV no es perfectamente lineal; se requieren técnicas de compensación de movimiento (MoComp) o autofoco.

---

## 4. Notas Adicionales

- El **DInSAR con dron** ofrece ventajas únicas: resolución espacial centimétrica, flexibilidad de despliegue, baja línea de base espacial controlada, y capacidad de revisita rápida.
- Las aplicaciones incluyen: monitoreo de deformación de estructuras (puentes, presas, edificios), subsidencia de terreno, deslizamientos, vulcanología a pequeña escala, y arqueología (detección de estructuras enterradas por deformación superficial diferencial).
- La operación **multiband** permite combinar resolución (banda X de corta longitud de onda) con penetración (banda L/P de larga longitud de onda), así como desambiguar la fase interferométrica usando longitudes de onda sintéticas $\Lambda$.
- Un desafío fundamental del DInSAR con dron es la **decorrelación temporal**: la vegetación o superficies blandas pueden perder coherencia entre pasadas separadas por horas o días.
- El **procesado BP** es el más adecuado para trayectorias no lineales de UAV; los algoritmos frecuenciales requieren MoComp previo.
- La precisión de posicionamiento del UAV (RTK GPS: ~3 cm) es comparable a la longitud de onda en banda X (~3 cm), lo que hace el control de línea de base crítico para interferometría en X.
