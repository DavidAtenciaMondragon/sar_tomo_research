# Improvement of GPR SAR-based Techniques for Accurate Detection and Imaging of Buried Objects

**Autores:** Marcos González-Díaz, María García-Fernández, Yuri Álvarez-López, Fernando Las Heras  
**Año:** 2019  
**Fuente/Publicación:** IEEE Transactions on Instrumentation and Measurement, vol. 69, no. 6, pp. 3126–3138. DOI: 10.1109/TIM.2019.2930159

---

## 1. Geometría del Sistema

El sistema es un GPR-SAR aerotransportado (airborne, Down-Looking) que utiliza un par de antenas UWB Vivaldi (Tx y Rx) montadas en un UAV, sobrevolando una cuadrícula plana de adquisición en configuración cuasi-monoestática. El dominio de reconstrucción es multilayer: capas horizontales paralelas con diferentes permitividades $\varepsilon_{r,p}$, donde el objeto enterrado se sitúa en la capa Q más profunda.

- **Geometría de referencia (Fig. 1):** el UAV vuela a altura $z = 0$; debajo hay P capas con profundidades $d_1, d_2, \ldots, d_P$ y permitividades $\varepsilon_{r,1}, \varepsilon_{r,2}, \ldots, \varepsilon_{r,Q}$; el objeto enterrado está en la capa Q.
- **Puntos de refracción** $\mathbf{r}_p$ en cada interfaz se calculan aplicando la Ley de Snell de forma iterativa.
- **Antenas:** dos Vivaldi UWB (RFSPACE), operando de 0.1 a 6 GHz; para validación experimental se usa la banda 0.6–6.5 GHz.
- **Adquisición de laboratorio:** rango XYZ de $70 \times 80$ cm con paso $\delta x = \delta y = 2$ cm (1476 puntos), a $z = 120$ cm sobre el suelo. Analizador vectorial de redes (VNA) N5244A PNA-X, calibrado de 0.4 a 7 GHz.
- **Validación realista (playa):** escáner manual a 140 cm de altura, barrido de 120 cm con paso 2 cm.

---

## 2. Ecuaciones de Resolución SAR/GPR

### Algoritmo Delay-And-Sum (DAS) multilayer — reflectividad en capa Q

$$\rho_Q(\mathbf{r'}) = \sum_{m=1}^{M} \sum_{n=1}^{N} E_{scatt}(\mathbf{r}_m, f_n) \prod_{p=1}^{q} e^{+j2\phi_{p,n}}$$

**Variables:**
- `r'` — punto de reconstrucción en la capa Q
- `r_m` — posición del m-ésimo punto de adquisición (antena Tx/Rx en cuasi-monoestático)
- `f_n` — n-ésima frecuencia de trabajo
- `E_scatt(r_m, f_n)` — campo dispersado medido en el m-ésimo punto a la n-ésima frecuencia
- `phi_{p,n}` — desfase debido a la propagación en la p-ésima capa
- `M` — número total de puntos de adquisición
- `N` — número de frecuencias

### Desfase por propagación en la p-ésima capa

$$\phi_{p,n} = k_{p,n} \cdot \|\mathbf{r}_p - \mathbf{r}_{p-1}\|_2$$

**Variables:**
- `k_{p,n}` — número de onda en la p-ésima capa a la n-ésima frecuencia
- `r_p` — punto de refracción en la interfaz de la p-ésima capa (para $p = 1, \ldots, q-1$); $\mathbf{r'} = \mathbf{r}_q$ es el punto de reconstrucción
- `r_{p-1}` — punto de refracción en la capa anterior (o posición de la antena para $p=1$)

### Puntos de refracción en escenario de dos capas (aproximación)

$$\mathbf{r}_1 = \mathbf{r}_{2,1} + \sqrt{\frac{\varepsilon_{r1}}{\varepsilon_{r2}}}(\mathbf{r}_{1,1} - \mathbf{r}_{2,1})$$

**Variables:**
- `r_{2,1}` — proyección del punto de reconstrucción $\mathbf{r}'$ sobre la interfaz (punto de refracción si $\varepsilon_{r2} \gg \varepsilon_{r1}$)
- `r_{1,1}` — proyección de la antena sobre la interfaz (punto de refracción si $\varepsilon_{r1} = \varepsilon_{r2}$)
- `epsilon_{r1}`, `epsilon_{r2}` — permitividades relativas de la capa superior e inferior

### Escenario de tres capas — puntos de refracción

$$\mathbf{r}_1 = \mathbf{r}_{2,1} + \sqrt{\frac{\varepsilon_{r1}}{\varepsilon_{r2}}}(\mathbf{r}_{1,1} - \mathbf{r}_{2,1})$$

$$\begin{cases} \tilde{\mathbf{r}}_1 = (\Delta_{1,2}\tilde{x}_2,\; \Delta_{1,2}\tilde{y}_2,\; d_1) \\ \tilde{\mathbf{r}}_2 = (\Delta_{2,3}\tilde{x}_3,\; \Delta_{2,3}\tilde{y}_3,\; d_2) \end{cases}$$

con los factores:

$$\Delta_{1,2} = 1 + \sqrt{\frac{\varepsilon_{r1}}{\varepsilon_{r2}}}\left(\frac{d_1 - d_2}{d_2}\right), \qquad \Delta_{2,3} = \frac{1 + \sqrt{\frac{\varepsilon_{r2}}{\varepsilon_{r3}}}\frac{d_2 - d_3}{d_3 - d_1}}{1 - \sqrt{\frac{\varepsilon_{r2}}{\varepsilon_{r3}}}\sqrt{\frac{\varepsilon_{r1}}{\varepsilon_{r3}}}\frac{(d_3-d_2)(d_1-d_2)}{d_2(d_3-d_1)}}$$

**Variables:**
- `d_1`, `d_2`, `d_3` — distancias desde la antena a cada interfaz en el eje z
- `epsilon_{r1}`, `epsilon_{r2}`, `epsilon_{r3}` — permitividades de cada capa

### Algoritmo Phase Shift Migration (PSM) — dominio frecuencia-número de onda

$$\rho_Q(z') = \sum_{n=1}^{N} \mathcal{F}_{xy}^{-1}\left\{ E_{scatt}(k_x, k_y, f_n) \cdot e^{+jd_1 k_{z,1,n}} \cdot \prod_{p=2}^{q-1} e^{+j(d_p - d_{p-1})k_{z,p,n}} \cdot e^{+j(z' - d_{q-1})k_{z,q,n}} \right\}$$

**Variables:**
- `k_x`, `k_y` — componentes del número de onda en x e y
- `k_{z,p,n}` — componente z del número de onda en la p-ésima capa a la n-ésima frecuencia
- `d_p` — profundidad de la p-ésima interfaz medida desde la antena
- `F_{xy}^{-1}` — transformada de Fourier inversa 2D en xy

### Número de onda vertical en la p-ésima capa

$$k_{z,p,n} = \sqrt{\left(\frac{2\pi f_n}{v_p}\right)^2 - k_x^2 - k_y^2}$$

**Variables:**
- `f_n` — n-ésima frecuencia
- `v_p` — mitad de la velocidad de propagación en la p-ésima capa ($v_p = c / (2\sqrt{\varepsilon_{r,p}})$ en configuración monoestática equivalente)
- `k_x`, `k_y` — números de onda transversales

### Ecualización de la respuesta en frecuencia (DAS normalizado)

$$\overline{\rho_Q}(\mathbf{r}', n) = \rho_Q(\mathbf{r}', n) \;/\; \max\{|\rho_Q(\mathbf{r}', n)|\}$$

Las N imágenes normalizadas se suman para obtener la imagen final. Esto neutraliza la mayor amplitud de las bajas frecuencias que enmascara las altas frecuencias.

### Resolución en rango (profundidad)

$$\Delta R = \frac{v_{prop}}{2 \Delta BW}$$

**Variables:**
- `v_prop` — velocidad de propagación en el medio ($= c/\sqrt{\varepsilon_r}$ para suelo sin pérdidas)
- `Delta_BW` — ancho de banda efectivo del sistema

### Celdas de dominio de reconstrucción (particionado por diagrama de antena)

$$\begin{cases} C_{x,n} = \dfrac{2h\tan(\phi_n/2)}{a} \\[6pt] C_{y,n} = \dfrac{2h\tan(\theta_n/2)}{b} \end{cases}$$

**Variables:**
- `h` — altura de las antenas sobre el suelo
- `theta_n`, `phi_n` — anchos de haz de la antena Vivaldi en los planos E y H a la n-ésima frecuencia
- `a`, `b`, `c` — dimensiones de la celda de reconstrucción en x, y, z (iguales a los pasos de muestreo $\delta x$, $\delta y$, $\delta z$)

---

## 3. Suposiciones del Modelo

1. El medio se modela como capas planares horizontales con permitividades constantes por capa; la conductividad introduce pérdidas de propagación dependientes de la frecuencia.
2. La configuración es cuasi-monoestática: Tx y Rx están muy próximos (prácticamente en la misma posición).
3. El algoritmo PSM asume que el campo EM se genera desde cada punto del dominio de investigación (modelo de reflector explosivo, ERM), y que el medio es homogéneo en el plano XY dentro de cada capa.
4. La velocidad de propagación se asume $v_p = c/2$ si la permitividad del suelo no puede estimarse (equivalente a $\varepsilon_r = 1$).
5. Para la ecualización, se considera que la respuesta en frecuencia de las antenas Vivaldi varía más de 7 dB en la banda de trabajo, dominando las frecuencias bajas y reduciendo la resolución efectiva.
6. El particionado del dominio de imagen (imaging domain partitioning) usa el diagrama de radiación de las antenas a cada frecuencia para restringir el número de puntos del dominio a evaluar, reduciendo el costo computacional sin perder calidad.
7. La validación por simulación usa el software de código abierto gprMax (FDTD 2D) con pulso gaussiano derivado, frecuencia central 3 GHz, ancho de banda 3.5 GHz.
8. El suelo de playa se modela con $\varepsilon_r = 3$ (obtenido de caracterización previa).

---

## 4. Notas Adicionales

- **Tres mejoras propuestas:** (1) ecualización de la respuesta en frecuencia de las antenas Tx/Rx; (2) procesado GPR-SAR por sub-bandas (sub-band processing) para combinar las ventajas de penetración de las bajas frecuencias con la resolución de las altas; (3) particionado del dominio de imagen según el ancho de haz de la antena a cada frecuencia.
- **Comparación DAS vs PSM:** PSM es 30–50 veces más rápido que DAS gracias al uso de FFT, pero requiere grilla de adquisición plana y uniforme, y medios con capas horizontales. DAS es más flexible ante geometrías irregulares.
- **Mejora de resolución con ecualización:** la resolución en rango (interfaz aire-arena) mejora de 3.9 cm (sin ecualización) a 2.7 cm (con ecualización) para la banda 1–5 GHz.
- **Sub-bandas:** bandas bajas (1–3 GHz) ofrecen mejor penetración pero menor resolución; bandas altas (4–6 GHz) dan mejor resolución pero menos penetración. No es posible obtener una sola imagen SAR que combine ambas ventajas simultáneamente con este método.
- **Escenarios de simulación:** Escenario A (arena seca, $\varepsilon_r = 2.5$, $\sigma = 0.02$ S/m); Escenario B (arena seca + arena húmeda, $\varepsilon_r = 3.5$); Escenario C (suelo arcilloso de Puerto Rico, modelo de Debye de dos términos, tres niveles de humedad: 2.5%, 5%, 10%).
- **Validación en entorno realista:** playa en coordenadas 43.545, -5.694; objetos enterrados: contenedor plástico vacío (diámetro 23 cm, 8 cm en aire) y disco metálico (diámetro 18 cm, 1 cm de espesor); profundidad de homogenización de la arena: 25 cm.
