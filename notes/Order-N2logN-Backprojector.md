# An Order N²log(N) Backprojector Algorithm for Focusing Wide-Angle Wide-Bandwidth Arbitrary-Motion SAR

**Autores:** John McCorkle, Martin Rofheart  
**Año:** 1996  
**Fuente/Publicación:** Proceedings of SPIE, Radar Sensor Technology, Vol. 2747, pp. 25–36. U.S. Army Research Laboratory / Soft Machine Resources.

---

## 1. Geometría del Sistema

Sistema SAR monoestático de **ultra-banda ancha (UWB)** con **trayectoria arbitraria**. Diseñado específicamente para sensores como el BoomSAR del ARL (apertura 100 m, muestreo 0.1 m, frecuencia 8 GHz). Modo **spotlight** implícito (apertura completa para cada píxel). La imagen es 2D ($N \times N$ píxeles). Aplicable también a GPR-SAR con cambio de velocidad de propagación.

El algoritmo es una **descomposición quadtree** del proceso de backprojection: en cada etapa se divide la imagen en 4 subimágenes y se reduce el número de proyecciones a la mitad.

---

## 2. Ecuaciones de Resolución SAR

### Pixel enfocado (Delay-and-Sum)

$$f_{q,r}(t) = \sum_{j=0}^{L-1} z_j \cdot s_j(t + T_{q,r,j})$$

**Variables:**
- `(q,r)` — índice del píxel en la imagen
- `j` — índice de la posición en la apertura ($j = 0 \ldots L-1$)
- `z_j` — peso de apertura (ventana de Hamming, etc.)
- `s_j(t)` — señal recibida con corrección de pérdida $R^2$ en la posición $j$
- `T_{q,r,j}` — desplazamiento temporal = round-trip time desde sensor $j$ al píxel $(q,r)$

### Número de operaciones: quadtree con estrategia halving

$$OPS = \mathcal{O}\!\left(\sum_{m=1}^{M} 2^{2m} \cdot 2^{M-m} \cdot 2 \cdot 2^{M-m} \cdot \sqrt{2}\right) = \mathcal{O}\!\left(N^2 \cdot \log N\right)$$

donde $N = 2^M$, imagen $N \times N$.

### Número total de subimágenes en la descomposición

$$\sum_{m=1}^{M} \alpha^{2m} = \frac{\alpha^2(\alpha^{2M}-1)}{\alpha^2-1}\bigg|_{\alpha=2} = \frac{2^{2(M+1)}-4}{3}$$

**Variables:**
- `M` — número de etapas de recursión, $N = 2^M$
- `\alpha = 2` — factor de partición (imagen dividida en cuadrantes)

### Centro de fase de la sub-apertura (caso simple)

$$x_{c,k} = \sum_{j \in J_{p,k}} w_{p,j,k} \cdot x_{p,j}, \quad k = 0 \ldots N_c - 1$$

donde $w_{p,j,k} = 1/L_c$ para la media simple.

**Variables:**
- `x_{c,k}` — centro de fase del k-ésimo elemento del nodo hijo $c$
- `J_{p,k}` — conjunto de índices del padre que forman la sub-apertura $k$
- `L_c` — longitud de la sub-apertura (elementos del padre sumados)

### Pesos de centros de fase desplazados (parámetro $\eta$)

$$w_{p,k,j} = \begin{cases} \frac{k\cdot\eta}{N_c-1} + \frac{1-\eta}{2} & j = \min(J_{p,k}) \\ 0 & \text{otro} \\ 1 - \frac{k\cdot\eta}{N_c-1} - \frac{1-\eta}{2} & j = \max(J_{p,k}) \end{cases}$$

**Variables:**
- `\eta` — parámetro de desplazamiento de centros (0 = media, 1 = extremos)

### Rango al primer sample del hijo

$$r_{c,k} = R_{c,k} - H_{c,k} \alpha_c$$

**Variables:**
- `R_{c,k}` — distancia del centro de fase hijo $x_{c,k}$ al centro de la subimagen $c$
- `H_{c,k}` — número de bins de rango antes del centro de la subimagen
- `\alpha_c` — espaciado de muestras de rango (m)

### Distancia padre→hijo (interpolación de rango)

$$\Gamma_{c,k,j,m} = \sqrt{\Psi_{c,k,m}^2 + L_{p,k,j}^2 - 2 L_{p,k,j} \Psi_{c,k,m} \cos\varphi_{c,k,j,m}}$$

**Variables:**
- `\Psi_{c,k,m} = r_{c,k} + m\alpha_c` — rango desde centro hijo al bin $m$
- `L_{p,k,j}` — distancia del elemento $j$ del padre al centro de fase $x_{c,k}$
- `\varphi_{c,k,j,m}` — ángulo entre la línea a $x_{c,k}$ y la línea al bin $m$

### Recursión: datos del hijo a partir del padre

$$s_{c,k}(m) = \sum_{j \in J_{p,k}} s_{p,j}\!\left(N(c,k,j,m)\right), \quad m=0\ldots M_{c,k}-1$$

donde $N(c,k,j,m) = (\Gamma_{c,k,j,m} - r_{p,j})/\alpha_c$ es el índice de punto flotante (requiere interpolación).

---

## 3. Suposiciones del Modelo

1. Trayectoria de vuelo completamente arbitraria; no se requiere movimiento uniforme.
2. Compensación de movimiento y variación de velocidad de propagación se manejan como simples desplazamientos de índice.
3. El algoritmo acepta señales UWB con ángulos de integración amplios (wide-angle).
4. Los artefactos son **locales** (confinados alrededor del blanco que los genera), a diferencia de los artefactos globales de los métodos FFT.
5. La estrategia halving/doubling: en cada etapa se dobla la longitud de la sub-apertura ($L_c = 2$) y se divide la imagen en 4 cuadrantes ($D_{p,a} = D_{p,r} = 2$).
6. La interpolación de tiempo introduce pérdidas; con espaciado de píxeles de 0.0625 m, la pérdida respecto al ideal es ~2.2 dB (quadtree) vs ~1.8 dB (delay-and-sum).
7. La imagen resultante puede tener artefactos entre lóbulos, controlables ajustando el factor holdoff $n$.

---

## 4. Notas Adicionales

- El factor holdoff $n$ controla el balance velocidad/calidad: el set completo de proyecciones se usa en las primeras $n$ descomposiciones.
- Validado con datos del BoomSAR (ARL): apertura 100 m, resolución cruzada medida 1.0 m (quadtree) vs 0.9125 m (delay-and-sum).
- Complejidad comparable al algoritmo $\omega k$, con la ventaja de aceptar trayectorias arbitrarias y propagación no uniforme.
- La implementación es nativamente paralelizable (árbol de fácil descomposición, comunicación interprocessor mínima).
- Referencia clave para los algoritmos FFBP posteriores (Ulander 2003, Góes 2020).
