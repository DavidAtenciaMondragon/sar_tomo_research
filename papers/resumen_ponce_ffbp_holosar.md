# Resumen crítico: Ponce et al. (IGARSS 2011 y TGRS 2016) desde la óptica de optimización de vuelo biestático

Ambos artículos son de Octavio Ponce (DLR) y tratan **SAR circular monoestático**, no biestático. Ninguno optimiza una trayectoria de Tx/Rx separados. Su valor para tu investigación está en dos bloques reutilizables: (1) un algoritmo de enfoque en dominio tiempo eficiente y preciso para trayectorias no rectilíneas (aplicable a tus trayectorias espirales Tx/Rx), y (2) un marco cuantitativo de "línea de base crítica" y compromiso resolución-PSLR para diseñar la separación entre pasadas/offsets, que es directamente trasladable a tu problema de elegir el offset del receptor (tus archivos `radarRx_espiral_offset0.json` / `offset180.json`).

---

## 1. Ponce et al., "Processing of Circular SAR Trajectories with Fast Factorized Back-Projection", IGARSS 2011

**Qué resuelve:** el Fast Factorized Back-Projection (FFBP) clásico de Ulander asume tramo de vuelo rectilíneo. Este trabajo lo adapta a trayectorias circulares (CSAR), aproximando el círculo como un polígono de segmentos rectos y definiendo, para cada subapertura, un sistema de coordenadas local tangente rotado un ángulo β respecto al sistema global. Como la ley de cosenos usada en el FFBP lineal ya no es válida (el eje x' local cambia con cada subapertura), sustituyen la interpolación por una transformación lineal + rotación (Ec. 1) y luego una interpolación polar-a-polar (P2P) usando producto punto para expresar cada subimagen en las coordenadas polares de la subapertura vecina antes de sumarlas coherentemente.

**Puntos cuantitativos relevantes:**
- Complejidad ~log₂N en vez de O(N) de la BP directa. Speedup medido: **~500× en CPU, ~1800× en GPU** frente a BP directa, para imágenes hasta 25K×25K píxeles a resolución λ/4.
- Requisito de muestreo angular (Ec. 3): la diferencia de distancia radial entre dos muestras angulares consecutivas Δr = R(α+Δα) − R(α) no puede superar c₀/(2B) sin causar aliasing; esto fija el paso angular máximo en función del tamaño de subapertura Δd y la distancia radial r.
- Hay un tamaño de subapertura óptimo (L_max): pasado cierto punto la curvatura de la subapertura curva el espectro y conviene pasar a interpolación polar-a-cartesiana en vez de seguir interpolando en polar — es una compensación entre exactitud numérica y eficiencia, no gratis.
- Validado con datos reales E-SAR (L-band-like, polarimetría completa) sobre Kaufbeuren: coherencia de amplitud media 0.99991, desviación de fase 4.7°, y tiempos de procesamiento de ~11 h (CPU) / ~3 h (GPU) por canal polarimétrico frente a ~121 días (CPU) / ~12 días (GPU) de la BP directa.

**Relevancia para tu trabajo:** tu pipeline biestático (Snell + TM multicapa) ya usa un enfoque tipo back-projection por trayectoria no ideal. Este paper es la referencia técnica de cómo acelerar ese tipo de enfoque cuando la trayectoria del sensor (aquí, Tx o Rx en espiral) no es rectilínea, sin perder la ventaja de la BP de manejar topografía real y trayectorias arbitrarias — algo que un procesador en dominio de Fourier no puede hacer. El criterio de muestreo angular (Ec. 3) es directamente aplicable para decidir cuántos puntos de la espiral necesitas muestrear para evitar aliasing radial en tu geometría.

---

## 2. Ponce et al., "First Airborne Demonstration of Holographic SAR Tomography (HoloSAR)", IEEE TGRS, vol. 54(10), 2016

**Concepto central:** HoloSAR = apertura sintética circular (azimut, 360°) **+** una segunda apertura sintética en elevación, formada por **múltiples pasadas circulares a distintas alturas** (líneas de base verticales ΔB_L). Es formalmente el análogo monoestático de lo que tú haces al comparar offsets de receptor (`offset0` vs `offset180`, o al plantear un `otimiza_tx_vs_planvoo`): varias pasadas separadas espacialmente que se combinan para ganar resolución en la dirección perpendicular a la línea de vista (LOS⊥).

### 2.1 Teorema de corte-proyección 3D y por qué importa la geometría de las pasadas
Cada pasada circular observa una porción del espectro 3D que es un **anillo** en (kₓ,k_y) inclinado según el ángulo de elevación θ_el (Fig. 4, Ecs. 3–6). Una sola pasada deja huecos espectrales en la dirección k_z → produce lóbulos secundarios en forma de cono en (x,y,z) (visible en Fig. 21a). Añadir pasadas a distintas alturas rellena esos huecos y reduce los lóbulos — pero **solo si la separación entre pasadas está bien elegida**.

### 2.2 Línea de base crítica (el resultado más trasladable a tu problema)
Definen la línea de base crítica en elevación como el límite de separación entre pasadas contiguas que **no** deja huecos espectrales (Ec. 24):

```
B⊥,crit = min( (λ/c)·BW·r·tan(θ_el − α_in),  (λ/c)·BW·r·sin(θ_el − α_in) )
```

donde BW es el ancho de banda del sistema, r la distancia oblicua, θ_el el ángulo de elevación y α_in la pendiente del terreno. Recomendación práctica del paper: la separación entre pasadas contiguas ΔB_L·cos(θ_el − α_tilt) debe acercarse a B⊥,crit **sin superarla** — superarla causa aliasing (huecos), quedarse muy por debajo desperdicia vuelo sin ganancia proporcional.

Con esto derivan la resolución vertical y el PSLR en función del ancho de banda y del número de pasadas (Ecs. 19–23, Fig. 6): p. ej. con 400 MHz y 4 pasadas obtienen δ_z ≈ 0.0846 m con PSLR ≈ −20 dB, frente a solo −8 dB con 1 pasada. La mejora es **no lineal** con el número de pasadas.

### 2.3 Validación experimental (campaña real, L-band, DLR F-SAR, Kaufbeuren)
- 19 pasadas circulares, ΔB_L medio de 15 m, línea de base vertical máxima L_H = 270 m (línea de base normal ΔB_L,⊥ = 223.72 m) — muy por debajo de la línea de base crítica calculada (B⊥,crit ≈ 143–256 m según dirección), es decir, **deliberadamente sobremuestreada** en elevación para maximizar supresión de lóbulos a costa de tiempo de vuelo.
- Un experimento de control muy útil: comparando 3 pasadas espaciadas cerca de la crítica (ΔB_L=135 m) contra 3 pasadas con separación irregular (ΔB_L=[45,195] m) contra 19 pasadas regulares, muestran que superar la línea de base crítica introduce artefactos de aliasing claramente visibles en el IRF 3D (Fig. 21b vs 21c vs 21d).
- Dato práctico interesante: **solo 5 pasadas bien distribuidas** a lo largo de la apertura vertical total (L_H=270 m, es decir, cerca de la separación crítica) dieron una calidad de imagen 2D (x,y) comparable a las 19 pasadas densamente espaciadas — es decir, para resolución horizontal pura, pocas pasadas bien separadas ≈ muchas pasadas juntas; pero para suprimir bien la ambigüedad/lóbulos en la dirección LOS⊥ si hace falta el muestreo denso.

### 2.4 Errores de navegación y autoenfoque (PGA) — relevante si algún día vuelas de verdad
La precisión típica de navegación aerotransportada (1–10 cm) es insuficiente para enfocar a λ/4 en L-band sobre 360°. Un resultado no obvio: en geometría **circular** (a diferencia de stripmap lineal), un error constante en la posición horizontal del track se manifiesta como una **señal sinusoidal** en la historia de rango (no una rampa lineal), mientras que un error constante en altura se ve como un **desplazamiento constante** (Ecs. 29–33, Fig. 15). Esto da una firma diagnóstica distinta para clasificar errores de navegación en trayectorias curvas/espirales.

### 2.5 Estrategias de procesamiento (coherente / incoherente / GLRT)
Para dispersores distribuidos (bosque), el paper compara sumar subaperturas coherentemente (mejor resolución en (x,y) pero lóbulos altos en blancos fuertes tipo bordes de bosque), incoherentemente (mejor supresión de lóbulos cónicos, peor resolución), y con GLRT + compressive sensing en elevación (mejor para blancos muy directivos, mal para dispersión de volumen semitransparente). Concluyen que el mejor supuesto para dispersores distribuidos es tratarlos como **isótropos**.

---

## Implicaciones concretas para tu optimización de vuelo biestático

1. **Sustituto directo de "elegir offsets a ojo":** en vez de fijar `offset0`/`offset180` de forma ad hoc, la fórmula de línea de base crítica (Ec. 24) —adaptada a tu geometría biestática, con r y θ definidos por la trayectoria Tx/Rx real y no por una única distancia oblicua monoestática— te da un criterio cuantitativo de cuánta separación entre pasadas (o entre offsets del receptor) necesitas para ganar resolución en profundidad sin caer en aliasing espectral, y cuántas pasadas necesitas realmente para un PSLR objetivo.
2. **FFBP circular (paper 1)** es el candidato natural si en algún momento el pipeline Snell/TM multicapa se vuelve el cuello de botella computacional para trayectorias espirales de alta resolución — mantiene la ventaja de la BP directa (topografía real, trayectoria no ideal) con complejidad log₂N.
3. **Ninguno de los dos aborda biestática real** (Tx y Rx en trayectorias distintas, refracción Snell en interfases): el teorema de corte-proyección y la línea de base crítica están derivados para geometría monoestática (r común de ida y vuelta). Si quieres usar estos resultados con rigor, la extensión pendiente es rederivar B⊥,crit y las Ecs. 19–23 para el caso r_Tx ≠ r_Rx con el path length compuesto (como ya haces en `precomputeTxDataMultislab.m` con `sum(Rseg)`), en vez de asumir que aplican literalmente.
4. **Diagnóstico de errores de trayectoria (§2.4)** es directamente reusable si en el futuro comparas tu plan de vuelo optimizado contra una trayectoria realmente volada: el patrón sinusoidal vs. desplazamiento constante te dice si el error es de posición horizontal o de altura sin necesitar inversión completa.
