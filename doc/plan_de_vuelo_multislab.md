# Optimización del Plan de Vuelo en Escenario Multicapa

## 1. Motivación

El script `run_plano_de_voo.m` optimiza la posición del receptor SAR bistático para maximizar la potencia recibida y/o la resolución 3D de targets subsuperficiales, asumiendo un medio subterráneo **homogéneo** (un único índice de refracción $n_2$). Esta hipótesis es insuficiente para suelos estratificados, donde cada capa tiene su propio índice de refracción y espesor.

`run_plano_de_voo_multislab.m` generaliza ese escenario a **N capas planas y horizontales**, describiendo el trayecto del rayo y los coeficientes TM mediante el formalismo matricial ABCD documentado en `models/multiple_slab_tm.md`.

---

## 2. Escenario simulado

### 2.1 Geometría

| Zona          | Índice de refracción | Espesor     |
|---------------|----------------------|-------------|
| Aire          | $n_0 = 1$            | semi-infinito (arriba) |
| Capa 1 (suelo superficial) | $n_1 = 2$ | $d_1 = 1\ \text{m}$ |
| Capa 2 (suelo profundo)    | $n_2 = 3$ | $d_2 = 1\ \text{m}$ |

- **Transmisor (Tx):** hélice cónica, $z \in [80, 120]\ \text{m}$, radio $\in [147.5, 172.5]\ \text{m}$, velocidad $120\ \text{m/s}$, 2 vueltas.
- **Targets:** volumen cúbico en capa 2 ($z \in [-1.1, -1.9]\ \text{m}$), centrado en el origen horizontal, $2 \times 2 \times 2 = 8$ puntos.
- **Receptor (Rx):** posición horizontal optimizada en cada instante; altitud = altitud instantánea del Tx.

### 2.2 Frecuencia

Se utiliza **L-band** ($f_c = 400\ \text{MHz}$, $\lambda \approx 0.75\ \text{m}$) para asegurar penetración subsuperficial razonable con medios dieléctricos ideales (sin pérdidas). El ancho de banda del chirp es $B = 50\ \text{MHz}$.

> **Nota física:** a 400 MHz, la penetración en suelo con pérdidas reales puede ser de varios metros; en suelo seco (desierto) es habitual. El modelo aquí es lossless (n puramente real), lo que constituye una premisa conservadora-favorable para la transmitancia.

---

## 3. Extensiones al modelo de propagación

### 3.1 Trayecto del rayo: invariante de Snell multislab

Para un rayo descendente desde el radar hasta el target, el invariante de Snell generalizado es:

$$p = n_0 \sin\theta_0 = n_1 \sin\theta_1 = n_2 \sin\theta_2 = \cdots = \text{cte}$$

donde $\theta_i$ es el ángulo de cada tramo con la normal a las interfaces (horizontales). La distancia horizontal total recorrida por el rayo es:

$$U(p) = \sum_{i=0}^{L} h_i \tan\theta_i(p)$$

con $h_i$ la altura/espesor del tramo $i$ (en aire: $h_0 = z_{\text{radar}}$; en cada capa completa: $h_i = d_i$; en la capa del target: $h_L$ = fracción de la capa hasta la profundidad del target).

Para que el trayecto conecte exactamente el radar $P$ con el target $T$, se impone $U(p) = \|\mathbf{P}_{xy} - \mathbf{T}_{xy}\|_2$. El parámetro de rayo $p$ se obtiene por **bisección** sobre $p \in [0,\ \min_i(n_i) - \varepsilon)$ (implementada en `common/calculateSlantRangeMultislab.m`).

La longitud física de cada segmento es:

$$R_i = \frac{h_i}{\cos\theta_i}$$

y la longitud total del trayecto (usada para el factor de dispersión esférica en la ecuación de radar) es:

$$R_{\text{total}} = \sum_{i=0}^{L} R_i$$

### 3.2 Coeficiente TM: formalismo ABCD multicapa

Los coeficientes de reflexión y transmisión para polarización TM a través de la pila de capas se calculan mediante la matriz de transferencia ABCD (ver `models/multiple_slab_tm.md`).

Para la cascada de $M$ slabs finitos entre el medio de entrada $n_{\text{in}}$ y el medio de salida $n_{\text{out}}$:

$$\mathbf{M}_{\text{tot}} = \prod_{j=1}^{M} \mathbf{M}_j, \qquad \mathbf{M}_j = \begin{pmatrix} \cos(k_{z,j} d_j) & jZ_j \sin(k_{z,j} d_j) \\ j/Z_j \cdot \sin(k_{z,j} d_j) & \cos(k_{z,j} d_j) \end{pmatrix}$$

donde $k_{z,j} = k_0 n_j \cos\theta_j$ y $Z_j = (\eta_0/n_j)\cos\theta_j$ es la impedancia TM del slab $j$.

El coeficiente de reflexión complejo de la pila completa:

$$\Gamma_{\text{TM}} = \frac{A Z_{\text{out}} + B - Z_{\text{in}}(C Z_{\text{out}} + D)}{A Z_{\text{out}} + B + Z_{\text{in}}(C Z_{\text{out}} + D)}$$

La transmitancia de potencia extremo a extremo:

$$T_{\text{TM}} = |\tau_{\text{TM}}|^2 \frac{Z_{\text{in}}}{Z_{\text{out}}}$$

que satisface $R_{\text{TM}} + T_{\text{TM}} = 1$ para medios sin pérdidas.

**Implementación:** `common/multislabTMcoef.m`. Se llama con:
```
n_in    = 1  (aire)
n_slabs = n_layers(1:L-1)   % capas intermedias entre aire y capa del target
d_slabs = d_layers(1:L-1)
n_out   = n_layers(L)       % capa que contiene al target (semi-infinita en el modelo)
```

Si $L = 1$ (target en la primera capa), `n_slabs = []` y la función reduce exactamente a las ecuaciones de Fresnel para una única interfaz.

---

## 4. Ángulo óptimo de incidencia (multicapa)

Para una única interfaz aire→$n_2$, la transmitancia TM es máxima en el **ángulo de Brewster**:

$$\theta_B = \arctan\!\left(\frac{n_2}{n_1}\right)$$

donde $|r_{\text{TM}}| = 0$ y $T_{\text{TM}} = 1$. Esta solución cerrada no existe en general para una pila multicapa (los efectos de interferencia tipo Fabry-Pérot entre capas desplazan el máximo de $T_{\text{TM}}$).

En el escenario simulado (n=2, d=1m sobre n=3):
- Brewster de la interfaz 1 (aire→n=2): $\theta_{B,1} = \arctan(2) \approx 63.4°$
- Brewster de la interfaz 2 (n=2→n=3) en términos del ángulo en n=2: $\theta_{B,2} = \arctan(3/2) \approx 56.3°$, que requeriría $\sin\theta_0 = 2\sin(56.3°) \approx 1.66 > 1$ (imposible desde el aire)

Por tanto, para ángulos de incidencia en aire alcanzables, **la segunda interfaz nunca puede alcanzar su Brewster propio**. El ángulo óptimo de la pila completa se encuentra mediante un barrido numérico:

```matlab
theta_sweep = linspace(0, pi/2 - 1e-4, 2000);
[~,~,~,T_sweep] = multislabTMcoef(f, theta_sweep, 1, n_slabs, d_slabs, n_layers(L));
[~, idx_opt] = max(T_sweep);
angulo_opt = theta_sweep(idx_opt);
```

Este ángulo óptimo reemplaza al ángulo de Brewster en la función `buildRxSearchGeometry` para guiar el arranque de la optimización.

---

## 5. Ecuación de radar bistática multicapa

La potencia recibida por el target $i$ es:

$$P_{r,i} = \frac{P_t \, G_{t,i} \, G_{r,i} \, \sigma \, T_{1,i} \, T_{2,i} \, \lambda^2}{(4\pi)^3 \, R_{T,i}^2 \, R_{R,i}^2}$$

donde:

| Símbolo | Descripción |
|---------|-------------|
| $P_t$ | Potencia transmitida (W) |
| $G_{t,i}, G_{r,i}$ | Ganancias de antena Tx y Rx hacia target $i$ |
| $\sigma$ | Sección eficaz del target (m²) |
| $T_{1,i}$ | Transmitancia TM multicapa del trayecto Tx→target$_i$ |
| $T_{2,i}$ | Transmitancia TM multicapa del trayecto target$_i$→Rx |
| $R_{T,i}$ | Longitud física total del trayecto Tx→target$_i$ (todos los segmentos) |
| $R_{R,i}$ | Longitud física total del trayecto Rx→target$_i$ |
| $\lambda$ | Longitud de onda |

**Nota sobre $T_1$ y $T_2$:** La transmitancia de retorno $T_{2,i}$ se evalúa llamando a `multislabTMcoef` con el ángulo de incidencia en aire del trayecto Rx→target, utilizando la misma convención (aire→capas→capa L). Para medios sin pérdidas, la reciprocidad de Lorentz garantiza que $T_{\text{up}} = T_{\text{down}}$ para el mismo parámetro de rayo $p$.

La potencia total sumada sobre todos los targets del volumen es:

$$P_r = \sum_i P_{r,i}$$

---

## 6. Función de costo combinada

Se minimiza la misma función adimensional que en el caso monocapa:

$$J = (1 - \alpha_{\text{res}}) \cdot \left[-\log_{10}(P_r)\right] + \alpha_{\text{res}} \cdot \log_{10}(\delta_{xy} \cdot \delta_z)$$

- $\alpha_{\text{res}} = 0$: solo maximiza potencia.
- $\alpha_{\text{res}} = 1$: solo minimiza resolución.
- Valor por defecto: $\alpha_{\text{res}} = 0.3$.

Las métricas de resolución $\delta_{xy}$ y $\delta_z$ se calculan con `calculateBistaticResolution`, usando $n_{\text{eff}} = n_{\text{layers}}(L_{\text{centroide}})$ como índice de refracción del medio del target.

---

## 7. Premisas del modelo

1. **Interfaces planas y horizontales.** El modelo de Snell y ABCD asume planaridad. No se consideran rugosidad ni heterogeneidad lateral.

2. **Medios sin pérdidas (dieléctricos ideales).** Los índices de refracción son puramente reales. En suelos reales, la parte imaginaria de $\varepsilon_r$ introduce atenuación adicional que reduciría la potencia recibida.

3. **Polarización TM (p-polarization).** Se optimiza para TM porque el ángulo de Brewster (nulo de reflexión) existe solo para esta polarización. Para TE no hay ángulo de Brewster, y la reflectancia crece monótonamente con el ángulo de incidencia.

4. **Target en la capa $L$, campo en el medio $n_L$ semi-infinito.** El target se modela como dispersor puntual dentro de la capa $L$. El trayecto Tx→target solo incluye las $L-1$ capas completas entre la superficie y la capa del target; el tramo parcial dentro de la capa $L$ contribuye a $R_{\text{total}}$ pero no a los coeficientes TM (que se calculan para interfaces planas entre slabs completos).

5. **Ecuación de radar en espacio libre para la propagación.** La dispersión esférica sigue la ley $1/R^2$ para el factor de propagación de potencia, con $R$ la longitud física del trayecto refractado. No se modela propagación en guía de onda ni efectos de Tierra curva.

6. **Reciprocidad del camino de vuelta.** Para medios sin pérdidas y polarización fija, la transmitancia de ida ($T_1$, aire→capa L) es igual a la transmitancia de vuelta ($T_2$, capa L→aire) para el mismo ángulo de incidencia en aire. Esta propiedad se verifica a partir de la simetría de la matriz ABCD y la conservación de energía.

7. **Bisección sobre el parámetro de rayo.** La tolerancia de bisección durante la optimización es $\varepsilon_{\text{opt}} = 10^{-4}\ \text{m}$ (vs. $10^{-10}\ \text{m}$ en el pipeline GS/PROC), suficiente para la búsqueda de posición del receptor con resolución de decenas de metros.

---

## 8. Diferencias con `run_plano_de_voo.m`

| Aspecto | Monocapa (`run_plano_de_voo.m`) | Multicapa (este script) |
|---------|--------------------------------|-------------------------|
| Parámetros del suelo | `IndiceRefracaoSolo` (escalar) | `CapasSuelo[{n,d}, ...]` (vector) |
| Solucionador de rayos | `calculateRefractionPointFermat` (Newton-Raphson, 1 interfaz) | `calculateSlantRangeMultislab` (bisección, N interfaces) |
| Coeficientes TM | `calculateTMcoef` (Fresnel, 1 interfaz) | `multislabTMcoef` (ABCD, N interfaces) |
| Precómputo Tx | `precomputeTxData` | `precomputeTxDataMultislab` |
| Función objetivo | `objective_function` | `objective_function_multislab` |
| Ángulo óptimo | $\arctan(n_2/n_1)$ (Brewster analítico) | barrido numérico de $T_{\text{TM}}(\theta)$ |
| Índice para resolución | `n2` (único) | `n_layers(L_centroid)` (capa del centroide) |
| Addpath necesario | `gs, proc, tools, flightpath` | + `common` |

---

## 9. Archivos del escenario

| Tipo | Archivo |
|------|---------|
| Parámetros sistema | `parametros/system_multislab_plano_voo.json` |
| Parámetros Tx | `parametros/radarTx_multislab_plano_voo.json` |
| Parámetros targets | `parametros/target_multislab_plano_voo.json` |
| Script principal | `simulations/run_plano_de_voo_multislab.m` |
| Precómputo Tx | `flightpath/precomputeTxDataMultislab.m` |
| Función objetivo | `flightpath/objective_function_multislab.m` |
| Solucionador de rayos | `common/calculateSlantRangeMultislab.m` |
| Coeficientes TM | `common/multislabTMcoef.m` |
| Salida de datos | `io/plan_vuelo_multislab/` |

---

## 10. Resultados esperados

Para el escenario con dos capas (n=2/n=3) y targets en capa 2 ($z \approx -1.5\ \text{m}$):

- **Ángulo óptimo multicapa** $\approx 63°$ (próximo al Brewster de la primera interfaz, pero desplazado hacia ángulos menores por el efecto de interferencia de la segunda capa).
- **Posición óptima del Rx:** horizontalmente en dirección opuesta al Tx respecto al centroide del volumen, a distancia $\approx (z_{\text{Rx}} + d_{\text{total}}) / \tan(\theta_{\text{opt}})$ del centroide.
- **Potencia recibida:** inferior al caso monocapa debido a las pérdidas adicionales de reflexión en la segunda interfaz (energía parcialmente reflejada en la interfaz capa1/capa2).
- **Resolución 3D:** similar al caso monocapa, ya que el modelo de resolución usa el índice efectivo de la capa del target ($n_{\text{eff}} = 3$) e ignora las capas superiores (aproximación válida cuando los ángulos en la capa del target son pequeños).
