# Sincronismo BiSAR por camino directo sin referencia externa

**Análisis de viabilidad y validación por simulación (31/10/2026)**

Escenario: el Radar 1 (Tx/Rx monostático en dron) emite chirps LFM; el Radar 2
es un receptor pasivo que graba, en ventanas discontinuas, el camino directo
Tx→Rx y los ecos bistáticos. No hay GPS, PPS, reloj común ni servidor de
tiempo. La pregunta no es solo si se pueden sincronizar, sino si la
sincronización basta para que el Backprojection (BP) bistático enfoque.

Parámetros de la simulación (heredados de `analise_1610/params/common.json`):
$f_c = 425$ MHz ($\lambda = 0.705$ m), $B = 50$ MHz, $T_p = 1.3\ \mu$s,
PRF = 400 Hz, Tx a 10 m/s, $y=-80$ m, $z=100$ m, apertura de 60 m (6 s);
Rx2 fijo en $(40,-70,20)$ m; cinco reflectores puntuales separados ≥ 32 m.

---

## 1. Modelo y parámetros a estimar

### 1.1 Relojes

El reloj del Radar 1 se toma como referencia de tiempo $t$. Su error absoluto
frente al tiempo físico solo escala $f_c$ y la PRF en unas pocas ppm, y no
afecta al problema relativo. El Radar 2 tiene un tiempo local

$$t_2 = t + \delta(t),\qquad
\delta(t) = \delta_0 + \varepsilon\,t + \tfrac12\gamma\,t^2 + w(t) + j_k ,$$

donde

| Símbolo | Significado | Valores típicos |
|---|---|---|
| $\delta_0$ | offset inicial (arranque en instantes distintos) | arbitrario (s) |
| $\varepsilon$ | error de frecuencia fraccional | OCXO 0.05 ppm · TCXO 1 ppm · XO 20 ppm |
| $\gamma$ | deriva de frecuencia (térmica/envejecimiento) | $10^{-12}$–$10^{-8}$ s⁻¹ |
| $w(t)$ | ruido estocástico del oscilador (FM blanca → paseo aleatorio, $\operatorname{Var}[w(t+\tau)-w(t)] = \sigma_y^2(1\,\mathrm{s})\,\tau$) | $\sigma_y(1\,\mathrm{s})$: $10^{-11}$–$10^{-8}$ |
| $j_k$ | jitter blanco de disparo/muestreo, común a la ventana $k$ | ps |

Si el LO y el reloj de muestreo del Rx2 salen del mismo oscilador (lo
habitual), la fase de demodulación es

$$\varphi_2(t) = 2\pi f_c\,\delta(t) + \phi_0 + \phi_{LO,k}.$$

El **offset de fase** $\phi_0$ es una constante desconocida. El **offset de
frecuencia** de portadora es $\Delta f = \varepsilon f_c$ (425 Hz con 1 ppm).
La **deriva temporal** del muestreo es $\varepsilon t$ (6 µs en 6 s con 1 ppm).
Las tres cosas que menciona el planteamiento son un único proceso
$\delta(t)$ visto desde la envolvente (retardo) y desde la portadora (fase).

### 1.2 Señales en el Rx2 (banda base, tiempo local)

Para el pulso $n$, emitido en $t_n$, y un camino de retardo geométrico $\tau$:

$$y(t_2) = A\,s\!\big(t_2 - \underbrace{[t_n + \tau + \delta]}_{a}\big)\,
\exp\!\big\{-j\,[2\pi f_c(\tau + \delta) + \phi_0 + \phi_{LO}]\big\}.$$

* Camino directo: $\tau_d(n) = \lVert \mathbf p_{Tx}(t_n) - \mathbf p_{Rx}\rVert / c$.
* Eco del píxel $\mathbf g$: $\tau_b(n,\mathbf g) = (\lVert \mathbf p_{Tx}-\mathbf g\rVert + \lVert \mathbf g-\mathbf p_{Rx}\rVert)/c$.

El BP bistático necesita, para cada pulso usado, la hora local
$t_n + \tau_b + \delta$ y la fase $2\pi f_c\delta + \phi_0$.

### 1.3 Vector de parámetros

| Parámetro | Tipo | Necesario para |
|---|---|---|
| $k$: índice del primer pulso observado (asociación) | entero | saber qué $t_n$ y qué $\mathbf p_{Tx}(t_n)$ corresponden a cada ventana |
| $\delta_0$ | continuo | alinear la envolvente |
| $\varepsilon$, $\gamma$ | continuos | extrapolar $\delta$ entre ventanas (puertas de grabación, modelo global) |
| $w(t_j)$, $j_k$ en cada ventana observada | estocásticos | coherencia de fase pulso a pulso |
| $\phi_0$ (+ $\phi_{LO,k}$) | fase | enfoque (la parte variable) y fusión coherente mono+bi (la constante) |
| Molestias: $A_d$, error de navegación del Tx, posición del Rx | — | se estiman o se acotan |

Lo que hay que intercambiar después del vuelo es solo el **log del Tx**:
instantes de emisión $t_n$ y navegación $\mathbf p_{Tx}(t_n)$. El Rx2 graba
ventanas en bruto con marcas de su propio reloj. Todo el sincronismo se hace
en post-proceso.

---

## 2. Observabilidad

El camino directo da, en cada ventana observada $j$ (pulso $n_j$):

$$\hat a_j = t_{n_j} + \tau_d(n_j) + \delta(t_j) + \nu_j,\qquad
\hat\varphi_j = -2\pi f_c\,[\tau_d(n_j) + \delta(t_j)] - \phi_0 - \phi_{LO,j} + \eta_j .$$

**Condición 1: geometría Tx–Rx conocida.** $\tau_d$ y $\delta$ entran *sumados*.
Un error de distancia Tx–Rx es indistinguible de un error de reloj (no
observabilidad fundamental). Para la fase hace falta conocer la distancia
Tx–Rx con error $\ll \lambda$ a lo largo de la línea Tx→Rx. Con referenciado
por el directo, el error que queda en el eco del blanco $P$ es
$k(\mathbf u_{TxP}-\mathbf u_{TxRx})\cdot\mathbf e$ en lugar de
$k\,\mathbf u_{TxP}\cdot\mathbf e$, así que **parte del error de navegación
se compensa**. La sección 6.6 lo cuantifica.

**Condición 2: índice de pulso (asociación).** Una sola ventana solo determina
$\delta$ módulo PRT. Desplazar la hipótesis $k$ pulsos cambia $t_{n}+\tau_d(n)$ en
$k\,\mathrm{PRT} + (s_{n+k}-s_n) + [\tau_d(n+k)-\tau_d(n)]$. Con PRF uniforme
($s_n=0$) y geometría casi lineal en el tiempo, ese cambio es casi afín en $t$
y lo absorben $\delta_0$ y $\varepsilon$: **$k$ es inobservable**. La condición
de observabilidad es que la secuencia de emisión tenga una firma no afín. Para
eso sirve un *stagger* (o un código por pulso) que el Tx registra en su log:
$t_n = n\,\mathrm{PRT} + s_n$, con $s_n$ pseudoaleatorio. Separación entre la
hipótesis correcta y las demás:
$\mathrm{rms}(s_{n+k}-s_n)\approx S\sqrt{2/3}\gg\sigma_t$, con al menos
$J\ge 3$ ventanas (el ajuste afín consume 2 grados de libertad).

**Condición 3: modelo de reloj.** $\varepsilon$ requiere ≥ 2 ventanas
separadas y $\gamma$ ≥ 3. Precisión aproximada:
$\sigma_\varepsilon \approx \sigma_t\sqrt{12/J}/T_{obs}$, que depende sobre
todo del *intervalo* $T_{obs}$ y poco de cuántas ventanas haya.

**Condición 4: fase por ventana.** $2\pi f_c\delta(t_j)+\phi_0$ es observable
**módulo $2\pi$** en cada ventana que contenga el directo, y módulo $2\pi$ es
justo lo que necesita el BP. Si el directo y los ecos están en la misma
ventana no hace falta desenrollar fase entre ventanas. El ruido de reloj entre
la llegada del directo y la de los ecos (< 1 µs) es despreciable:
$\sigma_y\sqrt{\Delta t}\,2\pi f_c\sim 10^{-3}$ rad incluso con XO.

**Condición 5: aproximaciones del modelo.** $\varepsilon f_c T_p \ll 1$ (Doppler
intra-pulso del offset; 0.01 ciclos con 20 ppm) y $\varepsilon\,T_{ventana}\ll 1/B$
(estiramiento del muestreo dentro de la ventana; 0.2 ps).

**Requisitos de hardware que se siguen de lo anterior:**

* línea de vista Tx→Rx2;
* el directo y los ecos caben en la misma ventana y en el rango dinámico del ADC (el directo es mucho más fuerte);
* mismo oscilador para el LO y el muestreo del Rx2;
* el Tx registra $t_n$ con stagger.

---

## 3. Precisión temporal necesaria para un BiSAR enfocado

Sea $e(t) = \hat\delta(t) - \delta(t)$ el error residual de la base de tiempos.
Actúa a la vez sobre la envolvente (el BP muestrea el dato comprimido en el
instante equivocado) y sobre la fase ($2\pi f_c e$). Como $f_c/B = 8.5$, **la
fase es siempre el requisito dominante**.

| Componente de $e(t)$ | Efecto | Ley | Tolerancia (este caso)¹ |
|---|---|---|---|
| Constante $\Delta\tau$ | desplazamiento a lo largo de $\nabla R_b$ | $\lvert\Delta\mathbf p\rvert \approx c\,\Delta\tau/\lvert\nabla_g R_b\rvert$ | **≈ 1.1 ns** (0.25 m) |
| Constante de fase | ninguno en una imagen; sesgo de fase absoluta en la fusión mono+bi | — | libre / < 0.3 rad para fusión |
| Lineal $\Delta\varepsilon\,(t-t_c)$ | fase lineal en tiempo lento → desplazamiento en azimut | $\lvert\Delta x\rvert \approx c\,\Delta\varepsilon\,R_{Tx}/v$ | **≈ 5·10⁻¹¹ (0.02 Hz en $f_c$)** |
| Cuadrática (QPE en el borde) | desenfoque, ensanchamiento | pérdida $=\lvert\int w(u)e^{jQu^2}du\rvert/\int w$ | **QPE ≈ π (Hann)**, π/4 sin ponderación |
| Fase blanca $\sigma_\phi$ | pérdida de pico, suelo de lóbulos | $L = e^{-\sigma_\phi^2/2}$ (−4.34 σ² dB) | **σφ ≈ 0.47 rad** (1 dB) |
| Jitter de tiempo $\sigma_t$ | ídem vía $\sigma_\phi = 2\pi f_c\sigma_t$ | ídem | **σt ≈ 175 ps** (1 dB) |
| Índice de pulso erróneo | geometría y tiempo de emisión erróneos | — | catastrófico (debe ser 0) |

¹ Valores de exp1 (Fig. `exp1_tolerancias`), ajustados a la teoría
(desviación < 5 %). Con $\nabla_g R_b = 1.54$, la resolución en rango sobre
el suelo es $c/(B\lvert\nabla_g R_b\rvert) = 3.9$ m (5.8 m con Hann) y la
resolución en azimut ≈ 2.2 m (Hann).

Conclusiones del cuadro:

* La envolvente es fácil: un error de ~1 ns ya es una fracción de la
  resolución, y el directo da ~0.3 ns a 30 dB.
* La fase es lo difícil: hacen falta **~175 ps rms de coherencia pulso a
  pulso** (0.47 rad a 425 MHz, unos 5 cm de camino eléctrico) y una
  frecuencia residual de **~5·10⁻¹¹**. Un TCXO libre (1 ppm) está 2·10⁴
  veces por encima, y su desplazamiento azimutal equivalente es de km.
* Estas tolerancias escalan con $1/f_c$: en banda X serían 23 veces más
  estrictas.

---

## 4. Cómo explotar el camino directo

Algoritmo implementado (`lib/simula_adquisicion.m`, `lib/sincroniza.m`):

1. **Adquisición.** La primera ventana graba una PRT completa. Se aplica el
   filtro adaptado y se detecta el *primer* máximo que supera el umbral
   (13 dB sobre el ruido) y es máximo local en ±20 muestras. Así se toma el
   directo, que llega primero, y no un eco ni un lóbulo de rango.
2. **Puertas de seguimiento.** El Rx2 solo necesita grabar alrededor de
   $\hat a_1 + m\,\mathrm{PRT}$, con un margen
   $2S + \varepsilon_{max}T + \tfrac12\gamma_{max}T^2 + 5\sigma_y\sqrt T$
   (±20 µs para TCXO). Este margen solo usa información del receptor.
3. **Estimación ML del directo en cada ventana.** Newton sobre $\lvert r(\tau)\rvert^2$,
   con $r(\tau)$ el filtro adaptado interpolado exactamente por Fourier, da
   $\hat a_j$ y $\hat r_j = r(\hat a_j)$. Cotas de Cramér-Rao:
   $$\sigma_t \ge \frac{1}{2\pi\beta\sqrt{2\,\mathrm{SNR}}},\qquad
     \sigma_\phi \ge \frac{1}{\sqrt{2\,\mathrm{SNR}}},$$
   con $\beta = 12.7$ MHz (ancho de banda rms del chirp con flancos). A 30 dB:
   $\sigma_t = 300$ ps y $\sigma_\phi = 0.024$ rad.
4. **Asociación.** Para cada hipótesis $k$ se ajusta
   $\hat a_j = \alpha + \beta\,[t_{k+m_j} + \tau_d(k+m_j)]$ por LS y se elige
   el residuo mínimo (`lib/asocia_pulsos.m`).
5. **Base de tiempos por ventana.**
   $\hat\delta_j = \hat a_j - t_{n_j} - \tau_d(n_j)$.
6. **Referencia de fase por ventana (E3).**
   $\hat\psi_j = -[\arg\hat r_j + 2\pi f_c\tau_d(n_j)] = 2\pi f_c\delta_j+\phi_0+\phi_{LO,j}$.
   Es el principio del canal de referencia en radar pasivo. Cancela también
   el ruido de fase del Tx y del LO del Rx, porque son comunes al directo y a
   los ecos de la misma ventana.
7. **Cancelación del directo** (se resta $\hat r_j\,s(t-\hat a_j)/E_s$),
   compresión con Hann y BP:
   $$I(\mathbf g) = \sum_j w_j\,\mathrm{RC}_j\big(t_{n_j}+\tau_b+\hat\delta_j\big)\,e^{\,j(2\pi f_c\tau_b + \hat\psi_j)}.$$

Comparación de estrategias (todas sin verdad terreno):

| | Retardo usado | Fase usada | Error de fase esperado |
|---|---|---|---|
| E0 | 1.ª ventana, constante | constante | $2\pi f_c\,\varepsilon\,t$ → falla |
| E1 | modelo $\delta_0+\varepsilon t+\gamma t^2/2$ | $2\pi f_c\hat\delta(t)$ | ruido del oscilador no modelado $\propto 2\pi f_c\sigma_y\sqrt{T}$ |
| E2 | medido en cada ventana | $2\pi f_c\hat\delta_j$ | $2\pi f_c\sigma_t = (f_c/\beta)/\sqrt{2\mathrm{SNR}}$ |
| E3 | medido en cada ventana | **fase portadora medida** | $1/\sqrt{2\mathrm{SNR}}$ |

E3 mejora a E2 en un factor $(f_c/\beta)^2 = 1124$ (**30.5 dB de SNR**). Para
σφ ≤ 0.47 rad, E2 necesita SNR_d ≳ 35 dB y E3 solo ≳ 4 dB. En la práctica,
para E3 la limitación es la detección (≳ 15 dB), no la precisión.

---

## 5. Pulsos perdidos (escucha discontinua)

Se separan tres efectos distintos:

1. **Sobre el sincronismo de E3: ninguno**, siempre que cada ventana grabada
   contenga el directo. Cada ventana se autorreferencia. No hay que desenrollar
   fase ni extrapolar el reloj a través de los huecos.
2. **Sobre la asociación y las puertas.** Hacen falta $J\ge3$ ventanas (en la
   práctica ≥ 4) y stagger. Durante un hueco $T_g$ el reloj se desplaza
   $\varepsilon_{max}T_g$, y la puerta debe cubrirlo o se re-adquiere grabando
   una PRT completa (XO de 20 ppm: 20 µs por segundo de hueco).
3. **Sobre los métodos basados en modelo (E1, o la predicción de la fase de
   una ventana sin directo).** El error de predicción en el centro de un hueco
   es un puente browniano:
   $$\sigma_\phi(T_g) \approx \sqrt{(2\pi f_c\sigma_y)^2\,T_g/6 + \sigma_w^2},\qquad
     T_{g,max} \approx 6\left(\frac{0.47}{2\pi f_c\sigma_y}\right)^2 .$$
   Con eso: OCXO ~ 10³ s, TCXO ~ 19 s, XO ~ 0.2 s y XO barato ~ 2 ms (menos
   que una PRT). Por encima aparecen saltos de ciclo.
4. **Sobre la imagen (no es sincronismo).** Con menos pulsos se pierden
   $10\log_{10}(J_{full}/J)$ dB de SNR integrada. El muestreo aleatorio sube
   el suelo de lóbulos (~$-10\log_{10}J$). El muestreo periódico 1/G produce
   réplicas azimutales a $\lambda R_{Tx}/(v\,\mathrm{PRT}\,G)$ (28 m con
   G = 128). Este problema no se arregla sincronizando mejor; por eso la
   simulación lo mide aparte (ideal con el mismo patrón frente a ideal
   continuo).

---

## 6. Resultados de la simulación

Todas las cifras proceden de `run_all.m`. Los CSV completos están en
`resultados/` y las figuras en `resultados/figs/`. Las métricas de imagen
se calculan sobre parches de ±12 m alrededor de cada blanco y se comparan
con la imagen formada con el **reloj verdadero y el mismo conjunto de
pulsos**. Así la pérdida medida se debe solo al sincronismo.

### 6.1 Precisión necesaria (exp1, `exp1_tolerancias`)

| Error inyectado | Umbral simulado | Teoría |
|---|---|---|
| Offset constante → 0.25 m de desplazamiento | 1.17 ns | $c\Delta\tau/\lvert\nabla R_b\rvert$: coincide (< 3 %) |
| Frecuencia residual → 0.25 m | 5.3·10⁻¹¹ (0.022 Hz) | sistema $[\nabla R_b;\nabla\dot R_b]$: coincide |
| QPE → 1 dB | 3.2 rad (≈ π, con Hann) | integral de Fresnel ponderada: coincide |
| Fase blanca → 1 dB | 0.47 rad | $e^{-\sigma^2/2}$: coincide |
| Jitter blanco → 1 dB | 177 ps | $\sigma_\phi = 2\pi f_c\sigma_t$: coincide |

La imagen ideal tiene PSLR −32.7 dB e ISLR −31.9 dB. La simulación
reproduce las cinco leyes de §3 en todo el rango, así que las tolerancias
pueden usarse como especificación.

### 6.2 Cadena completa (exp2, `exp2_reloj_asociacion`, `exp2_imagenes`)

TCXO (ε = −0.60 ppm), δ0 = 15.8 s, escucha aleatoria del 51 % (1222 de
2400 PRT), SNR_d = 30 dB.

* Asociación correcta (k = 161), con un margen de **10⁴** entre la mejor
  hipótesis y la segunda.
* ε estimado con un error de 1.4·10⁻¹⁰.

| Estrategia | σt (ps) | σφ (rad) | desplaz. (m) | pérdida (dB) | PSLR | ISLR | ρ | ¿Enfoca? |
|---|---|---|---|---|---|---|---|---|
| E0 offset | 2.1·10⁶ | 1.79 | 9.2 | −62.6 | −1.0 | 4.0 | 0.05 | **NO** |
| E1 modelo | 50 | 0.13 | 0.12 | −0.08 | −23.9 | −19.5 | 0.994 | SÍ |
| E2 tiempo | 305 | 0.81 | 0.10 | −2.95 | −21.8 | −16.0 | 0.991 | **NO** |
| E3 fase | 305 | **0.023** | 0.004 | **−0.00** | −25.7 | −20.2 | 1.000 | **SÍ** |
| Ideal | 0 | 0 | 0 | 0 | −25.8 | −20.3 | 1 | — |

La fusión coherente con la imagen monostática da una diferencia de fase
entre picos de 0–1° con E3, igual que con el reloj ideal. Es decir, el
referenciado por el directo ancla también la fase absoluta. Esto es
relevante para la fusión mono+bi de la tesis.

### 6.3 Precisión alcanzada frente a SNR y oscilador (exp3, `exp3_snr`, `exp3_oscilador`)

Medianas de 4 realizaciones; TCXO; escucha continua.

| SNR_d (dB) | σt medido / CRB (ps) | σφ E2 (rad) | σφ E3 (rad) | Pérdida E2 (dB) | Pérdida E3 (dB) | Cumple E1 / E2 / E3 |
|---|---|---|---|---|---|---|
| 10 | 4·10⁵ / 3001 | 1.80 | 1.74 | −40.4 | −33.8 | 0 / 0 / 0 |
| 15 | 1.5·10⁵ / 1688 | 1.80 | 0.66 | −29.3 | −1.2 | 0 / 0 / 0 |
| 20 | 952 / 949 | 1.76 | 0.076 | −25.2 | −0.04 | 0 / 0 / **1** |
| 30 | 300 / 300 | 0.79 | 0.024 | −2.75 | −0.01 | 0 / 0 / **1** |
| 40 | 95 / 95 | 0.25 | 0.008 | −0.28 | −0.00 | 0 / 0.5 / **1** |
| 50 | 30 / 30 | 0.08 | 0.002 | −0.03 | −0.00 | 0 / 1 / **1** |

Lectura de la tabla:

* El estimador ML del directo **alcanza la CRB** desde 20 dB, tanto en
  retardo como en fase.
* **Por debajo de ~15–20 dB falla la detección, no la estimación.** Los ecos
  (20 dB por pulso en esta escena) superan el umbral y el detector toma un
  eco por el directo. El error es entonces de cientos de ns, y la columna
  `det_ok` del CSV sigue a 1 porque algo cruzó el umbral. En el diseño real,
  el directo tiene que ser claramente más fuerte que cualquier eco, lo cual
  es lo habitual con línea de vista.
* E2 sigue exactamente $2\pi f_c\,\mathrm{CRB}_t$ y necesita ≥ 40–50 dB.
  E3 sigue $1/\sqrt{2\mathrm{SNR}}$ y cumple desde 20 dB. La diferencia
  medida es la prevista de ~30 dB.

Por clase de oscilador (SNR_d = 30 dB):

| Oscilador | σφ E1 (rad) | Pérdida E1 (dB) | ΔPSLR E1 (dB) | σφ E3 (rad) | Cumple E1 / E2 / E3 |
|---|---|---|---|---|---|
| OCXO | 0.023 | −0.00 | +0.2 | 0.024 | **1** / 0 / **1** |
| TCXO | 0.14 | −0.09 | **+8** | 0.024 | 0 / 0 / **1** |
| XO | 1.32 | −5.7 | +30 | 0.024 | 0 / 0 / **1** |
| XO barato | 1.78 | −23 | +31 | 0.024 | 0 / 0 / **1** |

Con TCXO, E1 casi no pierde pico. Aun así falla el criterio porque el error
de fase del paseo aleatorio es de **baja frecuencia** y sube los lóbulos
cercanos 7–10 dB. Un modelo paramétrico de reloj solo basta con OCXO. **E3
es independiente del oscilador.**

### 6.4 Pulsos perdidos (exp4, `exp4_pulsos_perdidos`)

Medianas de 3 realizaciones; TCXO; SNR_d = 30 dB.

| Patrón de escucha | J | Asociación | Pérdida E3 vs ideal (dB) | E3 cumple | ΔPSLR por muestreo (dB) | ρ ideal-patrón vs continuo |
|---|---|---|---|---|---|---|
| continuo | 2400 | ok | −0.005 | sí | 0 | 1.000 |
| 1 de 8 | 300 | ok | −0.012 | sí | +2.0 | 0.999 |
| 1 de 32 | 75 | ok | −0.007 | sí | +5.8 | 0.993 |
| 1 de 128 | 19 | ok | −0.017 | sí | +22.6 | 0.71 |
| aleatorio 50 % | 1184 | ok | −0.005 | sí | +5.6 | 0.996 |
| aleatorio 5 % | 119 | ok | −0.006 | sí | +19.8 | 0.87 |
| aleatorio 1 % | 27 | ok | −0.037 | sí | +25.0 | 0.60 |
| ráfagas 20/180 | 240 | ok | −0.006 | sí | +29.4 | 0.60 |
| ráfagas 5/495 | 25 | ok | −0.017 | sí | +32.4 | 0.32 |

* **El sincronismo E3 no se degrada con ningún patrón**, incluso con 19–25
  ventanas en 6 s. Todo el deterioro de la imagen se debe al muestreo de la
  apertura (columnas de muestreo): ruido, suelo de lóbulos con muestreo
  aleatorio, y réplicas y lóbulos de subapertura con patrones periódicos o
  ráfagas. E2 queda en torno a −3 dB en todos los patrones (su límite es la
  SNR, no los huecos).
* Puenteo de fase (ventanas sin directo utilizable) frente a la teoría del
  puente browniano:

| Oscilador | Hueco máximo simulado (σφ ≤ 0.47 rad) | Teoría $6(0.47/2\pi f_c\sigma_y)^2$ |
|---|---|---|
| OCXO | > 1.3 s (límite del ensayo) | 1.9·10³ s |
| TCXO | > 1.3 s (límite del ensayo; σφ = 0.12 rad a 1.28 s) | 18 s |
| XO | 0.08 s (32 PRT) | 0.19 s |
| XO barato | < 2 PRT (saltos de ciclo) | 2 ms |

  La curva simulada sigue a la teórica. El hueco máximo simulado es ~2 veces
  menor que el teórico porque la teoría usa el error medio del puente y no el
  del centro del hueco.

### 6.5 Asociación (exp5, `exp5_asociacion`)

P(asociación correcta) con 200 ensayos por casilla, σt = 300 ps y ventanas
aleatorias en 6 s.

| Stagger S | J = 3 | J = 4 | J = 5 | J = 8 | J = 16 |
|---|---|---|---|---|---|
| 0 (PRF uniforme) | 0.00 | 0.00 | 0.00 | 0.01 | 0.00 |
| 10 ns | 0.02 | 0.45 | 0.94 | 1.00 | 1.00 |
| 1 µs | 0.80 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5 µs | 0.93 | 1.00 | 1.00 | 1.00 | 1.00 |

(TCXO; con XO barato hacen falta S ≥ 1 µs y J ≥ 5, porque la deriva
cuadrática no modelada en el ajuste afín actúa como ruido.)

**Sin stagger, la asociación es imposible** aunque haya cientos de ventanas,
como predice §2. Con S = 5 µs (0.2 % de la PRT) bastan 4 ventanas.

### 6.6 Error de navegación del Tx (exp6, `exp6_navegacion`)

Pérdida de pico frente a la misma adquisición sin error de navegación
(AR(1), τc = 1 s, σ por eje).

| σnav | Bi E2 | **Bi E3** | Bi reloj ideal | Mono |
|---|---|---|---|---|
| 1 cm | −2.73 | −0.02 | −0.03 | −0.09 |
| 2 cm | −2.78 | −0.05 | −0.10 | −0.36 |
| 5 cm | −3.03 | **−0.27** | −0.54 | −2.05 |
| 10 cm | −3.61 | **−1.05** | −1.95 | −6.24 |

El referenciado por el directo **absorbe la parte del error del Tx
proyectada sobre la línea Tx→Rx**. Por eso E3 tolera ~1.5 veces más error de
navegación que el bistático con reloj perfecto y ~3 veces más que el
monostático. El error de navegación no es un obstáculo adicional para el
sincronismo: el cuello de botella sigue siendo el monostático
(σnav ≲ 2–3 cm a 425 MHz).

---

### 6.7 Conclusión de viabilidad

1. **Es viable.** Con el camino directo como única referencia, la estrategia
   E3 (retardo + fase portadora del directo en cada ventana) produce imágenes
   BiSAR indistinguibles de las del reloj ideal: pérdida < 0.05 dB,
   desplazamiento < 1 cm, ρ > 0.999. Lo hace para cualquier oscilador
   (OCXO a XO barato), arranque y offset arbitrarios y escucha discontinua
   (hasta el 1 % de las PRT).
2. **Condiciones necesarias:**
   * SNR del directo ≳ 20 dB por pulso y claramente por encima de los ecos;
   * stagger (o código) en la emisión, registrado en el log del Tx, de
     S ≳ 1 µs;
   * al menos 4 ventanas;
   * directo y ecos en la misma ventana, con el mismo oscilador para el LO y
     el muestreo del Rx2;
   * distancia Tx–Rx conocida con error ≪ λ (lo cubre la propia navegación
     que exige el monostático).
3. **El sincronismo solo de tiempo no basta.** El retardo del directo da
   ~300 ps a 30 dB, por encima de la tolerancia de ~177 ps a 425 MHz. Hace
   falta la fase portadora, que es (f_c/β)² ≈ 30 dB más precisa.
4. **Un modelo paramétrico de reloj (offset + frecuencia + deriva) solo
   basta con OCXO.** Con TCXO o peor, el ruido del oscilador exige corregir
   ventana a ventana.
5. **Los pulsos perdidos no afectan al sincronismo E3.** Sí afectan a la
   imagen por muestreo de la apertura, y eso no es un problema de
   sincronismo. Solo si una ventana carece de directo hay que puentear su
   fase, y entonces el hueco máximo lo fija el oscilador (TCXO ≳ 1 s,
   XO ≈ 0.1 s, XO barato < 1 PRT).

---

## 7. Métricas para demostrar experimentalmente el sincronismo

En un experimento real no hay "reloj ideal". Las métricas tienen que ser
auto-consistentes o apoyarse en reflectores de referencia.

**Nivel sincronismo (sin verdad terreno):**

* Margen de asociación: residuo de la 2.ª mejor hipótesis frente a la mejor (≫ 10).
* Residuo del retardo del directo tras el modelo E1 frente a la CRB a la SNR
  medida. Si es mayor que la CRB, hay ruido de oscilador o multicamino.
* **Validación cruzada:** se ajusta el reloj con las ventanas impares y se
  predicen retardo y fase de las pares (error de predicción frente a las
  tolerancias de §3).
* Desviación de Allan de $\hat\delta(t)$ frente a la hoja de datos del oscilador.
* Consistencia de $\hat\varepsilon$ entre subaperturas o pasadas.
* Dispersión de la fase del directo tras quitar la geometría (E3 da
  directamente su $\sigma_\phi$ residual).

**Nivel imagen** (reflectores triédricos con posición topográfica):

* Error de posición de los picos (< 0.1 de la resolución; aquí 0.25 m).
* Resolución −3 dB, PSLR e ISLR frente a la IRF teórica de la geometría
  (no frente a un ideal simulado).
* Pérdida de pico frente al presupuesto de enlace del reflector (calibración
  radiométrica), o frente a la imagen formada con el modelo E1 y con E3:
  si coinciden, el ruido de fase no domina.
* Coherencia entre subaperturas (dos mitades de la apertura): una
  $\lvert\rho\rvert$ alta indica que no hay error de fase de baja frecuencia.
* Contraste y entropía de la imagen: una mejora tras autofocus indica
  sincronismo insuficiente. Con un sincronismo correcto, el autofocus no
  debería mejorar la imagen.
* **Coherencia mono–bi:** diferencia de fase de los picos entre la imagen
  monostática y la bistática (E3 ancla la fase absoluta), y ganancia de la
  fusión coherente frente a la incoherente.
* Repetibilidad entre pasadas (registro de los reflectores).

**Criterios de éxito usados en la simulación** (`P.crit`, frente al reloj
ideal en el peor blanco): pérdida de pico ≤ 1 dB, desplazamiento ≤ 0.25 m,
ΔPSLR ≤ 3 dB, ΔISLR ≤ 3 dB, coherencia ρ ≥ 0.95.

---

## 8. Limitaciones del modelo simulado

* Blancos puntuales isótropos, sin diagramas de antena ni clutter
  distribuido.
* Sin multicamino en el camino directo. El rebote en el suelo del directo
  (dos rayos) sesgaría retardo y fase; en la práctica hay que separarlo en
  rango (la diferencia de camino debe ser ≫ c/B) o elegir alturas que lo
  eviten.
* Stop-and-hop; el offset de frecuencia intra-pulso no se modela (justificado
  en §2, condición 5).
* No se modelan la saturación ni la cuantización del ADC del Rx2 (directo
  fuerte + ecos débiles).
* El Rx2 está fijo y su posición se conoce. Si estuviera en movimiento, su
  error de navegación entraría igual que el del Tx.
* El oscilador se modela con FM blanca + deriva. El ruido flicker haría el
  modelo E1 algo peor; no afecta a E3.
