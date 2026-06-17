# Fundamento teórico: Reflexión y Transmisión TM en un Sistema Multicapa

Este documento describe el fundamento teórico, las ecuaciones y las premisas detrás de
`multislab_TM.m`, que calcula los coeficientes de reflexión (Γ) y transmisión (τ), así
como la reflectancia (R) y transmitancia (T), para una onda plana con polarización TM
(paralela) incidiendo de forma oblicua sobre un sistema de **N slabs dieléctricos**
apilados entre dos semiespacios de vacío.

---

## 1. Configuración del sistema y premisas

```
   Medio 0        Medio 1     Medio 2   ...   Medio N        Medio N+1
   (vacío)        (Slab 1)    (Slab 2)        (Slab N)       (vacío)
  n0, η0  θ1  →  n1,η1  θ2 → n2,η2  θ3 →  ... → nN,ηN  θN+1 → n0,η0  θ1
        z=0          z=-d1      z=-(d1+d2)   ...   z=-Σdi
```

**Premisas del modelo:**

1. Onda plana monocromática de frecuencia `f`, incidiendo desde un medio semi-infinito
   de **vacío** (medio 0: `n0 = 1`, `η0 = sqrt(μ0/ε0)`).
2. El sistema termina en otro semiespacio de **vacío** (medio N+1, idénticas
   propiedades al medio 0). Esto es clave para varias simplificaciones (Sección 6).
3. Cada slab `i = 1..N` es un medio dieléctrico, no magnético (`μr = 1`), sin pérdidas,
   caracterizado por:
   - Permitividad relativa `eps_r(i)`
   - Espesor `d(i)`
   - Índice de refracción `n(i) = sqrt(eps_r(i))`
   - Impedancia intrínseca `eta(i) = eta0 / n(i)`
4. Las interfaces son planas, infinitas y paralelas entre sí (perpendiculares al eje `z`).
5. Polarización **TM** (campo eléctrico contenido en el plano de incidencia, también
   llamada polarización "paralela" o "p").
6. Todo el análisis se realiza en régimen permanente (fasores, dependencia temporal
   `e^{jωt}` implícita).

---

## 2. Geometría de propagación: Ley de Snell generalizada

Para una onda plana oblicua, la componente del vector de onda **paralela** a las
interfaces (`k_x`) se conserva al atravesar cada interfaz. Esto es la **ley de Snell**:

```
n_prev * sin(theta_prev) = n_curr * sin(theta_curr)
```

Aplicada sucesivamente desde el medio de incidencia (vacío, `n0`, ángulo `θ1`) hasta el
interior del slab `i`, se obtiene el ángulo de refracción `θ_{i+1}` dentro de cada slab:

```
theta(1) = theta1                        (ángulo de incidencia, en vacío)
theta(2) = asin( n0    / n(1) * sin(theta(1)) )   (ángulo dentro del Slab 1)
theta(3) = asin( n(1) / n(2) * sin(theta(2)) )   (ángulo dentro del Slab 2)
   ...
theta(N+1) = asin( n(N-1) / n(N) * sin(theta(N)) ) (ángulo dentro del Slab N)
```

En el script esto corresponde al bucle:

```matlab
for i = 2:(N_slabs + 1)
    if i == 2
        n_prev = n0;
    else
        n_prev = n(i-2);
    end
    n_curr = n(i-1);
    arg = n_prev / n_curr * sin(thetas(i-1));
    thetas(i) = asin(arg);
end
```

### Invariante global (conservación de `n·sinθ`)

Como `n·sin(θ)` se conserva en **cada** interfaz, también se conserva de extremo a
extremo del sistema:

```
n0 * sin(theta1) = n(i) * sin(theta_{i+1})   para todo i = 1..N
```

Dado que el medio de salida (N+1) también es vacío (`n0`), se cumple:

```
n0 * sin(theta_salida) = n0 * sin(theta1)   =>   theta_salida = theta1
```

**Esta es la base de la simplificación `Z_out = Z_in` usada más adelante** (Sección 6).

### Reflexión interna total (TIR)

Si en algún paso `n_prev/n_curr > 1` (pasando de un medio más denso óptimamente a uno
menos denso) y el ángulo de incidencia supera el ángulo crítico, entonces:

```
arg = (n_prev/n_curr) * sin(theta_prev) > 1
```

y `asin(arg)` no tiene solución real → toda la energía se refleja:

```
Gamma_TM = 1,   Tau_TM = 0   (R = 1, T = 0)
```

El script detecta esta condición (`abs(arg) > 1`) y corta el cálculo para ese ángulo
(`TIR_flag`). Nótese que la transición del slab N hacia el vacío de salida **nunca**
produce TIR, porque por el invariante global `n0·sin(θ1) ≤ n0` siempre (ya que
`sin(θ1) ≤ 1`), de modo que siempre existe un `θ_salida` real.

---

## 3. Impedancia de onda y número de onda transversal (TM)

Para una onda plana viajando con ángulo `θ` respecto a la normal (eje `z`) en un medio
con impedancia intrínseca `η = sqrt(μ/ε)` e índice `n = sqrt(eps_r)`, se definen:

**Número de onda en la dirección de propagación normal (z):**

```
k0 = 2*pi*f / c                     (número de onda en el vacío)
k_z(i) = k0 * n(i) * cos(theta_{i+1})
```

**Impedancia de onda transversal para polarización TM (paralela):**

```
Z_TM(i) = eta(i) * cos(theta_{i+1})
```

> Convención utilizada (análoga a Pozar, *Microwave Engineering*): para TM,
> `Z_TM = η·cosθ`; para TE sería `Z_TE = η/cosθ`. Esta convención es consistente con el
> modelo de "línea de transmisión equivalente", donde cada slab se representa como un
> tramo de línea con impedancia característica `Z_TM(i)` y longitud eléctrica
> `k_z(i)·d(i)`.

En el script:

```matlab
k_z = (2*pi*f * n(i) / c) * cos(thetas(i+1));
Z_TM = eta(i) * cos(thetas(i+1));
```

---

## 4. Modelo de línea de transmisión equivalente (matriz ABCD)

Cada slab `i` se modela como un tramo de línea de transmisión sin pérdidas de:

- Impedancia característica: `Z_TM(i)`
- Longitud eléctrica: `β·l = k_z(i) · d(i)`

La matriz de transmisión (parámetros ABCD) que relaciona los fasores de
tensión/corriente equivalentes (campos transversales `E_x`, `H_y`) en la entrada y la
salida del tramo es la matriz clásica de una línea sin pérdidas:

```
        ┌                                      ┐
M_i  =  │  cos(k_z d_i)        j Z_TM sin(k_z d_i)  │
        │  j sin(k_z d_i)/Z_TM    cos(k_z d_i)       │
        └                                      ┘
```

Esto corresponde exactamente a:

```matlab
M_i = [cos(k_z * d(i)),       1j*Z_TM*sin(k_z * d(i));
       1j/Z_TM*sin(k_z * d(i)), cos(k_z * d(i))];
```

**Propiedad importante:** para un medio sin pérdidas, `det(M_i) = cos^2 + sin^2 = 1`.
Esto se usa más adelante para verificar la conservación de energía.

---

## 5. Matriz ABCD total del sistema (cascada de N capas)

La gran ventaja del formalismo ABCD es que la matriz total de una cascada de `N`
tramos es simplemente el **producto matricial** de las matrices individuales, en el
mismo orden físico en que la onda las atraviesa:

```
M_tot = M_1 · M_2 · ... · M_N  =  [ A  B ]
                                   [ C  D ]
```

En el script, esto se construye acumulando el producto desde la matriz identidad:

```matlab
M_tot = eye(2,2);
for i = 1:N_slabs
    ... % construir M_i
    M_tot = M_tot * M_i;
end
A = M_tot(1,1);  B = M_tot(1,2);
C = M_tot(2,1);  D = M_tot(2,2);
```

Como cada `M_i` tiene determinante 1 (sin pérdidas), también:

```
det(M_tot) = A*D - B*C = 1
```

---

## 6. Impedancias de los medios de entrada y salida

El sistema completo se ve, desde el punto de vista de la línea de transmisión
equivalente, como una red de dos puertos `M_tot` cargada en su salida con la
impedancia del medio de salida `Z_out`, y excitada desde el medio de entrada con
impedancia `Z_in`.

```
Z_in  = eta0 * cos(theta1)        (impedancia TM del vacío de entrada, ángulo theta1)
Z_out = eta0 * cos(theta_salida)  (impedancia TM del vacío de salida, ángulo theta_salida)
```

Por el invariante de Snell de la Sección 2, `theta_salida = theta1` (ambos extremos son
vacío con el mismo `n0`), por lo tanto:

```
Z_out = eta0 * cos(theta1) = Z_in
```

Esta igualdad **Z_in = Z_out** es válida siempre en este sistema (medio de entrada y
salida idénticos) y simplifica directamente las fórmulas de reflexión/transmisión y el
cálculo de la transmitancia (Sección 8).

```matlab
Z_in  = eta0 * cos(th1);
Z_out = eta0 * cos(th1);  % = Z_in, por simetría vacío-vacío
```

---

## 7. Coeficiente de reflexión Γ_TM a partir de ABCD

Para una red de dos puertos con parámetros `[A B; C D]`, cargada en el puerto de salida
con impedancia `Z_out` y vista desde el puerto de entrada con impedancia de referencia
`Z_in`, la impedancia de entrada de la red completa es:

```
Z_entrada_red = (A * Z_out + B) / (C * Z_out + D)
```

El coeficiente de reflexión a la entrada se define, como en cualquier línea de
transmisión, por la discontinuidad de impedancia respecto a la línea de alimentación
(`Z_in`):

```
Gamma_TM = (Z_entrada_red - Z_in) / (Z_entrada_red + Z_in)
```

Sustituyendo `Z_entrada_red` y multiplicando numerador y denominador por
`(C·Z_out + D)`:

```
                A*Z_out + B - Z_in*(C*Z_out + D)
Gamma_TM = ---------------------------------------
                A*Z_out + B + Z_in*(C*Z_out + D)
```

Como `Z_out = Z_in` (Sección 6), esta expresión se reduce a la forma implementada en el
script (sustituyendo `Z_out → Z_in`):

```
                A*Z_in + B - C*Z_in^2 - D*Z_in
Gamma_TM = -----------------------------------
                A*Z_in + B + C*Z_in^2 + D*Z_in
```

```matlab
num = A*Z_in + B - C*Z_in^2 - D*Z_in;
den = A*Z_in + B + C*Z_in^2 + D*Z_in;
Gamma_TM(idx) = num / den;
```

---

## 8. Coeficiente de transmisión τ_TM y conservación de energía

El coeficiente de transmisión de tensión/campo equivalente para la misma red de dos
puertos, terminada en `Z_out` y alimentada desde `Z_in`, es:

```
                       2 * Z_out
Tau_TM = -----------------------------------------
           A*Z_in + B + C*Z_in*Z_out + D*Z_out
```

```matlab
Tau_TM(idx) = 2*Z_out / (A*Z_in + B + C*Z_in*Z_out + D*Z_out);
```

### Reflectancia y transmitancia (potencia)

```
R_TM = |Gamma_TM|^2
T_TM = |Tau_TM|^2 * (Re{Z_in} / Re{Z_out})
```

Como el sistema **no tiene pérdidas** y `Z_in = Z_out` (medios de entrada y salida
idénticos, Sección 6), el factor de impedancias es 1, y la transmitancia se reduce a:

```
T_TM = |Tau_TM|^2
```

```matlab
R_TM = abs(Gamma_TM).^2;
T_TM = abs(Tau_TM).^2;     % válido porque Z_in = Z_out
```

### Conservación de energía

Para un sistema sin pérdidas (`det(M_tot) = 1`, todos los `Z_i` reales), se cumple
idénticamente:

```
R_TM + T_TM = 1     (para todo theta1 sin TIR)
```

Esta identidad es la verificación numérica fundamental del modelo: si
`R_TM + T_TM ≠ 1`, hay un error en la definición de impedancias, ángulos o en la matriz
ABCD. (Verificado numéricamente con el sistema de 3 slabs de ejemplo:
`max|R_TM + T_TM - 1| ≈ 6.7e-16`, es decir, error de redondeo de máquina.)

En el caso de **TIR** (Sección 2), se impone directamente `Gamma_TM = 1`, `Tau_TM = 0`,
por lo que `R_TM + T_TM = 1` también se cumple trivialmente.

---

## 9. Resumen del algoritmo (correspondencia con el código)

Para cada ángulo de incidencia `theta1`:

1. **Snell**: calcular `theta(i)` para `i = 1..N+1` (ángulo en vacío de entrada y en
   cada slab). Si `|n_prev/n_curr * sin(theta_prev)| > 1` en algún paso → TIR
   (`Gamma=1`, `Tau=0`, fin).
2. **Por cada slab `i = 1..N`**:
   - `k_z(i) = (2*pi*f*n(i)/c) * cos(theta_{i+1})`
   - `Z_TM(i) = eta(i) * cos(theta_{i+1})`
   - Construir `M_i` (matriz ABCD del tramo de línea equivalente al slab `i`)
3. **Cascada**: `M_tot = M_1 * M_2 * ... * M_N = [A B; C D]`
4. **Impedancias de los extremos**: `Z_in = Z_out = eta0 * cos(theta1)`
5. **Coeficientes**:
   - `Gamma_TM = (A*Z_in + B - C*Z_in^2 - D*Z_in) / (A*Z_in + B + C*Z_in^2 + D*Z_in)`
   - `Tau_TM   = 2*Z_out / (A*Z_in + B + C*Z_in*Z_out + D*Z_out)`
6. **Potencias**: `R_TM = |Gamma_TM|^2`, `T_TM = |Tau_TM|^2`
7. **Verificación**: `R_TM + T_TM = 1`

---

## 10. Variables principales del script y su significado físico

| Variable          | Significado                                                        |
|-------------------|---------------------------------------------------------------------|
| `f`               | Frecuencia de la onda incidente [Hz]                                |
| `c`               | Velocidad de la luz en el vacío [m/s]                               |
| `eps0`, `mu0`     | Permitividad y permeabilidad del vacío                              |
| `eta0`, `n0`      | Impedancia intrínseca e índice de refracción del vacío (= 1)       |
| `eps_r(i)`        | Permitividad relativa del slab `i`                                   |
| `n(i)`            | Índice de refracción del slab `i` = `sqrt(eps_r(i))`                |
| `d(i)`            | Espesor físico del slab `i` [m]                                     |
| `eta(i)`          | Impedancia intrínseca del slab `i` = `eta0/n(i)`                    |
| `theta1` / `th1`  | Ángulo de incidencia (en vacío), medido desde la normal             |
| `thetas(i)`       | Ángulo de propagación en el medio `i-1` (Snell)                     |
| `k_z(i)`          | Número de onda normal (eje z) dentro del slab `i`                   |
| `Z_TM(i)`         | Impedancia de onda TM dentro del slab `i` = `eta(i)*cos(theta_{i+1})` |
| `M_i`             | Matriz ABCD del slab `i`                                            |
| `M_tot` (`A,B,C,D`)| Matriz ABCD acumulada de todo el sistema                            |
| `Z_in`, `Z_out`   | Impedancias TM de los medios de entrada/salida (vacío) = iguales   |
| `Gamma_TM`        | Coeficiente de reflexión complejo                                   |
| `Tau_TM`          | Coeficiente de transmisión complejo                                 |
| `R_TM`, `T_TM`    | Reflectancia y transmitancia (potencia), `R+T=1`                    |

---

## 11. Limitaciones y alcance del modelo

- Válido únicamente para medios **sin pérdidas** (`eps_r` real, sin parte imaginaria).
  Si se introdujeran pérdidas (`eps_r` complejo), `det(M_i) ≠ 1` y `R+T < 1`
  (la diferencia representaría la potencia disipada).
- Válido para medios **no magnéticos** (`mu_r = 1` en todos los slabs).
- Asume medios de entrada y salida idénticos (vacío), lo que permite la simplificación
  `Z_in = Z_out`. Si se generalizara a medios de entrada/salida distintos, habría que
  recalcular `theta_salida` mediante Snell hasta el medio N+1 y usar las fórmulas
  generales de la Sección 7-8 sin la simplificación `Z_out = Z_in`.
- Polarización TM únicamente. Para TE, bastaría cambiar `Z_TE(i) = eta(i)/cos(theta_{i+1})`
  manteniendo el resto del formalismo ABCD idéntico.
