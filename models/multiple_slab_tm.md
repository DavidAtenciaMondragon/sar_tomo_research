# Fundamento teórico: Reflexión y Transmisión TM en un Sistema Multicapa

Este documento describe el fundamento teórico, las ecuaciones y las premisas detrás de
`multislabTMcoef.m`, que calcula los coeficientes de reflexión (Γ) y transmisión (τ), así
como la reflectancia (R) y transmitancia (T), para una onda plana con polarización TM
(paralela) incidiendo de forma oblicua sobre un sistema de **N slabs dieléctricos
finitos** apilados entre un semiespacio de vacío (entrada) y un **semiespacio de salida
dieléctrico** (la última capa, semi-infinita en z−).

---

## 1. Configuración del sistema y premisas

```
   Medio 0        Medio 1     Medio 2   ...   Medio N-1       Medio N
   (vacío)        (Slab 1)    (Slab 2)        (Slab N-1)      (semi-inf.)
  n0=1, η0  →    n1,η1  →    n2,η2  →   ...  → n_{N-1} →   nN, ηN
  θ_in       θ1           θ2            θ_{N-1}           θ_out
       z=0       z=-d1     z=-(d1+d2)         z=-Σ_{i=1}^{N-1} d_i
```

**Premisas del modelo:**

1. Onda plana monocromática de frecuencia `f`, incidiendo desde un semiespacio de
   **vacío** (Medio 0: `n0 = 1`, `η0 = sqrt(μ0/ε0)`).
2. El sistema termina en un **semiespacio semi-infinito** de índice `n_out` (Medio N),
   que representa la capa de suelo que contiene al target. **No existe aire ni reflexión
   en la cara inferior**: toda la energía transmitida se propaga hacia z− sin regresar.
3. Entre el vacío y el semiespacio de salida existen `N−1` capas intermedias finitas
   (Slabs 1 … N−1), cada una con:
   - Permitividad relativa `eps_r(i)`, índice `n(i) = sqrt(eps_r(i))`
   - Espesor físico `d(i)`
   - Impedancia intrínseca `eta(i) = eta0 / n(i)`
4. Todos los medios son dieléctricos sin pérdidas (`eps_r` real) y no magnéticos
   (`μr = 1`).
5. Las interfaces son planas, infinitas y paralelas (perpendiculares al eje `z`).
6. Polarización **TM** (campo eléctrico contenido en el plano de incidencia).
7. Análisis en régimen permanente (fasores, dependencia temporal `e^{jωt}` implícita).

> **Relación con el script:** en `GS_snell_multislab_test_script.m`, si el target está en
> la capa L de un stack de M capas, entonces `n_slabs = n_layers(1:L-1)` y
> `n_out = n_layers(L)`. La capa L actúa como semiespacio semi-infinito de salida.

---

## 2. Geometría de propagación: Ley de Snell generalizada

La componente del vector de onda **paralela** a las interfaces (`k_x`) se conserva en
cada interfaz. Para los medios de la cadena `n0 → n1 → … → n_{N-1} → n_out`:

```
n0 * sin(θ_in) = n(i) * sin(θ_i)   para todo i = 1..N-1
n0 * sin(θ_in) = n_out * sin(θ_out)
```

En el código, con `n_seq = [n_in, n_slabs(1), …, n_slabs(N_slabs), n_out]`:

```matlab
for i = 1:(N_slabs+1)
    arg = n_seq(i)/n_seq(i+1) * sin(thetas(i));
    if abs(arg) > 1
        TIR = true; break;
    end
    thetas(i+1) = asin(arg);
end
```

`thetas(1) = θ_in` (en vacío), `thetas(end) = θ_out` (en el semiespacio de salida).

### Invariante global

```
n0 * sin(θ_in) = n_out * sin(θ_out)
```

Dado que `n_out > n0 = 1` en la aplicación de suelo, se tiene `θ_out < θ_in`:
el rayo se refracta hacia la normal al penetrar en un medio más denso.
**A diferencia del caso simétrico (vacío–vacío), θ_out ≠ θ_in y Z_in ≠ Z_out.**

### Reflexión interna total (TIR)

TIR ocurre sólo cuando `n_seq(i)/n_seq(i+1) > 1` (transición de medio más denso a
menos denso) y el ángulo supera el crítico: `arg > 1`. En el escenario habitual
(n0 = 1 < n_slabs ≤ n_out), las transiciones son siempre de menos denso a más denso,
por lo que TIR no se produce para ningún ángulo de incidencia real.

---

## 3. Impedancia de onda y número de onda transversal (TM)

Para cada slab intermedio `i` (ángulo de propagación `θ_i = thetas(i+1)`):

```
k_z(i) = k0 * n(i) * cos(θ_i)          k0 = 2πf/c
Z_TM(i) = eta(i) * cos(θ_i)            eta(i) = eta0 / n(i)
```

Para los medios extremos:

```
Z_in  = eta0            * cos(θ_in)     (vacío de entrada)
Z_out = (eta0 / n_out)  * cos(θ_out)   (semiespacio de salida)
```

**Z_in ≠ Z_out** en general, ya que `n_out > 1` y `θ_out ≠ θ_in`.

```matlab
Z_in  = (eta0/n_in)  * cos(thetas(1));
Z_out = (eta0/n_out) * cos(thetas(end));
```

---

## 4. Modelo de línea de transmisión equivalente (matriz ABCD)

Cada slab finito `i` (de los `N_slabs` intermedios) se modela como un tramo de línea
sin pérdidas con impedancia `Z_TM(i)` y longitud eléctrica `k_z(i) · d(i)`:

```
        ┌                                          ┐
M_i  =  │  cos(k_z d_i)          j Z_TM sin(k_z d_i)  │
        │  j sin(k_z d_i)/Z_TM    cos(k_z d_i)         │
        └                                          ┘
```

```matlab
M_i = [cos(k_z*d(i)),        1j*Z_TM*sin(k_z*d(i));
       1j/Z_TM*sin(k_z*d(i)), cos(k_z*d(i))       ];
```

Para un medio sin pérdidas, `det(M_i) = 1`.

---

## 5. Matriz ABCD total del sistema

```
M_tot = M_1 · M_2 · … · M_{N_slabs}  =  [ A  B ]
                                          [ C  D ]
```

```matlab
M_tot = eye(2,2);
for i = 1:N_slabs
    M_tot = M_tot * M_i;
end
A = M_tot(1,1);  B = M_tot(1,2);
C = M_tot(2,1);  D = M_tot(2,2);
```

`det(M_tot) = A*D − B*C = 1` (conservación de energía, todos los slabs sin pérdidas).

---

## 6. Impedancias de los medios de entrada y salida

La red de dos puertos `M_tot` está:
- Excitada desde el vacío de entrada con `Z_in = eta0 · cos(θ_in)`
- Cargada en el semiespacio de salida con `Z_out = (eta0/n_out) · cos(θ_out)`

Como `n_out > n_in = 1` y `θ_out < θ_in`, en general:

```
Z_in  = eta0 · cos(θ_in)          > Z_out  (el vacío tiene mayor impedancia TM)
Z_out = (eta0/n_out) · cos(θ_out) < Z_in
```

**No existe simplificación Z_in = Z_out.** Se usan las fórmulas generales de la
red de dos puertos en las Secciones 7 y 8.

---

## 7. Coeficiente de reflexión Γ_TM

La impedancia de entrada de la red cargada con `Z_out` es:

```
Z_red = (A · Z_out + B) / (C · Z_out + D)
```

El coeficiente de reflexión respecto al medio de entrada (`Z_in`):

```
Γ_TM = (Z_red − Z_in) / (Z_red + Z_in)
```

Sustituyendo `Z_red` y simplificando:

```
             A·Z_out + B − Z_in·(C·Z_out + D)
Γ_TM = ─────────────────────────────────────────
             A·Z_out + B + Z_in·(C·Z_out + D)
```

```matlab
num = A*Z_out + B - Z_in*(C*Z_out + D);
den = A*Z_out + B + Z_in*(C*Z_out + D);
Gamma_TM(idx) = num / den;
```

---

## 8. Coeficiente de transmisión τ_TM y conservación de energía

```
             2 · Z_out
τ_TM = ─────────────────────────────────────
        A·Z_in + B + C·Z_in·Z_out + D·Z_out
```

```matlab
Tau_TM(idx) = 2*Z_out / (A*Z_in + B + C*Z_in*Z_out + D*Z_out);
```

### Reflectancia y transmitancia (potencia)

```
R_TM = |Γ_TM|²
T_TM = |τ_TM|² · (Z_in / Z_out)
```

El factor `Z_in/Z_out` corrige la diferencia de impedancias entre los semiespacios
de entrada y salida: la potencia incidente es proporcional a `1/Z_in` y la transmitida
a `1/Z_out`.

```matlab
R_TM = abs(Gamma_TM).^2;
T_TM = abs(Tau_TM).^2 .* (Z_in/Z_out);
```

### Conservación de energía

Para medios sin pérdidas (`det(M_tot) = 1`, todos los `Z_i` reales):

```
R_TM + T_TM = 1     (para todo θ_in sin TIR)
```

Esta identidad es la verificación numérica del modelo. Se puede comprobar
analíticamente para el caso de interfaz única (N_slabs = 0, M_tot = I):

```
Γ = (Z_out − Z_in)/(Z_out + Z_in)
τ = 2·Z_out/(Z_in + Z_out)
R + T = (Z_out−Z_in)²/(Z_out+Z_in)² + 4·Z_in·Z_out/(Z_out+Z_in)² = 1  ✓
```

---

## 9. Resumen del algoritmo (correspondencia con el código)

Para cada ángulo de incidencia `θ_in`:

1. **Snell**: calcular `thetas(i)` para `i = 1..N_slabs+2` usando `n_seq`.
   Si `|n_seq(i)/n_seq(i+1) · sin(thetas(i))| > 1` → TIR (`Γ=1`, `τ=0`, fin).
2. **Por cada slab finito `i = 1..N_slabs`**:
   - `k_z(i) = k0 · n(i) · cos(thetas(i+1))`
   - `Z_TM(i) = (eta0/n(i)) · cos(thetas(i+1))`
   - Construir `M_i` (ABCD del tramo equivalente)
3. **Cascada**: `M_tot = M_1 * … * M_{N_slabs} = [A B; C D]`
4. **Impedancias de los extremos**:
   - `Z_in  = eta0 · cos(thetas(1))`
   - `Z_out = (eta0/n_out) · cos(thetas(end))`
5. **Coeficientes**:
   - `Γ_TM = (A·Z_out + B − Z_in·(C·Z_out + D)) / (A·Z_out + B + Z_in·(C·Z_out + D))`
   - `τ_TM = 2·Z_out / (A·Z_in + B + C·Z_in·Z_out + D·Z_out)`
6. **Potencias**: `R_TM = |Γ|²`, `T_TM = |τ|² · Z_in/Z_out`
7. **Verificación**: `R_TM + T_TM = 1`

---

## 10. Variables principales del script y su significado físico

| Variable            | Significado                                                            |
|---------------------|------------------------------------------------------------------------|
| `f`                 | Frecuencia de la onda incidente [Hz]                                   |
| `c`                 | Velocidad de la luz en el vacío [m/s]                                  |
| `eta0`              | Impedancia intrínseca del vacío ≈ 377 Ω                                |
| `n_in`              | Índice del medio de entrada (vacío = 1)                                |
| `n_slabs`           | Índices de refracción de los slabs finitos intermedios [1 × N_slabs]  |
| `d_slabs`           | Espesores de los slabs finitos [m] [1 × N_slabs]                      |
| `n_out`             | Índice del semiespacio de salida (capa del target, semi-infinita)      |
| `theta_in`          | Ángulo de incidencia en vacío, desde la normal [rad]                   |
| `thetas(i)`         | Ángulo de propagación en el i-ésimo medio de la cadena                 |
| `k_z(i)`            | Número de onda normal (eje z) en el slab `i`                           |
| `Z_TM(i)`           | Impedancia de onda TM del slab `i` = `(eta0/n(i))·cos(thetas(i+1))`   |
| `Z_in`              | Impedancia TM del vacío de entrada = `eta0·cos(θ_in)`                 |
| `Z_out`             | Impedancia TM del semiespacio de salida = `(eta0/n_out)·cos(θ_out)`   |
| `M_i`               | Matriz ABCD del slab `i`                                               |
| `M_tot` (`A,B,C,D`) | Matriz ABCD acumulada de todos los slabs finitos                       |
| `Gamma_TM`          | Coeficiente de reflexión complejo                                      |
| `Tau_TM`            | Coeficiente de transmisión complejo                                    |
| `R_TM`, `T_TM`      | Reflectancia y transmitancia de potencia; `R + T = 1` (sin pérdidas)  |

---

## 11. Limitaciones y alcance del modelo

- **Medios sin pérdidas**: `eps_r` debe ser real. Con pérdidas (`eps_r` complejo),
  `det(M_i) ≠ 1` y `R + T < 1` (la diferencia es la potencia disipada).
- **Medios no magnéticos**: `μr = 1` en todos los slabs.
- **Última capa semi-infinita**: el modelo asume que no existe reflexión en la cara
  inferior del último medio. `d_slabs` sólo contiene los espesores de las capas
  **por encima** del semiespacio de salida. Si el target está más profundo que el
  espesor declarado de la última capa en el JSON, `calculateSlantRangeMultislab`
  asigna el target a esa última capa sin error (comportamiento semi-infinito).
- **Polarización TM únicamente**. Para TE, basta cambiar
  `Z_TE(i) = eta(i)/cos(θ_i)` manteniendo el formalismo ABCD idéntico.
- **Caso simétrico (vacío–vacío)**: si `n_out = n_in = 1`, se recupera automáticamente
  `Z_out = Z_in`, `θ_out = θ_in`, y las fórmulas reducen a la forma simplificada
  con `A·Z_in + B − C·Z_in² − D·Z_in` en el numerador de Γ, consistente con la
  literatura estándar (e.g. Pozar, *Microwave Engineering*).
