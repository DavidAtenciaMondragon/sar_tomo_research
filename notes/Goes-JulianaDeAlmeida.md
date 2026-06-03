# Techniques for High-Resolution 3D Images with Synthetic Aperture Radar

**Autora:** Juliana de Almeida Góes  
**Año:** 2022  
**Tipo:** Tesis Doctoral — Universidade Estadual de Campinas (UNICAMP), Faculdade de Engenharia Elétrica e de Computação  
**Orientador:** Prof. Dr. Hugo Enrique Hernandez Figueroa  
**Co-orientador:** Dr. Leonardo Sant'Anna Bins (INPE)  
**Banda:** P-band (λ = 70.54 cm, chirp BW = 50 MHz)  
**Sistema:** Drone-borne SAR (dron multirrotor, autonomía ~20 min, velocidad máx. 18 m/s)

> **Enfoque de este resumen:** El análisis está centrado en el **Capítulo 4 — Spiral SAR**, que contiene la contribución original de la tesis sobre resolución 3D para trayectorias espirales (cilíndricas y cónicas) con drones. El Capítulo 3 (algoritmo FFBP) se menciona brevemente al final.

---

## 1. Geometría del Sistema

### Descripción general
El sistema SAR opera desde un dron en una **trayectoria espiral (helicoidal)** alrededor de la escena. El sistema es **monoestático** (misma antena transmite y recibe). Dos variantes principales:

- **Espiral cilíndrica** (estado del arte): radio constante $\rho_0$, altura varía linealmente. Ángulo de inclinación $\beta = 90°$.
- **Espiral cónica** (propuesta nueva): radio varía linealmente desde $\rho_{top}$ (cima) hasta $\rho_{base}$ (base), altura decrece. Ángulo de inclinación $\beta < 90°$.

La geometría se basa en SAR Tomografía adaptada para apertura tomográfica inclinada (no vertical). La relación entre los parámetros clave es:

| Símbolo | Descripción |
|---------|-------------|
| $z_0$ | Altura media de la espiral |
| $\rho_0$ | Radio medio de la espiral |
| $\psi_0 = \tan^{-1}(\rho_0/z_0)$ | Ángulo de incidencia medio |
| $R_0 = \sqrt{z_0^2+\rho_0^2}$ | Distancia media radar–objetivo |
| $B$ | Apertura tomográfica total (longitud del vector espiral) |
| $B_\perp = B|\cos(\beta-\psi_0)|$ | Apertura tomográfica efectiva (⊥ a la LOS) |
| $\beta$ | Ángulo de inclinación de la espiral |
| $N_t$ | Número de vueltas |

---

## 2. Capítulo 4 — SPIRAL SAR: Análisis Completo de Resolución 3D

### 4.1 Contexto: Geometrías SAR para Imágenes 3D

#### 4.1.1 SAR Interferometría
$$\Delta\varphi = \frac{4\pi}{\lambda}\Delta R \tag{4.1}$$

**Variables:** $\lambda$ — longitud de onda; $\Delta R$ — diferencia de camino entre dos antenas interferométricas.

#### 4.1.2 SAR Tomografía
Resolución y máxima anchura ⊥ a la LOS:
$$\delta_{LOS_\perp} = \frac{\lambda R_0}{2B} \tag{4.2}$$
$$H_{LOS_\perp} \leq \frac{\lambda R_0}{2\Delta B} \tag{4.3}$$

**Variables:** $B$ — apertura tomográfica; $\Delta B$ — distancia de muestreo entre pistas; $R_0$ — distancia más corta radar–objetivo.

#### 4.1.3 Circular SAR
Resoluciones para objetivo isotrópico **en el eje** del círculo (Ishimaru et al. 1998):
$$\delta_{xy} = \frac{1.12\lambda}{2\pi\sin\psi} \tag{4.4}$$
$$\delta_z = \sqrt{\frac{\ln(2)}{\pi}\frac{c}{W\cos\psi}} \tag{4.5}$$

**Variables:** $\psi$ — ángulo de incidencia; $W$ — ancho de banda efectivo (tras compresión en rango); $c$ — velocidad de la luz. Nota: para $\psi = 45°$, $\delta_{xy} \approx \lambda/4$.

#### 4.1.4 Multi-Circular SAR
Combina SAR Circular + SAR Tomografía. Misma resolución en suelo (Ec. 4.4). La diversidad en ángulo de elevación amplía el ancho de banda efectivo vertical, mejorando $\delta_z$ respecto al SAR circular simple.

#### 4.1.5 Spiral SAR
La espiral cilíndrica resuelve el problema de distancia de muestreo crítica de Multi-Circular SAR (muestreo cuasi-continuo). La **espiral cónica** (nueva propuesta) combina la ventaja de la espiral con la apertura inclinada de SAR Tomografía para maximizar la resolución vertical.

---

### 4.2 El Desplazamiento de Número de Onda (Wavenumber Shift)

Concepto clave que explica la mejora de resolución al integrar múltiples perspectivas.

#### 4.2.1 Fundamentos del concepto

**Diferencia de ángulo de incidencia** entre dos posiciones de radar separadas $b_\perp$:
$$\Delta\psi \approx \frac{b_\perp}{R_0} \tag{4.6}$$

**Componentes del vector de onda** (terreno plano):
$$k_g = 2k\sin\psi, \quad -k_z = 2k\cos\psi \tag{4.7–4.8}$$

**Variables:** $k = 2\pi/\lambda$ — número de onda de señal; $k_g$, $k_z$ — componentes de suelo y vertical.

**Desplazamiento del número de onda en suelo:**
$$\Delta k_g \approx \frac{4\pi}{\lambda}\cos\psi\,\Delta\psi \approx \frac{4\pi b_\perp}{\lambda R_0}\cos\psi \tag{4.9–4.10}$$

**Desplazamiento espectral equivalente (que produce el mismo $\Delta k_g$):**
$$\Delta f = -\frac{cb_\perp}{\lambda R_0\tan\psi} \tag{4.12}$$

**Línea de base crítica** (decorrelación cuando los espectros no se solapan):
$$b_{\perp,c} = \frac{W\lambda}{c}R_0\tan\psi \tag{4.13}$$

**Resolución en suelo mejorada** al combinar dos espectros desplazados:
$$\delta_g = \frac{c}{2W_g}, \quad W_g = W\sin\psi + \Delta f_g, \quad \Delta f_g = \frac{cb_\perp}{\lambda R_0}\cos\psi \tag{4.14–4.17}$$

#### 4.2.1.5 Desplazamiento Vertical del Número de Onda — EXPRESIÓN REVISADA

> **Contribución original:** La literatura previa (Ponce et al. 2016) usaba $\partial\varphi/\partial z$ en lugar de $\Delta k_z$, produciendo valores incorrectos de $\delta_z$. La expresión correcta es:

$$\Delta k_z \approx \frac{4\pi}{\lambda}\sin\psi\,\Delta\psi \approx \frac{4\pi b_\perp}{\lambda R_0}\sin\psi \tag{4.18–4.19}$$

**Sensibilidad fase-altura** (≠ $\Delta k_z$):
$$\frac{\partial\varphi}{\partial z} = \frac{4\pi b_\perp}{\lambda R_0\sin\psi} \tag{4.21}$$

**Altura de ambigüedad:**
$$z_{2\pi} = \frac{\lambda R_0}{2b_\perp}\sin\psi \tag{4.22}$$

Relación entre sensibilidad y desplazamientos:
$$\frac{\partial\varphi}{\partial z} = \frac{\Delta k_g}{\tan\psi} + \Delta k_z \tag{4.23}$$

#### 4.2.2 Wavenumber Shift en Multi-Circular y Spiral SAR

Reemplazando la línea de base $b_\perp$ por la apertura tomográfica efectiva total $B_\perp$:

$$\Delta k_z(B_\perp) \approx \frac{4\pi B_\perp}{\lambda R_0}\sin\psi \tag{4.24}$$

**Ancho de banda efectivo vertical y su desplazamiento espectral:**
$$W_z = W\cos\psi + \Delta f_z, \quad \Delta f_z(B_\perp) = \frac{cB_\perp}{\lambda R_0}\sin\psi \tag{4.25–4.26}$$

### ★ ECUACIÓN CENTRAL: Resolución Vertical de Spiral SAR

$$\boxed{\delta_z = \sqrt{\frac{\ln(2)}{\pi}\frac{c}{W_z}}} \tag{4.27}$$

donde

$$W_z = W\cos\psi_0 + \frac{cB_\perp}{\lambda R_0}\sin\psi_0$$

**Variables:**
- `δ_z` — resolución vertical (−3 dB, half-power) [m]
- `W_z` — ancho de banda efectivo en la dirección vertical [Hz]
- `W` — ancho de banda del pulso tras compresión en rango [Hz]
- `ψ₀` — ángulo de incidencia medio de la espiral
- `B_⊥` — apertura tomográfica efectiva (⊥ a la LOS) [m]
- `R₀` — distancia media radar–objetivo [m]
- `λ` — longitud de onda [m]
- `c` — velocidad de la luz [m/s]

> **Interpretación física:** $\delta_z$ disminuye (mejora) al aumentar $B_\perp$. La espiral cónica con $\beta = \psi_0$ maximiza $B_\perp = B$ para un $B$ dado, dando la mejor resolución vertical posible con esa trayectoria.

**Distancia de muestreo crítica** (decorrelación si se supera):
$$\Delta B_{\perp,c} = \begin{cases} \dfrac{W\lambda}{c}R_0\tan\psi, & 0 \leq \psi < 45° \\[6pt] \dfrac{W\lambda}{c}\dfrac{R_0}{\tan\psi}, & 45° \leq \psi < 90° \end{cases} \tag{4.31}$$

**Altura de ambigüedad** (Spiral SAR con $N_t$ vueltas):
$$z_{2\pi} = \frac{\lambda R_0}{2\Delta B_\perp}\sin\psi \tag{4.32}$$
$$z_{2\pi} = \frac{N_t\lambda R_0}{2B_\perp}\sin\psi \tag{4.33}$$

> **Ventaja del Spiral SAR:** para objetivos isotrópicos, el muestreo es cuasi-continuo → $k_z$ varía de forma continua → decorrelación prácticamente inexistente, a diferencia de Multi-Circular SAR donde existen gaps en $k_z$.

---

### 4.3 Coordenadas y Parámetros de la Espiral Cónica

#### 4.3.1 Parametrización de la Trayectoria

Radio varía linealmente con velocidad radial $V_\rho$:
$$\rho(t) = \rho_{top} + V_\rho t, \quad V_\rho = \frac{\rho_{base}-\rho_{top}}{t_{max}} \tag{4.34–4.35}$$

Ángulo de azimut con velocidad tangencial constante $V_0$:
$$\alpha(t) = \frac{V_0}{V_\rho}\left[\ln(\rho(t)) - \ln(\rho_{top})\right] \tag{4.37}$$

Tiempo total de vuelo:
$$t_{max} = \frac{2\pi N_t}{V_0}\cdot\frac{\rho_{base}-\rho_{top}}{\ln(\rho_{base})-\ln(\rho_{top})} \tag{4.38}$$

**Caso cilíndrico** ($\rho_{base}=\rho_{top}=\rho_0$):
$$\alpha(t) = \frac{V_0}{\rho_0}t, \quad t_{max} = \frac{2\pi N_t\rho_0}{V_0} \tag{4.39–4.40}$$

Altura y coordenadas Cartesianas:
$$z(t) = z_{top} + V_z t, \quad V_z = -\frac{z_{top}-z_{base}}{t_{max}} \quad (V_z<0) \tag{4.41–4.42}$$
$$x(t) = \rho(t)\cos(\alpha(t)), \quad y(t) = \rho(t)\sin(\alpha(t)) \tag{4.43–4.44}$$

Parámetros derivados de diseño:
$$z_0 = \frac{z_{base}+z_{top}}{2}, \quad \rho_0 = \frac{\rho_{base}+\rho_{top}}{2} \tag{4.45–4.46}$$
$$B = \sqrt{(z_{top}-z_{base})^2+(\rho_{top}-\rho_{base})^2} \tag{4.47}$$
$$\beta = \tan^{-1}\!\left(\frac{z_{top}-z_{base}}{\rho_{base}-\rho_{top}}\right) \tag{4.48}$$
$$\psi_0 = \tan^{-1}\!\left(\frac{\rho_0}{z_0}\right), \quad R_0 = \sqrt{z_0^2+\rho_0^2} \tag{4.49–4.50}$$

**Variables:**
- `ρ_top`, `ρ_base` — radios en cima y base [m]
- `z_top`, `z_base` — alturas en cima y base [m]; el dron desciende ($z_{top} > z_{base}$)
- `V_ρ` — velocidad radial [m/s]
- `V_z` — velocidad vertical [m/s] (negativa)
- `V₀` — velocidad tangencial constante [m/s]
- `B` — apertura tomográfica = hipotenusa del triángulo ($\Delta z$, $\Delta\rho$) [m]
- `β` — ángulo de inclinación: $\beta=90°$ para cilíndrica, $\beta<90°$ para cónica

#### 4.3.2 Apertura Tomográfica Efectiva

$$B_\perp = B\,|\cos(\beta-\psi_0)| \tag{4.51}$$

**Variables:** $\beta$ — ángulo de inclinación; $\psi_0$ — ángulo de incidencia medio.

| Caso | $\beta$ | $B_\perp$ |
|------|---------|-----------|
| Cónica óptima | $\beta = \psi_0$ | $B$ (máximo) |
| Cilíndrica | $\beta = 90°$ | $B\sin\psi_0 < B$ |
| Cónica alineada con LOS | $\beta = \psi_0+90°$ | 0 (mínimo) |

**Corrección para objetivos fuera del eje** (efecto near-range):
$$\tilde{\psi}_0 = \tan^{-1}\!\left(\frac{\rho_0-|x_t|}{z_0}\right), \quad \widetilde{W} = \frac{W_z}{\cos(\psi_0)} \tag{4.52–4.53}$$

**Variables:** $x_t$ — posición del objetivo fuera del eje; $\widetilde{W}$ — ancho de banda equivalente corregido.

---

### 4.3.3 Resultados de Simulación: Efecto del Ángulo de Inclinación

**Parámetros de simulación (Tabla 4.1, $B_\perp$ constante):**

| Parámetro | Valor | Unidad |
|-----------|-------|--------|
| $B_\perp$ | 50 | m |
| $z_0$ | 95 | m |
| $\rho_0$ | 135 | m |
| $\psi_0$ | ~54.8° | — |
| $N_t$ | 10 | — |
| $V_0$ | 7.5 | m/s |
| $\lambda$ | 70.54 | cm |
| $W$ (tras comp.) | 20 | MHz |

**Conclusiones de simulación ($B_\perp$ fijo):**
- $\delta_x$, $\delta_y \approx 15$–18 cm son constantes con $\beta$ → solo dependen de $\lambda$ y $\psi_0$.
- $\delta_z \approx 1.2$ m es constante para todos los $\beta$ cuando $B_\perp$ se mantiene fijo. Valida Ec. (4.27).

**Simulación con $B$ fijo = 50 m, $\beta$ varía de 0° a 175° (Tabla 4.3):**

| $\beta$ | Geometría | $\delta_z$ (centro, sim.) |
|---------|-----------|--------------------------|
| 55° ($\approx\psi_0$) | Cónica óptima, B⊥LOS | ~1.25 m |
| 90° | Cilíndrica | ~2.0 m |
| 145° ($\psi_0+90°$) | Cónica, B ∥ LOS | >7 m |

**Spiral SAR vs. Multi-Circular SAR (Tabla 4.6):**

| $N_t$ | Tiempo Spiral | Tiempo Multi-Circ. | Diferencia |
|-------|--------------|-------------------|------------|
| 2 | 3 min 45 s | 5 min 39 s | −1 min 54 s |
| 4 | 7 min 31 s | 9 min 25 s | −1 min 54 s |
| 10 | 18 min 47 s | 20 min 44 s | −1 min 57 s |

- Spiral SAR converge más rápido: con $N_t = 4$ la imagen es casi idéntica a $N_t = 10$.
- En P-band con $R_0 \approx 165$ m, $\psi = 55°$: $\Delta B_{\perp,c} = 5.4$ m → Multi-Circular necesita pistas muy densas, problemático para la autonomía del dron.

---

### 4.4 Diseño de Trayectorias Espirales y Validación con Datos Reales

#### 4.4.1 Restricciones de Diseño

**Radio de iluminación constante** $\rho_{ci}$ (área siempre dentro del haz de la antena):

Condición 1 (base del vuelo, radio máximo + altura mínima):
$$\rho_{ci}^{base} \leq \frac{z_{min}}{\tan(\theta_{far})} - \rho_{max}, \quad \theta_{far} = \theta_{axis} - \frac{\theta_{3dB}}{2} \tag{4.54–4.55}$$

Condición 2 (cima del vuelo, radio mínimo + altura máxima):
$$\rho_{ci}^{top} \leq \rho_{min} - z_{max}\tan(\theta_{near}), \quad \theta_{near} = 90° - \theta_{axis} - \frac{\theta_{3dB}}{2} \tag{4.56–4.57}$$

$$\rho_{ci} = \min\!\left(\rho_{ci}^{base},\;\rho_{ci}^{top}\right) \tag{4.58}$$

**Variables:**
- `θ_axis` — ángulo de depresión del eje de la antena
- `θ_3dB` — ancho de haz de elevación (−3 dB)
- `θ_far`, `θ_near` — ángulos de depresión en far-range y near-range

#### 4.4.2 Trade-Off: $\delta_z$ vs. Área de Cobertura

Variando $\Delta\rho = \rho_{max} - \rho_{min}$ con alturas fijas:
$$\rho_{min} = \rho_{max} - \Delta\rho \tag{4.59}$$

Al aumentar $\Delta\rho$: $B\uparrow$ → $\delta_z$ mejora, pero $\rho_{ci}\downarrow$ (menor área cubierta). Para el sistema P-band con $z_{max}=114$ m, $z_{min}=84$ m, $\rho_{max}=118.5$ m:

| Trayectoria | $\beta$ | $\delta_z$ teórico | $\rho_{ci}$ | $B$ |
|-------------|---------|-------------------|------------|-----|
| Cilíndrica | 90° | 1.98 m | 72 m | 30 m |
| Cónica ($\beta=55°$) | 55° | 1.52 m | 51 m | 36.6 m |
| Cónica ($\beta=30°$) | 30° | — | ~21 m | 60 m |

#### 4.4.3 Campaña de Vuelo Real

- **Lugar:** Paulínia, São Paulo, Brasil — 30 julio 2021
- **Blanco:** Reflector de cuatro esquinas (quad-corner, 50×50 cm, 4 triedros)
- **Sistema:** Drone-borne SAR P-band

**Parámetros de vuelo reales (Tabla 4.7):**

| Parámetro | Cilíndrica | Cónica | Unidad |
|-----------|-----------|--------|--------|
| Radio mínimo | 112.6 | 98.1 | m |
| Radio máximo | 118.2 | 118.8 | m |
| Altura mínima | 83.2 | 84.4 | m |
| Altura máxima | 114.3 | 113.0 | m |
| $N_t$ | 4 | 4 | — |
| Velocidad media | 6.94 | 6.85 | m/s |
| Posiciones radar | 67,355 | 63,457 | — |

**Procesado FFBP** — volumen de imagen: $6.4\times6.4\times6.4$ m³, resolución $10\times10\times10$ cm³, centrado en $(-2, 25.5, 0)$ m.

**Ecualización de datos** (compensación de modulación angular del reflector):
1. Extraer señal del quad-corner para cada posición del radar.
2. Filtrar con dos medias móviles (ventanas 31 y 551 píxeles).
3. Dividir cada columna de la matriz SAR por la señal filtrada.

#### 4.4.5 Resultados Finales — Resolución 3D Medida

**Tabla 4.10 — Comparación valores teóricos vs. reales:**

| Resolución | Cil. teórico (0,0) | Cil. actual | Cón. teórico (0,0) | Cón. actual |
|-----------|-------------------|--------------|--------------------|-------------|
| $\delta_x$ [cm] | 16.6 | 16.6 | 17.1 | 17.3 |
| $\delta_y$ [cm] | 18.8 | 17.4 | 19.6 | 17.1 |
| $\delta_z$ [m] | **1.98** | **2.33** | **1.52** | **1.76** |

**Conclusión principal:**
- Espiral cónica ($\beta \approx 55°$): $\delta_z = \mathbf{1.76}$ m  
- Espiral cilíndrica ($\beta = 90°$): $\delta_z = \mathbf{2.33}$ m  
- Mejora: **factor 1.32×** con datos reales (1.30× predicho por la teoría — error < 1.5%).
- Ambas resoluciones en suelo son prácticamente idénticas: $\delta_x \approx \delta_y \approx 17$ cm.

---

## 3. Suposiciones del Modelo

1. **Objetivo isotrópico** para el cálculo teórico de resolución; objetivos no isotrópicos requieren ecualización de datos.
2. **Far-field:** $b_\perp \ll R_0$ para la linealización del desplazamiento de número de onda (Ec. 4.6).
3. **Terreno plano** para las expresiones de $k_g$ y $k_z$ (Ecs. 4.7–4.8); en la práctica se usa un DTM refinado.
4. **Velocidad tangencial constante** $V_0$; la velocidad radial y vertical son mucho menores y se desprecian en las Ecs. 4.36–4.40.
5. **Trayectoria espiral suave** con variación lineal de radio y altura; en la práctica el dron usa waypoints discretos.
6. **Objetivo en el eje** para las expresiones de Circular SAR (Ecs. 4.4, 4.5); para objetivos fuera del centro se aplica la corrección de ángulo equivalente $\tilde{\psi}_0$ (Ec. 4.52).
7. **Ancho de banda fraccional bajo** ($W\lambda/c \ll 1$): justifica la aproximación de primer orden del desplazamiento de número de onda.
8. El sistema es **monoestático** en todos los análisis.
9. La resolución en suelo $\delta_{xy}$ es **independiente del ancho de banda** en geometrías circulares/espirales; depende solo de $\lambda$ y $\psi_0$ (Ec. 4.4).
10. El número de vueltas $N_t$ **no afecta** $\delta_z$ en Spiral SAR (solo afecta $z_{2\pi}$ y PSLR).

---

## 4. Síntesis de Ecuaciones de Resolución 3D — Spiral SAR

### Resolución en plano (ground)
$$\delta_{xy} = \frac{1.12\lambda}{2\pi\sin\psi_0}$$

### Resolución vertical
$$\delta_z = \sqrt{\frac{\ln(2)}{\pi}\frac{c}{W_z}}, \quad W_z = W\cos\psi_0 + \frac{cB_\perp}{\lambda R_0}\sin\psi_0$$

### Apertura tomográfica efectiva
$$B_\perp = B\,|\cos(\beta-\psi_0)|$$

### Condición óptima ($\delta_z$ mínimo)
$$\beta = \psi_0 \implies B_\perp = B \implies \delta_z^{min} = \sqrt{\frac{\ln(2)}{\pi}\frac{c}{W\cos\psi_0 + \frac{cB}{\lambda R_0}\sin\psi_0}}$$

### Altura de ambigüedad
$$z_{2\pi} = \frac{N_t\lambda R_0}{2B_\perp}\sin\psi_0$$

---

## 5. Notas sobre el Algoritmo FFBP (Capítulo 3)

El Capítulo 3 desarrolla un algoritmo **Fast Factorized Back-Projection (FFBP)** en coordenadas Cartesianas 3D. Sus características originales:
- Usa una **curva de Morton modificada** (Z-order 3D) para dividir la imagen en un octree flexible.
- Complejidad: $\mathcal{O}(P^2\log P)$ para 2D, $\mathcal{O}(P^3)$ para 3D.
- Speedup máximo: **21×** para imágenes 3D (12× promedio) vs. BP directo, manteniendo error de fase $\sigma_{\Delta\phi} < 4°$.
- El FFBP fue indispensable para generar las 172 imágenes 3D del análisis de trayectorias del Cap. 4.

---

## 6. Aplicaciones y Trabajo Futuro

- **Tomografía subsuperficial:** detección de objetos enterrados, estimación de propiedades del suelo por capas.
- **Inventario forestal y biomasa** con SAR P-band.
- **Minería:** mapeo de galerías subterráneas.
- Trabajo futuro: incorporar el paso de apertura sintética selectiva al FFBP para mejorar imágenes fuera de $\rho_{ci}$; repetir análisis en otras bandas de frecuencia (L, C); investigar decorrelación en escenarios de mayor $R_0$.
