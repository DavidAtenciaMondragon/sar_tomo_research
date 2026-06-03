# SAR Image Formation Toolbox for MATLAB

**Autores:** LeRoy A. Gorham, Linda J. Moore  
**Año:** 2010  
**Fuente/Publicación:** Proc. SPIE Vol. 7699, 769906. Algorithms for Synthetic Aperture Radar Imagery XVII. Air Force Research Laboratory, WPAFB, OH.

---

## 1. Geometría del Sistema

Sistema SAR monoestático con trayectoria **arbitraria en 3D** (spotlight genérico). La posición de la antena en cada pulso es $\underline{r}_a(\tau) = [x_a(\tau), y_a(\tau), z_a(\tau)]^T$. Los píxeles de imagen también son 3D: $\underline{r} = [x, y, z]^T$. El paper presenta dos algoritmos implementados en MATLAB:

1. **Matched Filter (MF)** — $\mathcal{O}(N^4)$ para imagen 2D
2. **Backprojection (BP)** — $\mathcal{O}(N^3)$ para imagen 2D

Aplicado a 4 datasets públicos de AFRL: Backhoe Data Dome, 2D/3D Volumetric Challenge, SAR-GMTI Challenge, Civilian Vehicle Radar Data Domes.

---

## 2. Ecuaciones de Resolución SAR

### Posición de la antena

$$\underline{r}_a(\tau) = [x_a(\tau),\; y_a(\tau),\; z_a(\tau)]^T$$

**Variables:**
- `\tau` — apertura sintética (slow time)

### Distancia antena–escena origen

$$d_a(\tau) = \sqrt{x_a^2(\tau) + y_a^2(\tau) + z_a^2(\tau)}$$

### Distancia antena–objetivo

$$d_{a_0}(\tau) = \sqrt{(x_a(\tau)-x)^2 + (y_a(\tau)-y)^2 + (z_a(\tau)-z)^2}$$

**Variables:**
- `(x, y, z)` — posición del objetivo (o del píxel de imagen)

### Rango diferencial

$$\Delta R(\tau_n) = d_{a_0}(\tau_n) - d_a(\tau_n)$$

**Variables:**
- `\tau_n` — instante del n-ésimo pulso
- `d_a(\tau_n)` — distancia al punto de compensación de movimiento (scene origin)

### Señal Phase History (SAR data)

$$S(f_k, \tau_n) = A(f_k, \tau_n) \exp\!\left(\frac{-j4\pi f_k \Delta R(\tau_n)}{c}\right)$$

**Variables:**
- `f_k` — frecuencia del k-ésimo sample ($k = 1, \ldots, K$)
- `A(f_k, \tau_n)` — amplitud (relacionada con la RCS del blanco)
- `c` — velocidad de la luz

### Tamaño máximo de escena (alias-free)

$$W_r = \frac{c}{2\Delta f}, \qquad W_x = \frac{\lambda_{min}}{2\Delta\theta}$$

**Variables:**
- `\Delta f` — paso en frecuencia (Hz)
- `\lambda_{min} = c/f_K` — longitud de onda mínima
- `\Delta\theta` — paso angular en azimut (rad)

### Resolución en rango y en rango cruzado

$$\delta_r = \frac{c}{2B} = \frac{c}{2(K-1)\Delta f}$$

$$\delta_x = \frac{\lambda_c}{2\theta_a} = \frac{\lambda_c}{2(N_p-1)\Delta\theta}$$

**Variables:**
- `B = (K-1)\Delta f` — ancho de banda total
- `\lambda_c = c/f_c` — longitud de onda en la frecuencia central
- `\theta_a = (N_p-1)\Delta\theta` — ángulo de apertura total
- `N_p` — número de pulsos

### Matched Filter (imagen a una posición $\underline{r}$)

$$I(\underline{r}) = \frac{1}{N_p K} \sum_{n=1}^{N_p} \sum_{k=1}^{K} S(f_k, \tau_n) \exp\!\left(\frac{+j4\pi f_k \Delta R(\tau_n)}{c}\right)$$

### Perfil de rango (para BP, usando IFFT)

$$s(m, \tau_n) = K \cdot \text{fftshift}\!\left\{\text{ifft}(S(f_k, \tau_n))\right\} \cdot \exp\!\left(\frac{j2\pi f_1(m-1)}{N_{fft}\Delta f}\right)$$

**Variables:**
- `m` — índice del bin de rango ($m = 1 \ldots M$, con $m=1$ en zero-frequency tras fftshift)
- `N_{fft}` — longitud del IFFT (zero-padding factor, recomendado $N_{fft} = 10K$)
- `f_1` — frecuencia mínima del pulso

### Backprojection (imagen final)

$$I(\underline{r}) = \sum_{n=1}^{N_p} s_{int}(\underline{r}, \tau_n)$$

donde $s_{int}(\underline{r}, \tau_n)$ es el valor interpolado del perfil de rango en $\Delta R(\tau_n)$ con corrección de fase:

$$s_{int}(\underline{r},\tau_n) = \text{interp}(s(m,\tau_n),\, \Delta R(\tau_n)) \cdot \exp\!\left(\frac{+j4\pi f_1 \Delta R(\tau_n)}{c}\right)$$

---

## 3. Suposiciones del Modelo

1. El objetivo es estacionario durante toda la colección de datos.
2. Cada pulso está compensado en movimiento respecto al origen de la escena (zero-phase en el scene origin).
3. La amplitud $A(f_k, \tau_n)$ es constante (blanco isotrópico puntual) en el MF; variable en el BP general.
4. La interpolación usa `interp1(...,'linear')` de MATLAB (primer orden); para mejores resultados se recomienda sinc (zero-padding del IFFT).
5. El IFFT eficiente requiere $N_{fft}$ potencia de 2; recomendado $N_{fft} = 10K$.
6. La escena es plana ($z = 0$) en los ejemplos dados; el BP permite $z \neq 0$ si se dispone de DEM.
7. No se incluye compensación de la curvatura de frente de onda (adecuado para geometrías no demasiado cercanas).

---

## 4. Notas Adicionales

- El código MATLAB completo (funciones `mfBasic` y `bpBasic`) se incluye en los apéndices A.1 y A.2 del paper.
- El paper es de uso didáctico: orientado a investigadores SAR que no son expertos en procesado.
- Complejidades: MF → $\mathcal{O}(N^4)$; BP → $\mathcal{O}(N^3)$; Polar Format → $\mathcal{O}(N^2 \log N)$.
- Datasets aplicados: Backhoe (10 GHz, BW 5.9 GHz), 2D/3D Volumetric (100×100 m), GMTI (circular SAR), Civilian Vehicles (9.6 GHz, BW 5.35 GHz).
- Resoluciones de ejemplo: dataset 2D/3D → $\delta r = 0.24$ m, $\delta x = 0.23$ m.
