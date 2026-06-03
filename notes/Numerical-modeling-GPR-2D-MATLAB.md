# Numerical Modeling of Ground-Penetrating Radar in 2-D Using MATLAB

**Autores:** James Irving, Rosemary Knight  
**Año:** 2006  
**Fuente/Publicación:** Computers & Geosciences, vol. 32, pp. 1247–1258. DOI: 10.1016/j.cageo.2005.11.006. Geophysics Department, Stanford University.

---

## 1. Geometría del Sistema

El paper presenta códigos MATLAB de modelado numérico FDTD (Finite-Difference Time-Domain) 2D de GPR para dos configuraciones geométricas:

1. **Reflexión en superficie (surface-based reflection GPR):** transmisor y receptor en la superficie terrestre ($z = 0$), barriendo a lo largo de una línea (eje x). Las antenas son perpendiculares al plano $x$-$z$ de la sección. Se usa la formulación TM (Transverse Magnetic) con campos $\{H_x, H_z, E_y\}$.

2. **Sondeo en pozo (crosshole / VRP):** transmisor y receptor en pozos separados, dentro del plano de la sección. Se usa la formulación TE (Transverse Electric) con campos $\{E_x, E_z, H_y\}$.

- **Grilla de modelado:** 2D en el plano $x$-$z$ (posición horizontal y profundidad). Las fuentes y receptores son líneas infinitas en la dirección perpendicular $y$ (no modelada).
- **Condiciones de contorno:** capas absorbentes PML (Perfectly Matched Layer) en todos los bordes para evitar reflexiones artificiales de los bordes de la grilla.
- **Parámetros eléctricos del modelo:** permitividad dieléctrica $\varepsilon$, permeabilidad magnética $\mu$, conductividad eléctrica $\sigma$; todos espacialmente variables.
- **Ejemplo de reflexión (Fig. 4):** modelo con capa vadosa superior ($\varepsilon_r = 9$, $\sigma = 1$ mS/m) y capa saturada inferior ($\varepsilon_r = 25$, $\sigma = 5$ mS/m), separadas por una interfaz inclinada; tres bloques anómalos en la capa superior ($\varepsilon_r = 16$, $\sigma = 1$ mS/m). Fuentes y receptores a lo largo de la superficie cada 0.2 m. Pulso Blackman-Harris con frecuencia dominante 100 MHz. Discretización: $\Delta x = \Delta z = 0.04$ m, $\Delta t = 0.08$ ns.
- **Ejemplo crosshole (Fig. 8):** modelo de velocidades de EM obtenido de permitividades relativas entre 20 y 32 en la zona saturada. Pozos de fuentes a $x = 0.5$ m (de 0.5 a 10.5 m de profundidad, cada 0.25 m); pozos de receptores a $x = 5.5$ m. Pulso Blackman-Harris 100 MHz. $\Delta x = \Delta z = 0.025$ m, $\Delta t = 0.02$ ns.

---

## 2. Ecuaciones de Resolución SAR/GPR

### Ecuaciones de Maxwell en el dominio frecuencia (espacio estirado/stretched)

$$\nabla \times \mathbf{E} = -i\omega\mu\mathbf{H}$$

$$\nabla \times \mathbf{H} = \sigma\mathbf{E} + i\omega\varepsilon\mathbf{E}$$

**Variables:**
- `E`, `H` — vectores de campo eléctrico y magnético
- `omega` — frecuencia angular ($\omega = 2\pi f$)
- `epsilon` — permitividad dieléctrica del medio [F/m]
- `mu` — permeabilidad magnética del medio [H/m]
- `sigma` — conductividad eléctrica del medio [S/m]

### Operador nabla en espacio de coordenadas estiradas (para PML)

$$\nabla = \hat{x}\frac{1}{s_x}\frac{\partial}{\partial x} + \hat{y}\frac{1}{s_y}\frac{\partial}{\partial y} + \hat{z}\frac{1}{s_z}\frac{\partial}{\partial z}$$

con variables de estiramiento complejas:

$$s_k = \kappa_k + \frac{\sigma_k}{\alpha_k + i\omega\varepsilon_0}, \quad k = x, y, z$$

**Variables:**
- `s_k` — variable de estiramiento compleja en la dirección k (igual a 1 en el interior de la grilla)
- `kappa_k` — parámetro real de estiramiento (igual a 1 en el interior)
- `sigma_k` — conductividad de la PML en la dirección k (igual a 0 en el interior)
- `alpha_k` — parámetro de amortiguamiento adicional de la PML
- `epsilon_0` — permitividad del espacio libre

### Ecuaciones de modo TM en dominio frecuencia (reflexión en superficie)

$$i\omega\mu H_x = -\frac{1}{s_z}\frac{\partial E_y}{\partial z}, \quad i\omega\mu H_z = \frac{1}{s_x}\frac{\partial E_y}{\partial x}$$

$$\sigma E_y + i\omega\varepsilon E_y = \frac{1}{s_x}\frac{\partial H_z}{\partial x} - \frac{1}{s_z}\frac{\partial H_x}{\partial z}$$

**Variables:**
- `H_x`, `H_z` — componentes del campo magnético en x y z
- `E_y` — componente del campo eléctrico en y (la única componente E para modo TM en 2D)

### Ecuaciones FDTD en dominio tiempo — modo TM (tras CPML)

$$\mu\frac{\partial H_x}{\partial t} = -\frac{1}{\kappa_z}\frac{\partial E_y}{\partial z} - \zeta_z(t) * \frac{\partial E_y}{\partial z}$$

$$\mu\frac{\partial H_z}{\partial t} = \frac{1}{\kappa_x}\frac{\partial E_y}{\partial x} + \zeta_x(t) * \frac{\partial E_y}{\partial x}$$

$$\varepsilon\frac{\partial E_y}{\partial t} = \frac{1}{\kappa_x}\frac{\partial H_z}{\partial x} - \frac{1}{\kappa_z}\frac{\partial H_x}{\partial z} + \zeta_x(t)*\frac{\partial H_z}{\partial x} - \zeta_z(t)*\frac{\partial H_x}{\partial z} - \sigma E_y$$

**Variables:**
- `*` — convolución temporal (implementada mediante la técnica CPML recursiva)
- `zeta_k(t)` — función de convolución de la PML en la dirección k

### Ecuaciones de actualización FDTD discretizadas (TM, O(2,4))

$$H_x\big|_{i,j+1/2}^{n+1/2} = H_x\big|_{i,j+1/2}^{n-1/2} - D_{b_x}\big|_{i,j+1/2}\left[-E_y\big|_{i,j+2}^n + 27E_y\big|_{i,j+1}^n - 27E_y\big|_{i,j}^n + E_y\big|_{i,j-1}^n\right] - D_e\big|_{i,j+1/2}\left[\Psi_{H_{xz}}\big|_{i,j+1/2}^{n-1/2}\right]$$

$$E_y\big|_{i,j}^{n+1} = C_a\big|_{i,j} E_y\big|_{i,j}^n + C_{b_x}\big|_{i,j}\left[-H_z\big|_{i+3/2,j}^{n+1/2} + 27H_z\big|_{i+1/2,j}^{n+1/2} - 27H_z\big|_{i-1/2,j}^{n+1/2} + H_z\big|_{i-3/2,j}^{n+1/2}\right] + C_{b_z}(\ldots) + C_c\big|_{i,j}\left[\Psi_{E_{yx}}\big|_{i,j}^n - \Psi_{E_{yz}}\big|_{i,j}^n\right]$$

**Variables (coeficientes de actualización):**

$$C_a = \left(1 - \frac{\sigma\Delta t}{2\varepsilon}\right)\left(1 + \frac{\sigma\Delta t}{2\varepsilon}\right)^{-1}$$

$$C_{b_k} = \frac{\Delta t}{\varepsilon}\left(1 + \frac{\sigma\Delta t}{2\varepsilon}\right)^{-1}(24\kappa_k\Delta k)^{-1}$$

$$C_c = \frac{\Delta t}{\varepsilon}\left(1 + \frac{\sigma\Delta t}{2\varepsilon}\right)^{-1}$$

$$D_{b_k} = \frac{\Delta t}{\mu}(24\kappa_k\Delta k)^{-1}, \qquad D_e = \frac{\Delta t}{\mu}$$

**Variables:**
- `Delta_t` — paso temporal de la simulación [s]
- `Delta_k` — tamaño de celda espacial en dirección k [m]
- `epsilon`, `mu`, `sigma` — propiedades eléctricas locales (varían espacialmente)
- `kappa_k` — parámetro de estiramiento PML (= 1 en el interior de la grilla)
- Los superíndices $n$ indican el paso temporal; los subíndices $i,j$ la posición espacial

### Criterio de estabilidad numérica (CFL para esquema O(2,4))

$$\Delta t_{\max} = \frac{6}{7}\sqrt{\frac{\mu_{\min}\varepsilon_{\min}}{1/\Delta x^2 + 1/\Delta z^2}}$$

**Variables:**
- `mu_min`, `epsilon_min` — valores mínimos de permeabilidad y permitividad presentes en la grilla
- `Delta_x`, `Delta_z` — tamaños de celda espacial en x y z

### Criterio de dispersión numérica (mínimo muestreo espacial)

Para el esquema O(2,4), se requieren al menos **5 muestras por longitud de onda mínima** en cada dirección:

$$\Delta x, \Delta z \leq \frac{\lambda_{\min}}{5} = \frac{v_{\min}}{5 f_{\max}}$$

**Variables:**
- `lambda_min` — longitud de onda mínima en el medio más lento de la grilla
- `v_min` — velocidad EM mínima en la grilla ($= c/\sqrt{\varepsilon_{r,\max}\mu_{r,\max}}$)
- `f_max` — frecuencia máxima de la señal fuente

### Perfil de conductividad PML (escala polinomial)

$$\kappa_k = \begin{cases} 1 & \text{interior de la grilla} \\ 1 + \left(\frac{d}{\delta}\right)^m (\kappa_{k_{\max}} - 1) & \text{región PML} \end{cases}$$

$$\sigma_k = \begin{cases} 0 & \text{interior de la grilla} \\ \left(\frac{d}{\delta}\right)^m \sigma_{k_{\max}} & \text{región PML} \end{cases}$$

$$\sigma_{k_{\max}} = \frac{m+1}{150\pi\sqrt{\varepsilon_r}\Delta k}$$

**Variables:**
- `d` — distancia dentro de la región PML desde la frontera interior/PML [m]
- `delta` — espesor total de la región PML [m]
- `m` — exponente polinomial de la PML (valor por defecto: $m = 4$)
- `kappa_{k_max}` — valor máximo de $\kappa_k$ (por defecto: $\kappa_{k_{\max}} = 5$)
- `epsilon_r` — permitividad relativa del material en el borde interior de la grilla adyacente a la PML

### Velocidad de propagación EM (aproximación de bajas pérdidas)

$$v = \frac{c}{\sqrt{\varepsilon_r}}$$

**Variables:**
- `v` — velocidad de propagación de ondas EM en el medio [m/s]
- `c` — velocidad de la luz en el vacío ($\approx 0.3$ m/ns)
- `epsilon_r` — permitividad relativa del medio (asumiendo $\mu_r = 1$ y bajas pérdidas)

---

## 3. Suposiciones del Modelo

1. El modelo es estrictamente 2D: todas las fuentes y receptores son elementos de línea infinitos en la dirección y (perpendicular al plano de la sección). Los diagramas de radiación y la propagación geométrica (spreading esférico vs. cilíndrico) difieren del caso 3D real.
2. Las propiedades eléctricas del modelo ($\varepsilon$, $\mu$, $\sigma$) son independientes de la frecuencia (no dispersivas). Para modelar dispersión dieléctrica (e.g., agua) se requeriría modificar el código con modelos como Debye o Cole-Cole.
3. Las antenas se modelan como fuentes de campo $E_y$ (TM) o $H_y$ (TE) en un único punto de la grilla; no se modela la geometría real de la antena ni su diagrama de radiación.
4. El pulso fuente normalizado es la primera derivada de una función de ventana Blackman-Harris, que al propagarse por la grilla produce una forma similar a un pulso de Ricker.
5. La aproximación semi-implícita del término de corriente de conducción en la ecuación de $E_y$ (Eq. 8c) garantiza propiedades numéricas superiores a las expresiones unilaterales para medios con pérdidas.
6. El esquema leap-frog staggered-grid (Yee, 1966) actualiza alternadamente los campos eléctrico y magnético, desfasados $\Delta t/2$ en tiempo y $\Delta k/2$ en espacio.
7. Las matrices de propiedades eléctricas tienen el doble del tamaño de las matrices de campos (porque $\varepsilon$, $\mu$, $\sigma$ se requieren en cada ubicación de componente de campo), y se interpolan a la resolución de la mitad del paso espacial antes de la simulación.

---

## 4. Notas Adicionales

- **Códigos MATLAB disponibles:** `TM_model2d.m` (reflexión en superficie, modo TM); `TE_model2d.m` (crosshole/VRP, modo TE). Scripts de ejemplo: `TM_run_example.m`, `TE_run_example.m`. Utilidades auxiliares: `finddx.m` (determina $\Delta x$, $\Delta z$ máximos), `finddt.m` (determina $\Delta t$ máximo), `gridinterp.m` (interpola propiedades a la resolución de medio paso), `padgrid.m` (añade celdas PML), `blackharrispulse.m` (genera el pulso fuente). Disponibles en: http://www.iamg.org/CGEditor/index.htm.
- **Modo TM vs TE:** TM se usa para reflexión en superficie porque las antenas (dipolos) son perpendiculares al plano de la sección $x$-$z$; TE se usa para crosshole porque las antenas están dentro del plano de la sección.
- **Coste computacional:** al ser 2D, los códigos son una fracción del coste de algoritmos 3D completos, pero no capturan comportamientos fuera del plano (difracción 3D, spreading geométrico esférico, etc.).
- **Salida del código TM:** cubo de datos multi-offset (common-source gathers) de la componente $E_y$ en función del tiempo y posición del receptor, para cada posición de fuente. Los common-offset gathers se extraen trivialmente del cubo multi-offset.
- **Aplicación de tomografía crosshole:** los tiempos de primera llegada del dato sintético sin ruido se invierten con tomografía de rayos rectos (least-squares con regularización de segunda derivada) para recuperar el modelo de velocidades. El resultado es un tomograma de velocidad que concuerda bien con el modelo verdadero, con algo de suavizado en direcciones de rayos de alto ángulo.
- **CPML (Convolutional PML):** la implementación CPML de Roden y Gedney (2000) es independiente del medio (funciona igual para cualquier $\varepsilon$, $\mu$, $\sigma$), evita la división de componentes de campo de otras implementaciones PML, y absorbe tanto ondas propagantes como evanescentes. En los códigos se usa $m = 4$ y $\kappa_{k_{\max}} = 5$ como valores por defecto, y el espesor PML equivale al número de celdas PML especificado por el usuario.
- **Limitación principal:** el código no modela dispersión en propiedades eléctricas, por lo que no es directamente aplicable a suelos con alta dispersión (e.g., arcillas húmedas en frecuencias bajas de GPR).
