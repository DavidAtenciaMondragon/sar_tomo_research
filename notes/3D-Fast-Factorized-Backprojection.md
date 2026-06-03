# 3D Fast Factorized Back-Projection in Cartesian Coordinates

**Autores:** Juliana A. Góes, Valquiria Castro, Leonardo Sant'Anna Bins, Hugo E. Hernandez-Figueroa  
**Año:** 2020  
**Fuente/Publicación:** 2020 IEEE Radar Conference (RadarConf20). UNICAMP / INPE, Brasil.

---

## 1. Geometría del Sistema

Sistema SAR 3D con **trayectoria helicoidal** (aplicación principal: dron a alturas entre 80–120 m, radio 180 m, 5 vueltas, velocidad 6.5 m/s). El algoritmo es general y aplicable a cualquier trayectoria. El volumen 3D se divide en subimágenes mediante una **curva de Morton modificada** (Z-order curve extendida a 3D), formando un árbol octree flexible con partición $D_{px} \times D_{py} \times D_{pz}$.

Parámetros de simulación: $\lambda_0 = 0.75$ m (P-band), $B = 150$ MHz, resolución en rango = 1 m, resolución en plano (x,y) = 0.16 m (−3 dB), resolución en z = 1.53 m.

---

## 2. Ecuaciones de Resolución SAR

### Índice 3D de Morton modificado

$$I_{3D} = q_x + D_{px} q_y + D_{px} D_{py} q_z$$

**Variables:**
- `q_x, q_y, q_z` — secuencias de índice en cada dimensión
- `D_{px}, D_{py}, D_{pz}` — número de divisiones por recursión en cada eje

### Centros de fase de sub-apertura (caso L impar)

$$\mathbf{r}_c[k] = \mathbf{r}_0\!\left[(L^c-1)/2 + kL^c\right]$$

**Variables:**
- `\mathbf{r}_0` — array de posiciones SAR originales (nodo raíz)
- `L` — número de sub-aperturas padre combinadas para formar una hija
- `k = 0, \ldots, N_c - 1` — índice del elemento hijo
- `N_c` — número de elementos de apertura en el nodo hijo

### Centros de fase (caso L potencia de 2)

$$\mathbf{r}_c[k] = \boldsymbol{\rho}_0\!\left[L^c/2 - 1 + kL^c\right]$$

**Variables:**
- `\boldsymbol{\rho}_0` — puntos medios entre posiciones SAR consecutivas

### Datos SAR del nodo hijo a partir del padre

$$s_c[k,m] = \sum_{l \in \mathcal{A}_{p,k}} s_p[l, v_{c,k,l,m}] \cdot \Delta\Phi_{c,k,l,m}$$

**Variables:**
- `s_p[l,m]` — datos radar del nodo padre, pulso $l$, muestra de rango $m$
- `v_{c,k,l,m}` — índice de punto flotante en $s_p$ (interpolación lineal)
- `\mathcal{A}_{p,k}` — índices de los pulsos padre que componen la sub-apertura $k$ del hijo

### Término de compensación de fase

$$\Delta\Phi_{c,k,l,m} = \exp\!\left\{-j\frac{4\pi}{\lambda_0}\!\left(\Psi_{c,k,m} - \Gamma_{c,k,l,m}\right)\right\}$$

**Variables:**
- `\lambda_0` — longitud de onda portadora
- `\Psi_{c,k,m}` — rango desde el centro de fase hijo $\mathbf{r}_c[k]$ al bin de rango $m$
- `\Gamma_{c,k,l,m}` — rango desde el centro de fase padre $\mathbf{r}_p[l]$ al mismo bin de rango $m$

### Resolución simulada (3 dB)

| Dirección | Resolución |
|-----------|-----------|
| Plano (x,y) | 0.16 m |
| Dirección z | 1.53 m |
| PSLR (x,y) | −9.1 dB |
| PSLR (z) | −28.7 dB |

### Reducción de tiempo de procesado

$$\text{Speedup} = \frac{t_{BP}}{t_{3D\text{-}FFBP}} = \frac{39.09\text{ h}}{3.45\text{ h}} \approx 11.3\times$$

---

## 3. Suposiciones del Modelo

1. Trayectoria de vuelo arbitraria (sin restricción a línea recta).
2. La imagen 3D tiene resolución no uniforme: mejor en el plano (x,y) que en z, inherente a la geometría helicoidal.
3. Los centros de fase de las sub-aperturas se calculan por interpolación solo en el nodo raíz; en los demás nodos se usan las expresiones recursivas (Ecs. 4–7).
4. El número de bins de rango es constante para cada nodo (cobertura esférica completa).
5. La suma coherente de sub-aperturas requiere compensación de fase mediante $\Delta\Phi_{c,k,l,m}$ para corregir la diferencia entre la trayectoria real y la trayectoria aproximada de la sub-apertura.
6. Error de posición máximo del sistema drone-SAR: 7.4 mm → desviación estándar de fase 0.12 rad (< umbral $\pi/8$ recomendado).

---

## 4. Notas Adicionales

- Extensión 3D del algoritmo FFBP de Ulander et al. (2003), usando coordenadas **Cartesianas** en lugar de polares, lo que lo hace apto para trayectorias arbitrarias.
- La curva de Morton 3D permite implementar el árbol octree con acceso eficiente a datos 3D.
- El grado de coherencia entre 3D-FFBP y BP directo es 0.9993 (casi idénticos).
- Aplicación directa a SAR tomografía con drones (referenciado en los papers del mismo grupo de UNICAMP/INPE).
- Código implementado en MATLAB R2018a, ejecutado en Intel i7-7700 (3.60 GHz), 64 GB RAM.
