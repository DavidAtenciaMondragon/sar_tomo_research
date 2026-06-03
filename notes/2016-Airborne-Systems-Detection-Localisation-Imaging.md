# Airborne Systems and Methods for the Detection, Localisation and Imaging of Buried Objects and the Characterisation of the Composition of the Subsurface

**Autores:** Borja Gonzáles Valdés, Yuri Álvarez López, Ana Arboleya, Yolanda Rodríguez Vaqueiro, María García Fernández, Fernando Las-Heras Andrés, Antonio García Pino  
**Año:** 2017 (publicación PCT, prioridad 21/01/2016)  
**Fuente/Publicación:** Patente Internacional WO 2017/125627 A1 — PCT/ES2017/000006. Solicitantes: Universidad de Oviedo (ES) y Universidad de Vigo (ES)

---

## 1. Geometría del Sistema

La invención describe un **sistema aerotransportado para detección, localización y obtención de imágenes de objetos enterrados** (minas antipersona, tuberías, restos arqueológicos) mediante radar de apertura sintética (SAR) y GPR (Ground Penetrating Radar) embarcado en UAV (vehículo aéreo no tripulado, preferentemente multirrotor/octacóptero).

### Sistema Monoestático

Un único módulo aéreo (UAV) que emite y recibe la señal radar:

- **Módulo aéreo (1):** contiene unidad radar (11), sistema de posicionamiento y guiado (13), unidad de control aérea (15).
  - Unidad radar (11): antena transmisora (111) + antena receptora (112) + módulo radar (113).
  - Módulo radar: PulsOn P410, banda 3–5 GHz (UWB).
  - Antenas helicoidales con polarización circular ortogonal entre Tx y Rx; ganancia >10 dB, S11 < −15 dB en 3–5 GHz.
- **Estación terrena (2):** control de vuelo (21), procesado de señal radar (23), aplicación de visualización (24), base RTK (22).
- **Comunicaciones (3):** IEEE 802.11 (Wi-Fi), Bluetooth, 3G/4G, ZigBee, o transceptores 433 MHz.

### Sistema Multiestático

Dos módulos aéreos (UAVs): un módulo emisor (101) y un módulo receptor (102), separados espacialmente. La sincronización se realiza mediante el sistema de comunicación radar (120) integrado en el módulo PulsOn P410.

### Geometría de Vuelo y Precisión

- El módulo aéreo vuela sobre la zona bajo estudio describiendo trayectorias de exploración (raster, zigzag, etc.) programadas en el sistema de control automático.
- El sistema de posicionamiento combina: GPS global (131), IMU/inercia (132), RTK (133) y fotogrametría (134).
- **Precisión de posicionamiento: ≤ 3 cm** en 3D (requisito crítico para procesado SAR coherente).
- Para el procesado SAR, las posiciones de adquisición deben estar **separadas ≤ λ/2** a la frecuencia más alta de trabajo (a 5 GHz: λ/2 = 3 cm).
- Altitud de vuelo típica: 50 cm – varios metros sobre la superficie.
- El radar ilumina el suelo con incidencia **perpendicular (downward-looking)**, permitiendo máxima energía transmitida al subsuelo.

---

## 2. Ecuaciones de Resolución SAR

### 2.1 Resolución en Rango / Profundidad

$$\delta_r = \frac{c}{2W}$$

**Variables:**
- `c` — velocidad de la luz ($3 \times 10^8$ m/s)
- `W` — ancho de banda del pulso radar (hasta 5 GHz de BW → $\delta_r = 3$ cm)

> A máximo ancho de banda (5 GHz), la resolución en profundidad (en aire) es ~3 cm.

### 2.2 Condición de Muestreo Espacial SAR

Para procesado SAR coherente, el espaciado máximo entre posiciones consecutivas del módulo aéreo debe ser:

$$\Delta x \leq \frac{\lambda_{\min}}{2} = \frac{c}{2 f_{\max}}$$

**Variables:**
- `\lambda_{\min}` — longitud de onda a la frecuencia máxima de trabajo
- `f_{\max}` — frecuencia máxima (5 GHz → $\Delta x \leq 3$ cm)

### 2.3 Imagen SAR — Suma Coherente (Back Projection)

$$s(\mathbf{h}) = \sum_{k} S(k, R_k(\mathbf{h})) \cdot e^{+j\frac{4\pi}{\lambda} R_k(\mathbf{h})}$$

donde $R_k(\mathbf{h}) = \|\mathbf{r}_k - \mathbf{h}\|$ es la distancia entre la posición del UAV en el pulso $k$ y el punto $\mathbf{h}$ del subsuelo.

**Variables:**
- `S(k, R)` — señal radar comprimida en rango para el pulso $k$, evaluada en la distancia $R$
- `\mathbf{r}_k` — posición 3D del módulo aéreo en el pulso $k$
- `\mathbf{h}` — punto 3D del subsuelo a reconstruir

### 2.4 Estimación de Permitividad — Método de Distancias

A partir de la imagen SAR, la permitividad relativa del subsuelo se estima mediante la diferencia de distancias entre el eco en el suelo y el eco en un objeto metálico de calibración enterrado:

$$\varepsilon_r = \left(\frac{c \cdot \Delta t_{\mathrm{air}}}{2 \cdot \Delta z_{\mathrm{real}}}\right)^2$$

equivalentemente, comparando la distancia aparente (imagen SAR con velocidad $c$) entre la superficie y el calibrador versus la distancia real conocida:

$$\varepsilon_r = \left(\frac{d_{\mathrm{aparente}}}{d_{\mathrm{real}}}\right)^2$$

**Variables:**
- `d_{\mathrm{aparente}}` — distancia entre eco del suelo y eco del calibrador en la imagen SAR (procesada con velocidad $c$)
- `d_{\mathrm{real}}` — profundidad real del objeto de calibración
- `\varepsilon_r` — permitividad relativa del subsuelo

### 2.5 Estimación de Permitividad — Método de Amplitudes

Alternativamente, la permitividad se estima a partir de la diferencia de amplitud entre el eco en el suelo y el eco del objeto metálico de calibración, usando el coeficiente de reflexión de Fresnel para incidencia normal:

$$\Gamma = \frac{\sqrt{\varepsilon_r} - 1}{\sqrt{\varepsilon_r} + 1}$$

**Variables:**
- `\Gamma` — coeficiente de reflexión en la interfaz aire-subsuelo (incidencia normal)
- La amplitud del eco del suelo es proporcional a $|\Gamma|$; la del calibrador metálico es máxima (reflexión total).

### 2.6 Velocidad de Propagación en el Subsuelo

Una vez estimada $\varepsilon_r$:

$$v = \frac{c}{\sqrt{\varepsilon_r}}$$

Esta velocidad se usa para reconstruir correctamente la posición de los objetos enterrados en el procesado SAR.

### 2.7 Eliminación de Clutter — Proceso Iterativo

El algoritmo de eliminación de clutter es iterativo:

1. Aplicar SAR a la matriz de señales radar → imagen de reflectividad 3D.
2. Identificar la región del suelo (reflexión en $z = 0$) y generar máscara.
3. Calcular la señal radar que produce la región enmascarada.
4. Restar dicha señal de la matriz original: $S_{\mathrm{new}} = S_{\mathrm{orig}} - S_{\mathrm{clutter}}$.
5. Aplicar SAR a $S_{\mathrm{new}}$ → imagen mejorada del subsuelo.
6. Repetir hasta convergencia.

---

## 3. Suposiciones del Modelo

1. El sistema de posicionamiento RTK proporciona precisión ≤ 3 cm, suficiente para procesado SAR coherente a frecuencias de hasta 5 GHz.
2. El módulo aéreo es un UAV multirrotor (octacóptero) con carga útil ≤ 1.5 kg.
3. La frecuencia máxima de trabajo es 5 GHz (banda de 3–5 GHz del PulsOn P410), equilibrando resolución en profundidad y penetración en subsuelo.
4. La superficie del suelo es aproximadamente **plana** localmente en el área bajo estudio.
5. El subsuelo tiene permitividad relativa **homogénea** (estimada mediante objeto de calibración), o al menos localmente homogénea.
6. Las antenas helicoidales con polarizaciones circulares ortogonales (Tx izquierda, Rx derecha) maximizan la energía reflejada por objetos de geometría compleja y reducen el clutter aire-suelo.
7. La oscilación del UAV (causada por viento o perturbaciones) introduce pequeños errores de fase que pueden corregirse con **Phase Gradient Autofocus (PGA)**.

---

## 4. Notas Adicionales

- La invención resuelve el problema principal de sistemas GPR-UAV anteriores: la **precisión de posicionamiento insuficiente** (GPS convencional: ~1 m; RTK: ≤ 3 cm) que impedía el procesado SAR coherente.
- A 5 GHz, $\lambda/2 = 3$ cm, que coincide exactamente con la precisión RTK del sistema → el sistema opera en el límite de coherencia.
- El **sistema multiestático** (dos UAVs) aumenta la diversidad espacial de iluminación, mejorando la detección de objetos difícilmente detectables con un solo ángulo de iluminación.
- Las antenas helicoidales con polarización circular (LHCP Tx, RHCP Rx) tienen la ventaja de que la reflexión especular en el suelo invierte la polarización, lo que **reduce el clutter superficial** respecto a la señal de objetos enterrados.
- El sistema permite escanear **25 m² en 10 segundos** (comparable a otros sistemas UAV), frente a los sistemas terrestres robotizados (~0.5 m en 10 s).
- Aplicaciones principales: detección de **minas antipersona** (59–69 millones enterradas globalmente), inspección de tuberías, arqueología, detección de cavidades.
- El **algoritmo de detección de objetos enterrados** (45) analiza la imagen 3D SAR en busca de agrupaciones de alta reflectividad, comparándolas con una base de datos de patrones de entrenamiento.
- El método de caracterización del subsuelo es **independiente del sistema radar** específico (no requiere recalibración si se cambia el hardware).
- Implementación práctica (Ejemplo 1): octacóptero 6 kg MTOW, 1.5 kg payload, módulo radar PulsOn P410 (3–5 GHz), antenas helicoidales 12 dB, RTK GPS (dos unidades), Wi-Fi para comunicaciones, procesado en laptop con MATLAB.
