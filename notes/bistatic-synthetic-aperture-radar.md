# Bistatic Synthetic Aperture Radar

**Autores:** A. M. Horne, G. Yates  
**Año:** ~2002–2003 (presentado en conferencia)  
**Fuente/Publicación:** QinetiQ Malvern, U.K. Proceedings paper (conferencia de radar, UK MoD Corporate Research Programme TG9). Páginas 6–10 del documento.

---

## 1. Geometría del Sistema

El sistema es un SAR **biestático**: el transmisor y el receptor están ubicados en **plataformas separadas**.

- El transmisor y el receptor se mueven independientemente, con vectores de posición y velocidades distintas.
- Se define un marco de coordenadas cartesiano centrado en la escena, con ejes $u$ y $v$ orientados **perpendicular a iso-Dopplers** e **iso-rangos** respectivamente.
- Los ejes $u$ y $v$ no son en general ortogonales (excepto en el caso monoestático).
- Se definen los ángulos:
  - $\alpha_0$ — orientación del transmisor respecto a la escena en el centro de apertura
  - $\beta_0$ — orientación del receptor respecto a la escena en el centro de apertura
  - $\varphi$ — ángulo biestático: ángulo entre las líneas focales del transmisor y del receptor hacia la escena
- Las **iso-rangos biestáticas** son **elipses** con focos en el transmisor y el receptor (en lugar de círculos en el caso monoestático).
- Los **iso-Dopplers biestáticos** son hipérbolas.
- **Modo spotlight biestático:** la huella de la antena única ("single antenna footprint") se ilumina sin restricciones en el tiempo de apertura. Se concentra en este modo porque elimina la necesidad de técnicas de "pulse chasing" y permite control independiente de los tiempos de apertura.

---

## 2. Ecuaciones de Resolución SAR

### Orientación de los ejes de imagen biestáticos

$$\tan\theta_u = \frac{\Delta\alpha \cos\alpha_0 + \Delta\beta \cos\beta_0}{\Delta\alpha \sin\alpha_0 + \Delta\beta \sin\beta_0}, \qquad \theta_v = \frac{\alpha_0 + \beta_0}{2}$$

**Variables:**
- `\Delta\alpha` — tasa angular del transmisor respecto a la escena
- `\Delta\beta` — tasa angular del receptor respecto a la escena
- `\alpha_0, \beta_0` — orientaciones del transmisor y receptor respecto a la escena en el centro de apertura
- `\theta_u` — orientación del eje de resolución $u$ (perpendicular a iso-Dopplers)
- `\theta_v` — orientación del eje de resolución $v$ (perpendicular a iso-rangos biestáticos)

---

### Resoluciones espaciales biestáticas en los ejes $u$ y $v$

$$\rho_u = \frac{c}{f_0 K}, \qquad \rho_v = \frac{c}{2\Delta f \cos(\varphi/2)}$$

donde:

$$K^2 = (\Delta\alpha)^2 + (\Delta\beta)^2 + 2\Delta\alpha\,\Delta\beta\cos\varphi$$

**Variables:**
- `c` — velocidad de la luz (m/s)
- `f_0` — frecuencia portadora del radar (Hz)
- `\Delta f` — ancho de banda del sistema (Hz)
- `\varphi` — ángulo biestático entre las líneas focales del transmisor y receptor hacia la escena (rad)
- `K` — factor efectivo de apertura angular biestática
- `\rho_u` — resolución en la dirección $u$ (acimut biestático efectivo)
- `\rho_v` — resolución en la dirección $v$ (rango biestático efectivo)

---

### Efecto del ruido de fase en la función de punto extendido (PSF) biestática

$$\langle |h(u)|^2 \rangle \approx |h_i(u)|^2 + \frac{\Psi(f)}{2} \int_f \left\{|h_i\!\left(u - \frac{R_0 c}{2f_0 v}\right)|^2 + |h_i\!\left(u + \frac{R_0 c}{2f_0 v}\right)|^2\right\} df$$

**Variables:**
- `h_i(u)` — función de punto ideal (sin ruido de fase)
- `\Psi(f)` — espectro de ruido de fase unilateral (single-sideband)
- `R_0` — rango representativo (m)
- `v` — velocidad representativa de la plataforma (m/s)
- Para el sistema biestático, el espectro de ruido de fase visto es la **suma** de las contribuciones de los dos osciladores independientes (3 dB mayor que en el caso monoestático)

---

## 3. Suposiciones del Modelo

1. Se trabaja en el **modo spotlight biestático** (BSSAR), donde la huella de antena única permanece iluminada durante toda la apertura, eliminando la necesidad de "pulse chasing".
2. Se utiliza una imagen monoestática de la reflectividad compleja de la escena como medida para la síntesis de datos biestáticos simulados.
3. El procesado biestático se basa en una **integral doble** de los datos RAW a lo largo de la apertura sintética y a través del ancho de banda del sistema: computacionalmente intensivo pero exacto.
4. Para procesado práctico, se adaptan algoritmos monoestáticos (PFA — Polar Format Algorithm) al caso biestático mediante dos modificaciones:
   - Compensación de la apertura de la plataforma de recepción (motion compensation extendida).
   - Reproyección de los datos compensados en el espacio-K a un ángulo que bisecta el ángulo biestático.
5. Los osciladores del transmisor y del receptor son **independientes**, lo que hace al sistema biestático susceptible al ruido de fase de ambos; el ruido de fase efectivo visto es la suma de las potencias de ruido de los dos osciladores.
6. La sincronización temporal entre TX y RX requiere precisión mejor que 100 ns para geolocalización precisa; la sincronización en frecuencia es más exigente y depende de los requisitos de resolución y de la duración de la apertura.
7. Los dos enfoques para la sincronización biestática son: (a) osciladores atómicos independientes de alta precisión en cada plataforma, o (b) transferencia continua de referencias de tiempo y frecuencia entre plataformas (enlace RF directo, enlace GPS/GPS).
8. Los efectos del ángulo biestático sobre la resolución representan una reducción respecto al sistema monoestático con el mismo ancho de banda y el mismo tiempo de apertura sintética.

---

## 4. Notas Adicionales

- El paper es una **revisión de investigación en curso** en QinetiQ Malvern (ex-RSRE/DERA), financiada por el programa de investigación corporativo del MoD del Reino Unido (TG9).
- La principal ventaja operativa del SAR biestático es que la plataforma receptora (cara y con procesador) puede permanecer **encubierta** (covert), mientras que el transmisor (más vulnerable) puede ser una plataforma más barata o incluso expendable.
- La resolución biestática es **menor** que la monoestática con el mismo ancho de banda y tiempo de apertura sintética. La función de punto extendido biestática puede ser **no ortogonal** (los ejes de resolución $u$ y $v$ no son necesariamente perpendiculares), lo que puede afectar la interpretabilidad de la imagen.
- Se demuestran imágenes procesadas a partir de **datos biestáticos simulados** usando el método de integral doble (verificación del modelo teórico) y usando el **algoritmo PFA biestático** (para eficiencia computacional).
- El ruido de fase biestático tiene implicaciones distintas al monoestático: el ruido de baja frecuencia no se cancela porque los osciladores del TX y RX son independientes, afectando especialmente a los sistemas con largas integraciones (apertura larga).
- Los resultados experimentales muestran imágenes de resolución sub-métrica (0.3 m) simuladas a 100 km de rango con velocidad de plataforma de 200 m/s.
- La investigación futura se enfoca en combinar todos los componentes para demostrar un experimento SAR biestático de alta resolución completamente sincronizado.
