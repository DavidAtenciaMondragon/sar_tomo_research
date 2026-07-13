# Guía de Estructura para una Tesis Doctoral en Ingeniería (SAR/Radar)

Análisis basado en dos tesis de doctorado de UNICAMP (Brasil):
- **Góes (2022)** — *Techniques for High-Resolution 3D Images with Synthetic Aperture Radar*, ~179 páginas, 64 referencias, 1 apéndice. Estilo IEEE.
- **Flores (2019)** — *Multipath Influence Analysis on Detection Probability of Radars and Optimal Estimators of Velocity in ArcSAR*, ~84 páginas, ~37 referencias, 8 apéndices. Estilo ABNT.

---

## 1. Partes obligatorias (ante-texto)

```
Portada / Cover page
Ficha catalográfica
Dedicatoria (opcional)
Agradecimientos / Acknowledgements
Resumen (idioma nativo + inglés, cada uno ~250 palabras)
Lista de figuras
Lista de tablas
Lista de acrónimos / abbreviations
Lista de símbolos (con unidades)
Tabla de contenidos (Summary)
```

### Observaciones
- El **resumen en inglés** es requisito institucional aunque la tesis esté en otro idioma.
- La **lista de símbolos** debe incluir: símbolo, descripción y unidades. Ambas tesis la tienen de 2–3 páginas.
- Las listas de figuras y tablas se numeran por capítulo (e.g., Figura 2.1, Tabla 4.1).
- La tabla de contenidos llega hasta nivel 3 de subsección (2.2.1.1 en Góes).

---

## 2. Esqueleto de capítulos

```
1. Introducción            (~5–10 páginas)
2. Marco teórico           (~20–30 páginas)
3. Contribución principal 1 (~30–40 páginas)
4. Contribución principal 2 (~30–40 páginas)
5. Conclusiones            (~5–8 páginas)
Bibliografía
Apéndices
```

Cada capítulo de contribución tiene la siguiente anatomía interna:

```
X.1  Introducción del capítulo (contexto + qué hace este capítulo)
X.2  Estado del arte / revisión de literatura
X.3  Modelo teórico / derivación
X.4  Resultados de simulación
     X.4.1  Escenario de simulación (tabla de parámetros)
     X.4.2  Caso 1 (e.g., superficie plana)
     X.4.3  Caso 2 (e.g., superficie real)
     X.4.4  Discusión
X.5  Validación con datos reales (si aplica)
     X.5.1  Descripción de la campaña / dataset
     X.5.2  Resultados procesados
     X.5.3  Comparación teórico vs. real (tabla)
     X.5.4  Discusión
X.6  Conclusiones del capítulo
```

> **Regla de oro:** simulación primero, datos reales después. Los datos reales validan lo que la simulación predijo.

---

## 3. Capítulo 1 — Introducción

### Contenido y orden recomendado

1. **Contexto amplio** (2–3 párrafos): panorama del área, importancia del problema, por qué es relevante hoy.
2. **Sistema o plataforma concreta** (si aplica): describe el hardware/dataset real que usarás. Muestra una foto o diagrama. Justifica por qué ese sistema.
3. **Planteamiento del problema**: qué gap existe en el estado del arte, qué limitación resuelve esta tesis.
4. **Contribuciones** (subsección explícita, con lista): enumera 3–5 aportaciones originales en viñetas. Sé específico y cuantitativo donde sea posible.
5. **Esquema de la tesis** (outline): un párrafo por capítulo explicando qué contiene y por qué viene en ese orden.
6. **Publicaciones derivadas**: lista de artículos publicados/enviados como resultado de esta tesis.

### Patrón de escritura
- El primer párrafo establece el campo con 1–2 referencias históricas fundacionales.
- Cada afirmación importante va respaldada por una cita `[N]` (IEEE) o `(AUTOR, AÑO)` (ABNT).
- La "brecha" de conocimiento se señala explícitamente con frases como *"However, existing models are only valid for..."* o *"To the best of our knowledge, there are no similar results in the literature."*

---

## 4. Capítulo 2 — Marco teórico / Fundamentals

### Propósito
Construir desde primeros principios el conocimiento que el lector necesita para entender los capítulos de contribución. No es un "survey" exhaustivo: es selectivo y pedagógico.

### Patrón de construcción teórica

```
Ley/fenómeno físico fundamental
  → Ecuación básica (numerada, derivada paso a paso)
  → Aproximaciones y límites de validez
  → Ejemplo numérico con números concretos
  → Figura ilustrativa
  → Conexión con el capítulo siguiente
```

**Ejemplo de Góes (Cap. 2):**
```
Radar equation (2.6) → azimuth resolution (2.2)–(2.3) → chirp pulse (2.14)
→ range compression (2.15)–(2.19) → Stripmap geometry → Back Projection algorithm
```
Cada concepto nuevo se ancla en el anterior. El lector nunca da un salto conceptual sin apoyo.

### Profundidad esperada
- 3 niveles de subsección máximo (2.2.1.1).
- Cada ecuación que aparezca en los capítulos de contribución debe estar definida aquí.
- Tablas de parámetros de referencia (e.g., bandas de frecuencia SAR, Tabla 2.1 en Góes).
- Figuras de geometría con coordenadas y variables etiquetadas.

---

## 5. Capítulos de Contribución (3, 4, ...)

### 5.1 Subsección de introducción del capítulo

Cada capítulo abre con:
1. **Qué problema resuelve** este capítulo (1 párrafo).
2. **Breve estado del arte** relevante para este capítulo en particular (puede extenderse a 1–2 páginas si no hay capítulo de revisión separado).
3. **Qué propone este capítulo** y cómo se organiza.

### 5.2 Desarrollo teórico / modelo propuesto

- Introduce el modelo con un **diagrama de escenario** primero (figura geométrica con variables).
- Define variables en la figura antes de las ecuaciones.
- Deriva las ecuaciones principales paso a paso, numerándolas todas.
- Señala las **premisas y limitaciones** del modelo explícitamente.
- Cita el trabajo previo que sirve de base con honestidad (*"Following the same principle of [REF]..."*).

### 5.3 Resultados de simulación

#### Formato estándar:
1. **Tabla de parámetros del escenario** (siempre): incluye todos los valores numéricos usados.
2. **Figura de resultado** con descripción completa en el caption.
3. **Párrafo de interpretación** inmediatamente después de cada figura: qué muestra, qué significa, cómo se compara con la predicción teórica.
4. **Subsección de discusión** al final: sintetiza los hallazgos, señala limitaciones, anticipa la validación con datos reales.

#### Tipos de figuras de simulación observadas:

| Tipo | Uso |
|------|-----|
| PSF (Point Spread Function) / Isosurface | Demostrar resolución del algoritmo |
| Curvas μ ± σ vs. SNR | Comparar estimadores estadísticos |
| Mapas 2D con colorbar | Distribución espacial de una métrica |
| Histogramas + curva teórica superpuesta | Validar distribución analítica con MC |
| Curvas comparativas (varias líneas) | Comparar métodos o escenarios |
| Tablas de comparación teórico vs. real | Confirmar predicciones |

#### Validación analítica vs. simulación (Flores):
- Se usan 10⁶ realizaciones de Monte Carlo.
- La curva analítica se superpone al histograma simulado.
- Si son indistinguibles, las premisas del modelo son válidas.

### 5.4 Resultados con datos reales

#### Estructura cuando existe campaña de vuelo/experimento:
1. **Descripción del experimento**: plataforma, parámetros del sistema, fecha/lugar, blanco usado (e.g., corner reflector).
2. **Pre-procesamiento**: pasos de igualización, filtrado, upsampling (menciona por qué son necesarios).
3. **Imágenes o resultados procesados**: figuras con ejes etiquetados y unidades.
4. **Tabla de comparación cuantitativa**: valores teóricos vs. valores medidos, con porcentaje de error o ratio.
5. **Discusión**: ¿confirman los datos la simulación? Si hay discrepancia, explica por qué. Señala limitaciones del experimento.

**Ejemplo de Góes (Tabla 4.10):**
```
| Resolución | Teórico (0,0) | Teórico (-2, 25.4) | Valor real |
|------------|--------------|---------------------|------------|
| δx [cm]    | 16.6         | 16.8                | 16.6       |
| δy [cm]    | 18.8         | 19.6                | 17.4       |
| δz [m]     | 1.98         | 1.52                | 2.33 / 1.76|
```

### 5.5 Subsección "Discussion" dentro del capítulo

Aparece al final de cada bloque de resultados (tanto simulación como datos reales):
- Sintetiza los hallazgos en 3–5 párrafos.
- Conecta con las hipótesis planteadas al inicio.
- Reconoce limitaciones con honestidad (*"...could not be confirmed nor denied. The reason is that..."*).
- Señala mejoras para trabajos futuros relacionadas con ese tema específico.

### 5.6 Conclusiones del capítulo (sección X.N al final)

Un párrafo corto (o lista de viñetas como en Flores) que resume:
- Qué se propuso.
- Qué se demostró.
- Cuál es la implicación práctica.

No repite lo mismo que la sección de discusión: es más sintético y orientado al impacto.

---

## 6. Capítulo de Conclusiones

### Estructura (Góes: prosa densa / Flores: viñetas por capítulo + subsección de trabajo futuro)

1. **Párrafo de apertura**: qué ofrece esta tesis en términos generales.
2. **Resumen por contribución**: un párrafo o viñeta por cada capítulo de contribución, con métricas cuantitativas específicas.
3. **Síntesis del impacto**: por qué estas contribuciones son importantes, aplicaciones habilitadas.
4. **Trabajo futuro** (puede ser subsección separada como en Flores §5.1): lista de extensiones naturales, experimentos pendientes, hipótesis no verificadas.

### Tono
- Afirmativo para lo demostrado: *"This thesis demonstrates that..."*, *"The results confirm..."*
- Condicional para lo pendiente: *"The hypothesis could neither be confirmed nor denied..."*, *"Future work should repeat this analysis for other frequency bands."*
- Las conclusiones citan números concretos: *"The speed-up factor is up to 21× for 3D images"*, *"Vertical resolution of 1.76 m"*.

---

## 7. Bibliografía

### Estadísticas observadas

| | Góes (2022) | Flores (2019) |
|--|--|--|
| Total de referencias | 64 | ~37 |
| Artículos de revista | ~35 (55%) | ~20 (54%) |
| Artículos de congreso | ~20 (31%) | ~10 (27%) |
| Libros / handbooks | ~5 (8%) | ~5 (14%) |
| Reportes / online | ~4 (6%) | ~2 (5%) |

### Estilos de citación
- **IEEE (Góes):** `[N]` en el texto, ordenadas por orden de aparición. Formato: `[N] Autores, "Título," Revista, vol., no., pp., mes año, DOI.`
- **ABNT (Flores):** `(APELLIDO, AÑO)` en el texto, ordenadas alfabéticamente. Incluye anotación *"Citado N vezes nas páginas X, Y, Z."* — útil para auditar la tesis.

### Revistas más citadas (área SAR/Radar, UNICAMP)
- IEEE Trans. on Geoscience and Remote Sensing (TGRS)
- IEEE Trans. on Aerospace and Electronic Systems (TAES)
- IEEE Geoscience and Remote Sensing Letters (GRSL)
- Remote Sensing (MDPI)
- Conferences: IGARSS, RadarConf, SBrT (Brasil)

---

## 8. Apéndices

### Criterio para mover algo a un apéndice
- Derivaciones matemáticas largas que interrumpirían el flujo del capítulo (Flores Apéndice D: angular error, Apéndice E: distribución del factor multipath).
- Implementaciones de código con anotaciones (Góes Apéndice A: ~12 páginas de MATLAB comentado con las ecuaciones correspondientes).
- Pseudocódigos de algoritmos (Flores Apéndice C).
- Lista de publicaciones generadas (Flores Apéndice A).
- Permisos de reutilización de figuras publicadas (Flores Apéndice B).

### Organización
- Los apéndices se numeran con letras: A, B, C, ...
- Se referencian desde el texto del capítulo: *"see Appendix D"*, *"as described in Appendix C"*.
- El MATLAB code appendix en Góes es particularmente valioso: muestra cómo mapear código a ecuaciones, línea por línea.

---

## 9. Convenciones de formato observadas

### Ecuaciones
- Numeradas `(X.Y)` donde X es el capítulo.
- Siempre referenciadas en el texto: *"as shown in Eq. (2.4)"*, *"from (3.11)"*.
- Las derivaciones intermedias también se numeran si se referencian después.
- Las variables se definen inmediatamente después de la ecuación en la que aparecen por primera vez.

### Figuras
- Numeradas `Figure X.Y` (Góes) o `Fig. X.Y` (Flores).
- Caption debajo de la figura, completo y autoexplicativo (el lector no debería necesitar leer el texto para entender qué muestra).
- Fuente incluida si es de otro trabajo: *"Source: Modified from [N]"* o *"Source: Author"*.
- Subfiguras etiquetadas (a), (b), (c) con descripción individual en el caption.
- Ejes siempre etiquetados con magnitud y unidades.

### Tablas
- Numeradas `Table X.Y`.
- Caption arriba de la tabla (a diferencia de las figuras).
- Siempre referenciadas en el texto antes de aparecer.

### Párrafos y prosa
- Párrafos de 4–6 líneas, alineación justificada.
- Cada párrafo tiene una idea central (no se mezclan ideas).
- No se usa lista de viñetas dentro del texto corriente (salvo en las secciones de contribuciones, conclusiones y trabajo futuro).
- Transiciones explícitas entre subsecciones: *"After discussing X in Section 2.2, the following section covers..."*

---

## 10. Concatenación de conocimientos (cómo se construye el argumento)

Ambas tesis siguen el mismo patrón lógico de tres niveles:

```
Nivel 1 – FUNDAMENTO
   "El fenómeno X existe y se rige por la ecuación Y."
   (Capítulo 2: marco teórico, citas a textbooks)

Nivel 2 – BRECHA
   "Sin embargo, el modelo existente Z no considera el caso W,
   lo que causa el problema P en la práctica."
   (Inicio de cada capítulo de contribución, citas a literatura reciente)

Nivel 3 – CONTRIBUCIÓN
   "Proponemos el método M que resuelve P.
   La simulación demuestra X, y los datos reales confirman Y."
   (Cuerpo del capítulo, resultados propios)
```

Esto se repite para cada contribución. Cada capítulo parte donde terminó el anterior:
- Góes Cap. 3 usa la formulación de BP del Cap. 2 para proponer FFBP.
- Góes Cap. 4 usa el FFBP del Cap. 3 para analizar geometrías espirales.
- Flores Cap. 3 usa el algoritmo de reflexiones del Cap. 2 para calcular probabilidad de detección.
- Flores Cap. 4 usa el modelo de señal del Cap. 3 para diseñar estimadores de velocidad.

> **Nunca hay capítulos "flotantes"**: cada uno justifica su existencia como paso necesario para el siguiente.

---

## 11. Publicaciones derivadas (aspecto importante no siempre documentado)

Ambas tesis listan las publicaciones como resultado directo de su trabajo. Góes tiene 3 publicaciones (2 artículos de revista + 1 congreso + 1 código en Zenodo). Flores tiene 5 publicaciones (3 artículos de revista/congreso + 2 en proceso de revisión).

**Patrón:** cada capítulo de contribución major normalmente generó al menos 1 publicación en congreso o revista. La tesis es esencialmente una versión extendida y unificada de esos artículos, con más contexto, derivaciones completas y validación adicional.

---

## 12. Aspectos que no son siempre obvios al escribir

### Equilibrio entre teoría y resultados
- Cap. 2 (teoría pura): 0% resultados propios.
- Cap. 3–4 (contribuciones): ~40% teoría/modelo, ~40% simulación, ~20% datos reales.
- Las tesis con solo simulación (Flores) son válidas si la derivación analítica es la contribución principal.

### Cómo reportar cuando la realidad difiere del modelo
Ambas tesis son honestas: cuando los datos reales no coinciden exactamente con la teoría, se explica la causa (igualización imperfecta, efectos de near-range, frecuencias no analizadas) y se lista como trabajo futuro. No se fuerza la interpretación.

### El rol de las limitaciones
Declarar las limitaciones del modelo no debilita la tesis: las fortalece. Demuestra que el autor entiende el alcance real de su contribución. En Góes: *"The hypothesis that decorrelation does not affect Spiral SAR could neither be confirmed nor denied."* En Flores: *"Validations of the results through field campaigns are required."*

### Longitud apropiada
- Tesis enfocada en ingeniería de señales/algoritmos: 80–180 páginas de contenido.
- No hay beneficio en extender artificialmente: la densidad de contenido importa más que el volumen.
- Los apéndices no se cuentan en el límite de páginas pero deben ser útiles (no relleno).

---

## 13. Lista de verificación antes de entregar

- [ ] Cada ecuación tiene número y cada variable definida en la primera aparición.
- [ ] Cada figura tiene caption completo, ejes etiquetados y fuente.
- [ ] Cada tabla tiene caption arriba, referencia en texto y fuente.
- [ ] El capítulo de contribución tiene: introducción, modelo, simulación, (datos reales), discusión, conclusiones.
- [ ] La sección de discusión señala limitaciones con honestidad.
- [ ] Las conclusiones contienen métricas numéricas concretas.
- [ ] La lista de símbolos cubre todo símbolo que aparezca en las ecuaciones.
- [ ] Las publicaciones derivadas están listadas (en apéndice o en la introducción).
- [ ] La bibliografía está en el estilo correcto (IEEE o ABNT) y todas las referencias están citadas en el texto.
- [ ] Los apéndices contienen solo material que interrumpiría el flujo del capítulo.
- [ ] El resumen (abstract) en inglés refleja las contribuciones y los resultados cuantitativos principales.
