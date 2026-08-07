# Propuestas de artículos para completar el "Estado da Arte"

Documento de trabajo — **no** es parte de la tesis. Reúne artículos reales, encontrados por búsqueda web, para llenar las cuatro secciones `\section{Introdução e Estado da Arte}` que hoy están vacías (`[Estado da arte a ser desenvolvido.]`) en `dissertacao_tomo_sar.tex`: Capítulos 2, 3, 4 y 5. Ver `revision.md`, hallazgo 2.1.

**Importante antes de citar:** solo verifiqué título, autores y venue por búsqueda web. Confirmá volumen/número/páginas/DOI exactos en el enlace de cada entrada antes de agregarla al `\begin{enumerate}` de `Referências` (formato IEEE que ya usa la tesis).

Ya hay ~20 referencias en la lista actual de la tesis; el TODO de la propia lista pide llegar a 40–50 (55% artículos de revista, 25% libros, 20% conferencias). Lo que sigue cubre sobre todo artículos de revista/conferencia recientes — todavía falta buscar libros de texto adicionales para completar esa proporción.

---

## Capítulo 2 — Fundamentos Teóricos / motivación general (§1.1, §2.1)

### Griffiths (2003) — bistatic radar, principios y práctica
- **Cita:** H. D. Griffiths, "From a different perspective: principles, practice and potential of bistatic radar," *Proc. RADAR 2003*, Adelaide, Australia, 3–5 Sep. 2003.
- **Relevancia:** es el "Griffiths 2003" que el propio TODO original de §5.1 (Cap. 5) ya pedía citar como referencia de fundamentos de radar biestático. Confirmado que existe con ese título/venue.
- **Dónde citarlo:** §2.1 (Radar de Abertura Sintética: Princípios e Histórico) o §1.1 (Contexto e Motivação), junto a Willis (1995) y Horne & Yates (2002), ya presentes en la lista.

### Bistatic Landmine and IED Detection Combining Vehicle and Drone Mounted GPR Sensors
- **Cita:** *Remote Sensing* (MDPI), vol. 11, n. 19, 2299, 2019.
- **URL:** https://www.mdpi.com/2072-4292/11/19/2299
- **Relevancia:** geometría biestática Tx-en-vehículo / Rx-en-dron para sensoriamento subsuperficial — motivación directa y concreta para la configuración biestática de la tesis.
- **Dónde citarlo:** §1.1 (Contexto e Motivação), junto a la mención de "minas antipessoais, tubulações, camadas estratigráficas".

### A Lightweight and Low-Power UAV-Borne Ground Penetrating Radar Design for Landmine Detection
- **Cita:** *Sensors* (MDPI), vol. 20, n. 8, 2234, 2020.
- **URL:** https://www.mdpi.com/1424-8220/20/8/2234
- **Relevancia:** refuerza la motivación de plataformas VANT para GPR subsuperficial — un antecedente aplicado directo del "sistema de sensoriamento proposto" (§1.2).
- **Dónde citarlo:** §1.1 o §1.2 (O Sistema de Sensoriamento Proposto).

---

## Capítulo 3 — Resolução Espacial 3D (tomografía SAR circular/helicoidal, k-espacio)

### A review of SAR tomography (2025)
- **Cita:** *Geo-spatial Information Science* (Taylor & Francis), 2025. DOI: 10.1080/10095020.2025.2510365.
- **URL:** https://www.tandfonline.com/doi/full/10.1080/10095020.2025.2510365
- **Relevancia:** revisión amplia y muy reciente de TomoSAR — buena ancla de "estado del arte actual" para abrir §3.1, y permite mostrar que el campo sigue activo (menciona ALOS-2/4, SAOCOM, NISAR, ROSE-L, Tandem-L, BIOMASS).
- **Dónde citarlo:** primera(s) línea(s) de §3.1 (Introdução e Estado da Arte).

### The Status of Technologies to Measure Forest Biomass and Structural Properties: State of the Art in SAR Tomography of Tropical Forests
- **Cita:** *Surveys in Geophysics*, 2019. DOI: 10.1007/s10712-019-09539-7.
- **URL:** https://link.springer.com/article/10.1007/s10712-019-09539-7
- **Relevancia:** review de TomoSAR con foco en penetración de la señal en medios volumétricos (bosque) — análogo conceptual útil al caso subsuperficial (solo) de la tesis.
- **Dónde citarlo:** §3.1, junto a Reigber & Moreira (2000) y Banda et al. (2016), ya presentes en la lista.

### Corrección a la cita "Nolan & Fahnestock (2004)" del TODO original
El comentario TODO de §3.1 pide citar "Nolan & Fahnestock (2004)", pero no encontré ese trabajo. Lo que sí existe y es ampliamente citado en tomografía SAR e inversión de trayectorias arbitrarias es:
- **C. J. Nolan y M. Cheney**, "Synthetic aperture inversion for arbitrary flight paths and non-flat topography," *IEEE Trans. Geosci. Remote Sens.*, vol. 41, n. 5, 2003.
- **C. J. Nolan y M. Cheney**, "Microlocal analysis of synthetic aperture radar imaging," *J. Fourier Anal. Appl.*, vol. 10, n. 2, 2004.

Revisá si el original quería decir "Nolan & Cheney" (posible mezcla de nombres con otro autor de la lista de referencias de vuelo, como Fahnestock en literatura de glaciología) antes de decidir cuál de las dos citar, o si citar ambas.

---

## Capítulo 4 — Transmissão em Meios Estratificados (matrices de transferencia, medios estratificados)

### Subsurface characterization of stratified media via co-polarized GPR refracted wave inversion (2026)
- **Cita:** *ScienceDirect* (revista a confirmar en el enlace), 2026. DOI: 10.1016/j.jappgeo.2026.S0926985126001047 (verificar).
- **URL:** https://www.sciencedirect.com/science/article/abs/pii/S0926985126001047
- **Relevancia:** **el más cercano a tu propio enfoque** — onda refractada en medio estratificado, inversión co-polarizada. Es casi obligatorio citarlo y diferenciar explícitamente tu contribución (formalismo ABCD con última capa semi-infinita e impedâncias assimétricas) frente a este trabajo.
- **Dónde citarlo:** §4.1 (Introdução e Estado da Arte), como referencia más cercana/reciente, con una frase que marque la diferencia de enfoque.

### Fast and Rigorous Modeling of Antenna–Medium Interactions Above Planar Stratified Media via the Generalized Scattering Matrix (2025)
- **Cita:** arXiv:2504.12613, 2025.
- **URL:** https://arxiv.org/pdf/2504.12613
- **Relevancia:** formalismo alternativo (matriz de scattering generalizada, no ABCD) para el mismo tipo de problema — buen contraste metodológico para justificar la elección de ABCD en la tesis.
- **Dónde citarlo:** §4.1, junto a Pozar (2011).

### Application of ground penetrating radar methods in soil studies: A review
- **Cita:** *Geoderma* (ScienceDirect), 2018.
- **URL:** https://www.sciencedirect.com/science/article/abs/pii/S0016706118303823
- **Relevancia:** review general de aplicaciones de GPR en suelos estratificados — buen contexto aplicado para complementar Daniels (2004), ya en la lista.
- **Dónde citarlo:** §4.1.

### Ultra-broad-band electrical spectroscopy of soils and sediments — a combined permittivity and conductivity model
- **Cita:** *Geophysical Journal International*, vol. 210, n. 3, pp. 1360–1375, 2017.
- **URL:** https://academic.oup.com/gji/article/210/3/1360/3861100
- **Relevancia:** doble uso — sostiene el modelo Cole-Cole que la tesis ya propone como *Perspectiva Futura* (§7.3, "extensão a meios dissipativos") y da contexto cuantitativo a la limitación de "meios sem perdas" (§7.2).
- **Dónde citarlo:** §4.1 y opcionalmente como referencia de apoyo en §7.3.

---

## Capítulo 5 — Otimização Multi-Objetivo da Trajetória (optimización de trayectoria en SAR biestático)

### Drozdowicz & Samczyński — Drone-Based 3D Synthetic Aperture Radar Imaging with Trajectory Optimization (2022)
- **Cita:** J. Drozdowicz y P. Samczyński, "Drone-Based 3D Synthetic Aperture Radar Imaging with Trajectory Optimization," *Sensors* (MDPI), vol. 22, n. 18, 6990, 2022. DOI: 10.3390/s22186990.
- **URL:** https://www.mdpi.com/1424-8220/22/18/6990
- **Relevancia:** **el trabajo más directamente comparable** al Capítulo 5 — dron, optimización de trayectoria, SAR 3D. Ellos optimizan PSLR/ISLR/tiempo de vuelo, **no** ángulo de Brewster ni resolución vía k-espacio biestático — excelente para enunciar la brecha explícitamente.
- **Dónde citarlo:** primera línea de §5.1 (Introdução e Estado da Arte), como referencia principal de brecha.

### Path Planning for GEO-UAV Bistatic SAR Using Constrained Adaptive Multiobjective Differential Evolution (2016)
- **Cita:** *IEEE Trans. Geosci. Remote Sens.*, 2016.
- **URL:** https://ieeexplore.ieee.org/document/7516590/
- **Relevancia:** optimización multiobjetivo explícita de trayectoria en un sistema SAR biestático — buen precedente metodológico para justificar la formulación logarítmica de $J(\vect{r}_{Rx};\alpha)$ (§5.2–5.3).
- **Dónde citarlo:** §5.1, junto a García-Fernández et al. (2019), ya en la lista.

### Joint 3D Trajectory and Power Allocation for HAPs–UAV Bistatic ISARAC in Low-Altitude Networks (2026)
- **Cita:** arXiv:2606.04600, 2026.
- **URL:** https://arxiv.org/html/2606.04600
- **Relevancia:** muy reciente — trayectoria 3D conjunta + biestático + integración sensado-comunicación. Útil para mostrar que el problema de trayectoria biestática sigue siendo un área de investigación activa en 2026.
- **Dónde citarlo:** §5.1, como referencia de actualidad.

### Patente — uso del ángulo de Brewster en GPR de banda milimétrica
- **Cita:** "Method and apparatus for using collimated and linearly polarized millimeter wave beams at Brewster's angle of incidence in ground penetrating radar to detect objects located in the ground," patente US.
- **URL:** https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/7893862
- **Relevancia:** antecedente directo del uso del ángulo de Brewster en GPR, pero **para maximizar acoplamiento en una geometría fija**, no como criterio de optimización de una trayectoria continua del receptor. Es la pieza que mejor sostiene la afirmación de brecha ya presente en la Introdução ("nenhum trabalho utilizou o ângulo de Brewster como guia geométrico para otimização de trajetória do receptor biestático").
- **Dónde citarlo:** §5.1 y/o §1.1, explícitamente contrastado con la contribución de la tesis.

---

## Resumen y siguientes pasos

| Capítulo | Artículos propuestos | Cubre el TODO original |
|---|---|---|
| Cap. 2 / Introdução | 3 | Parcial (motivación biestática/subsuperficial) |
| Cap. 3 (Resolução) | 2 + corrección de 1 cita | Parcial (falta más SAR circular/helicoidal puro) |
| Cap. 4 (Multicamadas) | 4 | Bueno (ya tenía Pozar/Born&Wolf/Daniels/Ishimaru/Banda) |
| Cap. 5 (Otimização) | 4 | Bueno (cubre la brecha central de la tesis) |

Con esto más las ~20 referencias ya presentes, la lista quedaría cerca de 30–33. Para llegar a las 40–50 que pide el TODO de `Referências` (con 25% libros y 20% conferencias), todavía falta una ronda adicional enfocada en **libros de texto** (electromagnetismo aplicado, procesamiento de señales de radar) y **actas de conferencia** (IGARSS, EuRAD, RADAR) — puedo hacer esa búsqueda si querés.
