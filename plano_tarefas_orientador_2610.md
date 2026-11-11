# Plan de trabajo: correcciones del orientador (Plano de Conclusão, 26/10/2026)

> **Documento de control.** Se avanza **una fase a la vez**. Cada fase tiene una
> *puerta de salida*. No se abre la siguiente fase hasta cerrar la puerta (o
> hasta registrar en §8 que el orientador aprobó saltarla).
> Regla del orientador: como máximo **dos frentes a la vez**, uno
> computacional/experimental y uno de escritura.
>
> Fuente: `Plano_Conclusao_Doutorado_David_2026-2027.docx.pdf` (17 pp.).
> Tesis revisada: `doc/dissertacao_tomo_sar.tex` (24/10/2026; 111 pp., 12 caps.).
> Fecha de este plan: **10/11/2026**.

---

## 0. Estado de partida

### 0.1 ¿La tesis tiene todos los estudios del repositorio?

| Estudio en el repositorio | ¿Está en la tesis? | Dónde / observación |
|---|---|---|
| Resolución en forma cerrada con refracción (`explicacion_resolucion`, `bisar/resolucion.py`) | Sí | Cap. 3 y Apéndice A |
| ABCD multicapa (`common/`, `plan_de_vuelo_multislab`) | Sí | Cap. 4 |
| Criterio de potencia / Brewster (`otimiza_tx`, `verificacao_numerica`, `relatorio_2809`) | Sí, **sobredimensionado** | Cap. 5 (~290 líneas) + parte de los caps. 2 y 4 |
| Motor `analise_1610` (Exp. A, B, D2) | Sí | Cap. 7 |
| Exp. D, D3, E, E2 (blancos reales, polarización) | Sí | Cap. 8 |
| Exp. C y C2 (superficie, nulo VV) | Sí | Cap. 9 |
| Exp. F y F2 (trayectoria, calibración de fase) | Sí | Cap. 10 |
| **Sincronismo por camino directo, 31/10** (`simulations/analise_sincronismo_3110`) | **NO** | Posterior a la tesis. Está **sin commit** en git. Su informe y sus figuras están en español |
| Datos reales | No (no existen todavía) | Cap. 11 = protocolo + `TODO` + `[a definir]` |

**Conclusión:** la tesis está al día **salvo el estudio de sincronismo del
31/10**, que además es justo el insumo que el orientador pide para el
"Pacote de sincronização" (§6 del plan) y para H3. Ese estudio también tiene
límites que deben declararse: cielo libre (sin refracción), reflectores
puntuales sobre el suelo, órbita circular común (no la hélice de la tesis) y sin
$H(f)$ ni efectos térmicos.

### 0.2 Plazos del orientador que ya vencieron (hoy es 10/11)

| Plazo | Entrega | Estado → acción |
|---|---|---|
| 31/10 | Inscripción como alumno regular 2027-1 | [ ] Confirmar que se hizo y guardar el comprobante. Si **no**, es la prioridad absoluta |
| 03/11 | Tabla del sistema real (con el equipo de hardware) | [ ] **Bloqueada por el laboratorio**: enviar una plantilla (§5.1) |
| 08/11 | Protocolo piloto v1 y criterios de aceptación | [ ] Fase 3 |
| 10/11 | Presupuesto de sincronización y plan de banco | [ ] Fase 2 (gran parte ya está en `analise_sincronismo_3110`) |
| 13/11 | Ensayo completo sin vuelo (del dato bruto a las métricas) | [ ] Fase 3, con datos sintéticos |
| 15/11 | Estructura del Artículo 1 e informe de campaña | [ ] Fase 4 |

### 0.3 Impacto del experimento real (sensores en coordinación)

El plan del orientador supone una **campaña piloto en noviembre** (P0–P3). Con
los sensores aún en coordinación, lo más probable es que se aplique su propio
plan de riesgo: *"Dados de novembro falham → dezembro vira integração; repetir
em janeiro"*. Este plan separa entonces todas las tareas en:

- 🟢 **Independiente de los sensores**: se puede hacer ya.
- 🟡 **Preparable**: se deja lista con datos sintéticos y se completa cuando llegue el hardware.
- 🔴 **Bloqueada**: requiere sensores o vuelo.

**Acción inmediata:** [ ] comunicar al orientador la situación del laboratorio y
proponer que P0/P1 (banco + referencia) pasen a la primera fecha posible de
diciembre o enero, sin mover el Artículo 1 (31/03).

---

## Fase 0: Trámites críticos (esta semana) 🟢

Son camino crítico de la defensa; no dependen de la ciencia.

- [ ] Verificar el estado de la inscripción 2027-1 en el sistema FEEC/DAC (selección del 04/11 al 30/11; posible entrevista).
- [ ] Auditar el historial de alumno especial: listar asignaturas, créditos y conceptos (solo cuentan A/B).
  - [ ] Calcular TC = créditos ÷ 16 × 24 meses. **≥12 créditos** es el mínimo para agosto; lo ideal son 16.
  - [ ] Dejar llenado (sin enviar) el formulario de *aproveitamento de estudos* del curso 61, para protocolarlo en **marzo de 2027**.
- [ ] Revisar el requisito de proficiencia en inglés (lectura y escritura) y cómo se acredita.
- [ ] Anotar en el calendario: 18/12 (resultado), 05–23/01 (confirmación), 01–03/02 (matrícula).

**Puerta F0:** inscripción confirmada y conteo de créditos conocido → informar al orientador.

---

## Fase 1: Reestructuración de la tesis sin datos reales 🟢 (frente de escritura)

Objetivo: dejar la tesis en el formato que pide el orientador **antes** de que
lleguen los datos, para que el Cap. 11 solo tenga que llenarse.

### 1.1 Marco científico congelado: H1–H3 (cambio estructural principal)

- [x] Sustituir las hipótesis experimentales **E1–E5** del Cap. 11 (`dissertacao_tomo_sar.tex` ≈ l. 3332–3353) por **H1 Imagen, H2 Clutter y H3 Coherencia**, con la afirmación y el "teste decisivo" tal como están en el §3 del plan.
  - [x] **Eliminar E1, la invariancia de potencia**, como hipótesis experimental: el orientador la da por cerrada ("não entra no plano de voo"). Queda como resultado de simulación en el Cap. 5.
  - [ ] Reetiquetar las estrategias de sincronismo del estudio del 31/10 (**E0–E3**) como **S0–S3**, para que no choquen con las hipótesis. Así se resuelve también la ambigüedad del "E3 duplicado". Revisar además el salto de página entre las pp. 74 y 75 del PDF.
- [x] Añadir al Cap. 1 los **seis principios ya establecidos** (§3 del plan) y la **matriz mínima de evidencia** (afirmación × evidencia numérica × evidencia real × regla).
- [x] Reescribir **O8** (≈ l. 494): "Projetar, executar e analisar a validação experimental…", ya no solo "definir o protocolo".
- [x] Reescribir el objetivo general (≈ l. 462) para incluir la validación experimental.
- [x] Hacer los objetivos O1–O8 **medibles** (métrica + umbral) y vincular cada uno a H1/H2/H3.

> **Avance 10/11/2026 (1.1):** nueva §1.4 *Hipóteses e Escopo* (principios (i)–(vi), Tab. 1.1 H1–H3 y Tab. 1.2 de la matriz de evidencia); §1.5 *Objetivos* con el objetivo general y O1–O8 vinculados a H1–H3 y con criterios verificables. En el Cap. 11, H1–H3 tienen criterios de decisión y E1 queda fuera del plan experimental. Las referencias E1–E4 del Cap. 11 se cambiaron a H1–H3. El PDF compila (113 pp.).
> Pendiente en otras subfases: renombrar S0–S3 (1.5, al integrar el sincronismo); el Resumen (l. ~197) y el Cap. 12 (l. ~3796, 3823) siguen hablando de "definir o protocolo" (1.2); la lista del plan de vuelo en el Cap. 11 sigue con 45/90/120/140/180 (3.4).

### 1.2 Nivel de evidencia en cada afirmación

- [ ] Definir una convención: `[A]` analítico, `[S]` simulado, `[E]` experimental, por ejemplo con una tabla al final del Cap. 1 y una columna en las tablas de síntesis.
- [ ] Revisar el **Resumen/Abstract** (l. 150–260) y el **Cap. 12**: quitar el lenguaje definitivo sobre resultados que solo son simulados ("demonstra-se" → "a simulação indica…").
- [ ] Separar la lista de contribuciones (§1.5, l. 497) por nivel de evidencia.
- [ ] Cambiar "trajetória ótima/proposta" por "**selecionada sob restrições**" (regla de la matriz).

### 1.3 Intervención por capítulo (tabla §2.4 del plan)

| Cap. | Decisión | Tareas concretas | ✓ |
|---|---|---|---|
| 1 | Revisar | Lo de 1.1 y 1.2; actualizar "Organização da Tese" | [ ] |
| 2 | Reducir | Mantener solo lo que se usa después. Recortar Brewster ideal y pseudo-Brewster (l. 695–800) a lo mínimo que necesita el Cap. 9 | [ ] |
| 3 | Mantener | Ampliar la comparación con la literatura (ver la matriz del estado del arte, 1.4) | [ ] |
| 4–5 | Condensar | Fundir en **un solo capítulo** (o reducir el 5 a ~⅓) que conserve las tres evidencias de la refutación y la discusión. Las figuras restantes de potencia van a un apéndice | [ ] |
| 6 | Mantener | Añadir **trazabilidad**: tabla de parámetros → configuración, versión del motor (hash de git) y pruebas de regresión. Sin nombres de scripts en el texto: va en un apéndice de reproducibilidad | [ ] |
| 7–8 | → Artículo 1 | Marcar qué figuras y tablas alimentan el A1 | [ ] |
| 9 | → Artículo 2, condicionado | Declarar el nulo VV como **cota optimista de primer orden** mientras no se valide (Fase 2.2) | [ ] |
| 10 | Dividir por función | Trayectoria → A1. Sincronización y calibración → A2. **Integrar aquí el estudio del 31/10** (1.5) | [ ] |
| 11 | Prioridad máxima | Reestructurar ya con las secciones finales vacías pero con su forma definitiva (Fase 3) | [ ] |
| 12 | Reescribir después del 15/06 | Por ahora no tocarlo, salvo para quitar las afirmaciones definitivas (1.2) | [ ] |

> **Avance 10/11/2026 (estructura, a pedido de David):** se quitaron los recuadros verdes de síntesis y las celdas coloreadas (ahora el mejor valor va en negrita). Los títulos-conclusión con `\paragraph` y las enumeraciones (i)–(iv) en negrita pasaron a prosa. Cada capítulo de resultados (7–10) termina con *Conclusões do Capítulo*. El Cap. 9 ya no tiene la tabla de preguntas y respuestas. El Cap. 12 se reescribió como *Síntese dos Resultados* por hipótesis (potencia, resolución, H1, H2, H3), sin tabla de preguntas y respuestas. "Trajetória proposta" pasó a "selecionada". **El Cap. 4 se reescribió** como *Propagação e Imageamento em Solos Estratificados*, con un análisis nuevo (`simulations/analise_multicamada/`): nueva Proposição 4.2 (resolución en un medio estratificado), error de focalización con un modelo homogéneo y transmisión TE/TM con pérdidas. El PDF compila (122 pp.).
> Siguen pendientes en 1.3: condensar el Cap. 5 y reducir el Cap. 2.

### 1.4 Ajustes editoriales (§2.3 del plan)

- [ ] Portada: "Campinas, 2026" → **2027** (l. 127).
- [ ] Figuras con rótulos en **español** dentro de la imagen: regenerar todo `figures/relatorio_1610/` (A_…, B_…, C_…, D_…, E2_…, F_…) y `figures/relatorio_2809/` con los ejes en portugués (o en inglés, pero en un solo idioma). Lo más simple es añadir un parámetro de idioma a los scripts `exp_*_fig*.py` / `exp_*.py` de `analise_1610`.
- [ ] Crear la **matriz del estado del arte** (tabla en el Cap. 1 o al inicio del Cap. 3): SAR biestático, TomoSAR, GPR aéreo, fusión multiestática y sincronización coherente × (geometría, banda, subsuelo sí/no, sincronización, fusión, validación real). Partir de `revision_literatura_bistatica.tex` y `doc/propuestas_estado_arte.md`.
- [ ] Secciones de publicaciones vacías (l. 568–571 y el apéndice de l. 3805): asignar desde ya los resultados y figuras de cada manuscrito (A1/A2), aunque todavía no estén enviados.
- [ ] `TODO` de dedicatoria y agradecimientos (l. 135, 144): se pueden dejar para el final.

### 1.5 Integrar el estudio de sincronismo del 31/10

- [ ] Hacer commit de `simulations/analise_sincronismo_3110/` (hoy no está en git).
- [ ] Nueva sección en el Cap. 10 (o en el futuro capítulo de sincronización): modelo del camino directo, CRB de retardo y fase, estrategias S0–S3, osciladores (OCXO/TCXO/XO), asociación ventana↔pulso con *stagger*, escucha discontinua y acoplamiento con la navegación.
  - Mensajes clave: **S3 (referencia de fase por ventana) ≈ reloj ideal** con SNR del directo ≳ 20 dB. La sincronización solo por retardo pierde ≈ 3 dB (300 ps contra una tolerancia de 176 ps). El modelo de reloj solo funciona con OCXO. Con PRF fija no se pueden asociar los pulsos. El error de navegación tolera ≈ 3 cm.
- [ ] Regenerar las figuras (`fig_informe.m`, `exp*.m`) con rótulos en portugués.
- [ ] Declarar sus límites: sin refracción, sin subsuelo, órbita circular común y sin $H(f)$. Esto es lo que justifica la Fase 2.1.
- [ ] Actualizar la tabla de tolerancias (`tab:traj_tolerancias`) con las tolerancias del estudio (jitter de 176 ps rms, etc.) y conciliar ambas fuentes.

**Puerta F1:** PDF compilado con H1–H3, O8 nuevo, Caps. 4–5 condensados,
portada 2027, sincronismo integrado y matriz del estado del arte → enviar al
orientador para su revisión de diciembre ("revisão dos capítulos 1–10").

---

## Fase 2: Simulaciones de cierre 🟢 (frente computacional)

Solo las que **cambian la evidencia de H1, H2 o H3** (regla del §3). Cada una
debe indicar antes de empezar qué tabla o figura va a modificar.

### 2.1 H3: presupuesto de error completo $\varphi(f,t)$ / $H(f,t)$ → Artículo 2

- [ ] Nota técnica del **presupuesto de error** (semana 1 del paquete de sincronización): relojes, disparo, $\varphi_0$, offset, deriva, jitter, posición RTK/IMU y $H(f)$. Para cada término: modelo, tolerancia, valor esperado y cuál domina.
- [ ] Extender el simulador de sincronismo (o portar sus estimadores al motor `analise_1610/bisar`) a la **geometría de la tesis**: hélice, refracción y blanco enterrado.
- [ ] Añadir **$H(f)$ relativa entre cadenas** (rizado de amplitud y fase en banda) y su calibración con un reflector de referencia más la señal directa.
- [ ] Inyectar errores conocidos y recuperarlos con los estimadores $\Delta\tau(t)$, $\varphi_0(t)$ y $H(f)$, con gráficos de error contra CRB (semana 2).
- [ ] Degradación controlada: curva de métrica de imagen frente a la fase residual, que sirva de referencia para el experimento.
- [ ] Verificar las metas del plan: fase residual < 30°, deriva ≪ 1 ns en la apertura, posición < 5 cm, $H(f)$ caracterizada.
- [ ] Plan de reserva documentado: decorrelación medida frente a la teórica, y fusión incoherente o decisional.

> Si Cristian Columna participa, esta subfase es su ventana de tres semanas
> (semana 1: presupuesto; semana 2: banco y estimadores; semana 3: integración y procedimiento operativo v1).

### 2.2 H2: robustez del nulo VV → Artículo 2

- [ ] Sustituir o complementar el SPM de 1.er orden con un **modelo de orden superior (SPM2 o IEM/AIEM)** para $\sigma^0_{VV}$ biestático.
- [ ] Aumentar el número de realizaciones de superficie y dar intervalos de confianza para el SCR/PSCR.
- [ ] Criterio de decisión: si el mínimo se **rellena**, se presenta como cota optimista y se retira la afirmación fuerte (riesgo previsto por el orientador).

### 2.3 H1: reevaluar las trayectorias finalistas → Artículo 1

- [ ] Reevaluar los finalistas de la frontera de Pareto (Exp. F) con **blancos anisótropos** (biblioteca D), **clutter** y **subaperturas parciales**, no solo con el blanco puntual.
- [ ] Comparar **mono / fija ≈105° / barrido 60°–150°** con presupuesto de muestras constante, que es el teste decisivo de H1.
- [ ] Monte Carlo de polarización con RCS anisótropa (TM co-polar frente a dual alineada).
- [ ] Factibilidad de vuelo del barrido 60°–150°: si no es factible, elegir otro candidato de la frontera de Pareto.

### 2.4 Trazabilidad y regresión (Cap. 6)

- [ ] Pruebas de regresión en `analise_1610/tests` que fijen las métricas clave (resolución, PSLR/ISLR y SCR de las configuraciones de referencia).
- [ ] Archivo de configuración versionado por figura (hipótesis → configuración → código → métrica → incertidumbre → conclusión, como en las reglas del §9).

**Puerta F2:** cada subfase cerrada con una figura o tabla nueva incorporada a la
tesis o a un artículo (regla: "toda figura nova entra imediatamente").

---

## Fase 3: Preparación experimental 🟡 (mientras se coordina con el laboratorio)

### 3.1 Insumos que se necesitan del laboratorio (enviar ya como plantilla)

Se completa la tabla del Cap. 11 (≈ l. 3372–3381), que hoy está toda en `[a definir]`:

- [ ] Radar del transceptor: $f_c$, $B$, forma de onda, PRF, potencia, ganancias, ADC y saturación.
- [ ] Receptor remoto: número de canales, **polarizaciones disponibles** (¿dual?), y si es coherente o no.
- [ ] Referencia de tiempo: oscilador (OCXO/TCXO/XO), GPSDO sí/no, PPS y enlace.
- [ ] Antenas: diagramas medidos, polarización y montaje.
- [ ] Navegación: GNSS-RTK, IMU, tasa y sistema de coordenadas.
- [ ] Plataformas: autonomía, velocidad y precisión de trayectoria (¿es factible la hélice? ¿y el barrido 60°–150°?).
- [ ] Formato de datos I/Q y metadatos que entrega el sistema.
- [ ] Fechas posibles para el banco (P0) y para el primer vuelo.

> Regla del orientador: **completar los parámetros reales antes de repetir las
> simulaciones de cierre**. Si los parámetros del laboratorio difieren mucho
> (banda, polarización), reejecutar solo las figuras que se vean afectadas.

### 3.2 Protocolo piloto v1 (P0–P3) y criterios de aceptación

- [ ] **P0 Banco:** reloj, disparo, latencia, $H(f)$ y estabilidad térmica. Salida: retardo y fase estimables, y formato de datos congelado.
- [ ] **P1 Referencia:** señal directa, reflector de referencia y trayectoria corta. Salida: foco correcto, y mono y bi alineados.
- [ ] **P2 Blanco enterrado:** mono simultáneo + 90°/105° + 135°/140°.
- [ ] **P3 Repetibilidad:** ≥3 repeticiones; ninguna conclusión con una sola pasada.
- [ ] Lista de **datos obligatorios por vuelo** (§5.1 del plan) convertida en un *checklist* de campo.
- [ ] Plantilla de **bitácora** (configuración, versiones, incidencias y criterios de descarte).
- [ ] Política de datos: los brutos son de solo lectura; los productos calibrados y procesados llevan versión.

### 3.3 "Ensayo completo sin vuelo" (entregable del 13/11)

- [ ] Generar con el simulador un **conjunto sintético con el formato del sistema real** (I/Q, metadatos y trayectorias RTK con ruido).
- [ ] Pipeline completo: lectura → sincronización (S3) → calibración $H(f)$ → BP mono y bi → fusiones (coherente igual, ponderada, incoherente y bi aislado) → métricas con incertidumbre.
- [ ] Cuando llegue el formato real, solo hay que cambiar el lector.

### 3.4 Reestructurar el Cap. 11 a su forma final

- [ ] Plan de vuelo nuevo según la matriz mínima (§5.3): **mono; 105°; 60°–150°; 135°–140°**. Opcionales: 45°, 90°, 120°, 150° y 180°. Sustituye la lista actual de 45/90/120/140/180 (≈ l. 3416).
- [ ] Blancos: referencia, compacto, tubo anisótropo y dos RCS. Superficie: zona limpia y zona con blanco, con el suelo medido.
- [ ] Subaperturas extraídas del mismo vuelo; ≥3 repeticiones (meta: 5).
- [ ] Secciones con su forma definitiva: Resultados / Incertidumbre / Fallos / Divergencias entre simulación y medida.

**Puerta F3:** protocolo v1 aprobado por el orientador, plantilla enviada al
laboratorio y pipeline sintético funcionando de punta a punta.

---

## Fase 4: Artículo 1 (TGRS) 🟢 → 🔴 (con un conjunto real calibrado)

*Geometry-Driven Coherent Monostatic–Bistatic Fusion for Subsurface UAV SAR Tomography.*

- [ ] Estructura e índice de figuras (entregable del 15/11). Fuentes: Caps. 3, 7, 8 y la parte de trayectoria del Cap. 10, más la Fase 2.3.
- [ ] Borrador completo con simulación, dejando una sección reservada para el conjunto real.
- [ ] Incorporar ≥1 conjunto real calibrado (campaña de enero o febrero).
- [ ] **Manuscrito completo: 28/02/2027. Envío: 31/03/2027.**
- [ ] Si el dato real no llega a tiempo: decidir con el orientador entre enviar a TGRS con la validación limitada que haya o cambiar a J-STARS.

---

## Fase 5: Campañas y Cap. 11 🔴 (depende de los sensores)

| Período (plan original) | Foco | Entrega | ✓ |
|---|---|---|---|
| Piloto (nov → probablemente dic/ene) | P0–P3 | Primer conjunto utilizable (G1) | [ ] |
| Enero | Integración: cerrar la sincronización; 90°/105° frente a 135°/140° | Cadena repetible (G2) | [ ] |
| Febrero | Imagen: 105° frente a 60°–150°, subaperturas y blancos iso/aniso | Figuras del A1 | [ ] |
| Marzo | Repetir los huecos críticos del A1 | Envío del A1 (G3) | [ ] |
| Abril | Clutter: barrido 120°–150°, suelo medido, polarizaciones | Mapa SCR/PSCR | [ ] |
| Mayo | Repeticiones pareadas, varios objetos y validación de la calibración | Cap. 11 + A2 v1 (G4) | [ ] |
| 1–15 jun | Contingencia (solo los fallos ya identificados) | **Data lock 15/06** | [ ] |

Para cada campaña: [ ] simulación informada por los datos reales, [ ] divergencias reportadas, [ ] resultados negativos incluidos.

---

## Fase 6: Artículo 2, calificación y defensa

- [ ] **Artículo 2 (TAES):** *Synchronization-Aware Bistatic Clutter Suppression…* Base: Fases 2.1, 2.2 y campañas de abril–mayo. Completo el 31/05; envío el 30/06.
- [ ] **Monografía de calificación** (≤20 pp. sin referencias): **30/04**. Pedido en el e-DAC: 14/05. Examen: 28/05. Antes: declaración de ética (1.er envío). Junto con el examen: proficiencia.
- [ ] Protocolar el aprovechamiento de créditos (marzo) y obtener **confirmación escrita de la CPG** (feb/mar) sobre los créditos, el plazo mínimo y la habilitación para agosto.
- [ ] **Cap. 12 reescrito** después del 15/06: qué se demostró y qué queda como limitación.
- [ ] 16–30/06: tesis completa, A2 enviado, banca definida (25/06).
- [ ] 1–9/07: Turnitin, ética (2.º envío), comprobantes de los artículos con afiliación FEEC, registro de externos y **pedido de defensa antes del 09/07**.
- [ ] 23–27/08: defensa. Hasta 60 días después: versión final (IN CCPG 06/2025, ficha catalográfica, DOI REDU).

---

## 7. Cadencia de seguimiento (§9 del plan)

- **Viernes:** memoria de una página (entrega concluida, figura nueva, bloqueo, decisión pedida y próxima meta).
- **Sábado:** reunión de 45–60 min.
- **Cada 2 semanas:** un *sprint* demostrable.
- **Fin de mes:** revisión de puerta (continuar / reducir / abandonar cada línea).

Orden sugerido de sprints, respetando las dos frentes:

| Sprint | Escritura | Computacional / experimental |
|---|---|---|
| S1 (10–23/11) | F1.1 + F1.2 (H1–H3, O8, nivel de evidencia) | F0 + F3.1 (plantilla al laboratorio) + F2.1 (presupuesto de error) |
| S2 (24/11–07/12) | F1.3 Caps. 4–5 + F1.5 sincronismo | F3.2–3.3 (protocolo + ensayo sintético) |
| S3 (08–21/12) | F1.4 editorial + matriz del estado del arte → **Puerta F1** | F2.1 ($H(f)$ + hélice) |
| S4 (ene, 1.ª quincena) | F4 estructura y borrador del A1 | F2.3 (finalistas aniso/clutter) + piloto si hay sensores |
| S5 (ene, 2.ª quincena) | Borrador del A1 | F2.2 (orden superior) + campaña de integración |
| S6 (feb) | **A1 completo, 28/02** | Campaña de imagen |

---

## 8. Registro de decisiones y bloqueos

| Fecha | Tema | Decisión / estado | Quién |
|---|---|---|---|
| 10/11/2026 | Sensores para el experimento | En coordinación con el laboratorio; el piloto de noviembre está en riesgo | David |
| | ¿Se mueve el piloto P0–P3 a dic/ene? | Pendiente | Orientador |
| | ¿Se funden los Caps. 4 y 5 en uno? | Pendiente | Orientador |
| | Idioma de las figuras (pt o en) | Pendiente | Orientador |
| | ¿Participa Cristian (3 semanas)? | Pendiente | Orientador |
| | ¿Habrá receptor de doble polarización? | Pendiente | Laboratorio |
