# Revisión de tesis doctoral — informe de jurado/asesor

**Documento revisado:** `dissertacao_tomo_sar.tex` (60 páginas compiladas, 7 capítulos + 2 apéndices)
**Título:** *Otimização da Trajetória do Receptor em Radar de Abertura Sintética Biestático com Voo Helicoidal*
**Fecha de la revisión:** 2026-08-03
**Método:** lectura completa del `.tex` (todas las secciones, ecuaciones, tablas), verificación algebraica y numérica de las derivaciones centrales, comparación cruzada de cada `\includegraphics` contra `doc/figures/`, e inspección visual de las figuras clave (renderizadas a PNG).

> Nota de transparencia: en esta misma sesión de trabajo yo (el asistente) reordené los Capítulos 4↔5, agregué la Sección 5.6 (Algoritmo Numérico) y la Sección 5.8 (estudio de sensibilidad α/grade). Esos agregados también fueron auditados aquí con el mismo rigor que el resto — no se les dio un trato especial — y dos de los hallazgos críticos (§1.2 y la recomendación de rango de α en 5.8) involucran texto que yo mismo escribí y que ahora debe corregirse.

---

## Resumen ejecutivo

La tesis tiene una **columna vertebral teórica sólida y bien encadenada**: Maxwell → Fresnel TM → Brewster → Fermat/Snell → hélice → radar biestático (Cap. 2); k-espacio → δz/δxy (Cap. 3); ABCD multicapa (Cap. 4); costo multiobjetivo + optimizador (Cap. 5); validación por backprojection (Cap. 6). La idea central (explotar el ángulo de Brewster como guía geométrica para la trayectoria del Rx, sujeto a una restricción de resolución) es clara, original y está bien motivada.

Sin embargo, en su estado actual **no está lista para banca**. Encontré:

- **3 inconsistencias numéricas/matemáticas verificables** que un jurado técnico detectará en minutos (una de ellas es un error real en una fórmula "boxed" del Capítulo 4, confirmado con contraejemplo numérico).
- **Estado del arte ausente** en los 4 capítulos de desarrollo (2, 3, 4, 5) — son placeholders `[Estado da arte a ser desenvolvido.]`.
- **Dos capítulos completos (4 y 5) sin validación** contra el simulador de backprojection — el propio Capítulo 6 lo admite en su conclusión.
- La lista de contribuciones de la Introducción y la de las Conclusiones **no coinciden entre sí**.

Ninguno de estos problemas es difícil de corregir individualmente, pero juntos son exactamente el tipo de cosas que un jurado usa para cuestionar el rigor general del documento. Prioriza la sección "Hallazgos críticos" antes que cualquier otra cosa.

---

## 1. Hallazgos críticos (bloquean cualquier defensa)

### 1.1 — Error matemático confirmado en el formalismo ABCD multicapa (Capítulo 4)

**Dónde:** Ec. (4.5) `eq:ABCD_matrix` y Ec. (4.10)–(4.11) `eq:Zin_Zout`, que definen
$$Z_i = \frac{n_i}{\sqrt{n_i^2-p^2}}$$
y de ahí $Z_{\text{in}}$, $Z_{\text{out}}$, usadas en $\Gamma_{\text{TM}}$ (4.12) y $T_{\text{TM}}$ (4.15).

**El problema:** con esa definición de $Z_i$, el **Lema 4.1** (§4.5, `lema:reducao_fresnel`) — que afirma que para $L=1$ (una sola interfaz) el modelo ABCD se reduce exactamente al coeficiente de Fresnel TM del Capítulo 2 — **no se cumple**. Verifiqué esto de tres formas independientes:

1. **Álgebra a mano:** sustituyendo $Z_i=n_i/\sqrt{n_i^2-p^2}$ en $\Gamma=(Z_{\text{out}}-Z_{\text{in}})/(Z_{\text{out}}+Z_{\text{in}})$ para $L=1$, se obtiene $\Gamma=(\cos\theta_i-\cos\theta_t)/(\cos\theta_i+\cos\theta_t)$ — **no** es $r_{\text{TM}}=(n_2\cos\theta_i-\cos\theta_t)/(n_2\cos\theta_i+\cos\theta_t)$ como afirma la prueba del Lema (falta el factor $n_2$).
2. **Contraejemplo numérico** ($n_1=1$, $n_2=2$, $\theta_i=30°$): $r_{\text{TM}}$ estándar (Fresnel, Cap. 2) da $R_{\text{TM}}=0{,}0800$. La fórmula del Cap. 4, tal como está escrita, da $R_{\text{TM}}=0{,}0031$ — **26 veces más chica**, y con signo de $\Gamma$ opuesto.
3. **Comparación con el código real** (`simulations/common/multislabTMcoef.m`): el código usa $Z_i=\eta_0\cos\theta_i/n_i$ (equivalente a $n_i/\cos\theta_i$ salvo por una dualidad recíproca), que **sí** reproduce $R_{\text{TM}}=0{,}0800$ correctamente, tanto para el caso de 1 interfaz como para un caso de 2 capas que probé aparte.

**Buena noticia:** el código MATLAB usado para todas las simulaciones (incluido el estudio de sensibilidad del §5.8) está **correcto** — el error está únicamente en cómo quedó transcrita la fórmula en el `.tex`. La corrección sugerida es
$$Z_i = \frac{n_i}{\cos\theta_i} = \frac{n_i^2}{\sqrt{n_i^2-p^2}}$$
(falta un factor $n_i$ en el numerador de la Ec. 4.5). Con esta corrección, verifiqué numéricamente que el Lema 4.1 sí se cumple exactamente.

**Por qué es crítico:** el Lema 4.1 es la única verificación formal de consistencia que el Capítulo 4 ofrece para su resultado central ($\Gamma_{\text{TM}}$, $T_{\text{TM}}$ vía ABCD), y tal como está escrita, **la prueba no cierra**. Cualquier jurado con formación en electromagnetismo que intente reproducir la prueba del Lema (es un cálculo de 4 líneas) lo va a notar. Hay que corregir $Z_i$ en las Ecs. (4.5), (4.10)–(4.11), reverificar que $\Gamma_{\text{TM}}$/$T_{\text{TM}}$ (4.12)–(4.15) sigan siendo internamente consistentes, y revisar el Apéndice A.3 (que actualmente solo remite a la Ec. 4.5 sin re-derivar $Z_i$ de forma independiente, por lo que hereda el mismo error).

### 1.2 — La Fig. 5.3 (frontera de Pareto) contradice al texto que la describe

**Dónde:** §5.4 "Análise da Fronteira de Pareto" (`sec:pareto`), Fig. `fig:pareto` (`pv_fig07_pareto.pdf`).

El texto afirma: *"O joelho da curva de Pareto [...] define o intervalo ótimo de compromisso $\alpha^*\in[0{,}2,\,0{,}4]$."* Pero al renderizar la figura que se incluye justo debajo, esta marca explícitamente — con un rombo en el panel (a) y una línea vertical punteada en el panel (b), ambos con la etiqueta **"$\alpha^*\approx0{,}87$"** — un punto óptimo completamente distinto. La figura y el texto que la acompaña **no dicen lo mismo**, y no hay ninguna nota que explique la diferencia (¿son criterios distintos — "codo de máxima curvatura" vs. "saturación de la ganancia de resolución"? Si es así, hay que decirlo explícitamente).

**Esto se propaga:** el §5.8 (estudio de sensibilidad multicapa que yo agregué esta sesión) cita esta misma figura para justificar que su recomendación de $\alpha\in[0{,}3,\,0{,}6]$ es *"consistente com o intervalo $\alpha^*\approx0{,}2$–$0{,}4$"* — heredando la inconsistencia sin haberla verificado contra la figura real. Esto también hay que corregirlo.

**Acción sugerida:** decidir cuál es el criterio real usado para identificar el "codo" (¿máxima curvatura sobre $(P_r,\,\delta_{xy}\delta_z)$? ¿otro?), regenerar la figura o el texto para que coincidan, y solo entonces volver a redactar la comparación en §5.8.

### 1.3 — La Tabla 6.2 y la Figura 6.2 reportan errores de validación distintos para los mismos casos

**Dónde:** §6.4.2, Tabla `tab:validation` vs. Figura `fig:validation` (`fig12_validation.pdf`).

| Caso | Tabla 6.2 | Figura 6.2 |
|---|---|---|
| $\delta_z$, $\Delta\phi=180°$ eixo | 0,05% | 0,1% |
| $\delta_z$, $\Delta\phi=90°$ eixo | 0,05% | 0,1% |
| $\delta_z$, off-axis | 7,1% | 6,6% |
| $\delta_{xy}$, on-axis | 7,8% | 8,4% |
| $\delta_{xy}/\delta_x$, off-axis | 11,6% | **13,2%** |

Ningún par coincide. El caso más delicado es el último: el **Resumo/Abstract** afirma en su primera página *"erros inferiores a [...] 12% na resolução horizontal $\delta_{xy}$"* — una cifra que toma prestada de la Tabla (11,6% < 12%), pero que **la propia Figura 6.2 del mismo capítulo contradice** (13,2% > 12%). Un jurado que lea el resumen y luego mire la figura va a encontrar que el resultado más citado de toda la tesis no se sostiene de forma inequívoca con los propios datos presentados.

**Acción sugerida:** regenerar tabla y figura desde la misma corrida de datos (parece que se generaron por separado, quizás en momentos distintos del desarrollo), y ajustar la cifra del Resumo/Abstract al valor verificado y consistente.

---

## 2. Hallazgos mayores (necesarios para un documento defendible)

### 2.1 — "Estado da arte" ausente en todos los capítulos de desarrollo

Los Capítulos 2, 3, 4 y 5 tienen todos una sección `\section{Introdução e Estado da Arte}` que es literalmente `\textit{[Estado da arte a ser desenvolvido.]}`, con un comentario TODO listando referencias por citar. Esto no es un detalle menor: sin revisión bibliográfica, el jurado no tiene cómo evaluar la brecha que la tesis dice llenar (la introducción general sí la enuncia — "nenhum trabalho utilizou o ângulo de Brewster como guia geométrico..." — pero eso debe sostenerse capítulo por capítulo, con citas). Ya hay ~20 referencias en la lista final; el trabajo pendiente es integrarlas (y las ~20-30 que faltan, según el propio TODO de la lista de referencias) en cada `Estado da arte`, no generarlas desde cero.

### 2.2 — Secciones de resultados sin completar en capítulos centrales

- §2.4.2 "Modelo de Sinal e Forma de Onda" y §2.4.4 "Algoritmo de *Backprojection* Biestático" son TODOs vacíos — pero el Capítulo 6 completo depende de un "simulador de backprojection" cuyo algoritmo nunca se describe en el cuerpo de la tesis. Un jurado preguntará cómo funciona el simulador que produjo la Tabla 6.2/Fig. 6.2 si el capítulo de fundamentos no lo explica.
- §3.7 "Resultados de Simulação" (Cap. 3, resolução espacial) es un TODO vacío, aunque su contenido ya existe de facto en el Capítulo 6 (Tabla 6.2). Conviene decidir: o se llena §3.7 con un adelanto/resumen y referencia cruzada a Cap. 6, o se elimina la sección duplicada.
- §4.6 "Resultados de Simulação" (Cap. 4, multicamadas) es un TODO vacío: **no hay ninguna curva $T_{\text{TM}}$ vs. $\theta_0$, ni verificación numérica de $R+T=1$, en toda la tesis** para el modelo multicapa — a pesar de que este modelo es la Contribución #4 explícita de la Introducción. Esto es exactamente el tipo de verificación que habría detectado el error de §1.1 antes de escribir el Lema.

### 2.3 — El propio Capítulo 6 admite que los Capítulos 4 y 5 no están validados

§6.5 "Validação do Modelo de Otimização de Trajetória" (Cenário B) y §6.6 "Validação do Modelo Multicamadas" (Cenário C) son ambos TODOs, y la §6.8 "Conclusões do Capítulo" lo dice explícitamente: *"As validações dos modelos de otimização (Cenário B) e multicamadas (Cenário C) são identificadas como requisitos de trabalho futuro."* En términos de peso de la tesis, esto significa que **2 de los 5 capítulos de desarrollo (Caps. 4 y 5) —y las Contribuições #3 y #4 de la Introdução— no tienen ninguna validación independiente** contra el simulador de backprojection, solo consistencia interna del propio optimizador. Esto está directamente conectado con 1.1: si Cenário C se hubiera ejecutado, el error de $Z_i$ probablemente se habría detectado por discrepancia con la simulación.

### 2.4 — Las listas de "contribuições" de la Introdução y de las Conclusões no coinciden

| # | Introdução §1.3 | Conclusões §7.1 |
|---|---|---|
| 1 | Fresnel TM + Brewster | Fresnel TM + Brewster ✓ |
| 2 | Resolución 3D | Resolución 3D ✓ |
| 3 | Función de costo log. | **Formulação de Fermat/Newton-Raphson** (no estaba en la lista de §1.3) |
| 4 | ABCD $N$ camadas | Función de costo log. |
| 5 | **Validación por backprojection** (con cifras de error) | ABCD $N$ camadas |

La Conclusão reemplaza la contribución de "validación por backprojection" (que sí tiene resultados concretos, Tabla 6.2) por "formulação de Fermat" (que es más un fundamento metodológico que una contribución original independiente, y que en el Resumo se presenta como uno de los "quatro pilares", no como una "contribuição" en sí). Además, **ninguna de las dos listas** menciona el algoritmo de optimización numérica multi-arranque (§5.6) ni el estudio de sensibilidad α/grade multicapa (§5.8) — ambos añadidos sustanciales de contenido que hoy no están acreditados como contribución en ningún lado del documento. Recomiendo unificar ambas listas (deberían decir *exactamente* lo mismo, es estándar en una tesis) y evaluar si el algoritmo numérico y el estudio de sensibilidad ameritan ser una 6ª contribución explícita.

### 2.5 — Capítulo 2 muy escasamente ilustrado pese a tener material gráfico ya generado

Encontré **12 figuras generadas** (`fig01_radar_basics`, `fig02_real_vs_synth`, `fig03_circular_sar`, `fig04_helical_params`, `fig05_snell_geometry`, `fig08_full_system`, `fig09_psf_comparison`, `fig11_offaxis_correction`, `pv_fig03_snell_1d`, `pv_fig04_bistatic_radar`, `pv_fig05_cost_function`, `pv_fig09_alpha_sensitivity`) que existen en `doc/figures/` pero **no están referenciadas en ningún lugar del `.tex`**. Mientras tanto, secciones enteras del Capítulo 2 (§2.1 historia del SAR, §2.3 traçado de Fermat, partes de §2.4) no tienen ni una figura. Vale la pena revisar cada una de las 12 y decidir: insertarla donde corresponda, o borrarla del repositorio si quedó obsoleta (algunas, como `pv_fig09_alpha_sensitivity`, muestran el mismo `α*≈0,87` del hallazgo 1.2 y podrían ayudar a resolver esa inconsistencia si se entiende primero de dónde sale ese número).

---

## 3. Hallazgos menores / pulido

| # | Hallazgo | Dónde |
|---|---|---|
| 3.1 | Citas manuales `[N]` en vez de `\cite{}`+BibTeX. Hoy los ~20 números coinciden con la lista (los verifiqué), pero el TODO de la bibliografía pide expandir a 40–50 referencias — renumerar a mano 40+ citas dispersas en 2200 líneas es una fuente segura de errores futuros. | En todo el documento |
| 3.2 | El acrónimo **VANT** se usa en el cuerpo del texto (§5.7, "trajetórias de VANT") pero no está en la Lista de Acrônimos. | `Lista de Acrônimos`, §5.7 |
| 3.3 | Proposição de §5.4 dice que para $\alpha=1$ la solución tiende a *"mínimo $\psi_0$"* — pero $\delta_{xy}\propto 1/\sin\psi_0$ (Ec. 3.19), por lo que minimizar $\psi_0$ empeoraría $\delta_{xy}$, no lo mejoraría. Puede haber una razón geométrica válida (p. ej. relacionada a $\delta_z$/$R_0$, no a $\delta_{xy}$ directamente) pero el texto no la explica — vale la pena revisar o aclarar. | §5.4, Proposição |
| 3.4 | La constante del criterio $-3$dB de $J_0$ se cita como $u_{-3}\approx1{,}20$; un recálculo rápido por interpolación da $u\approx1{,}13$. Diferencia pequeña (~6%) pero afecta el coeficiente "0,60" de la Ec. 3.20 — vale la pena resolver $J_0(u)=1/\sqrt2$ numéricamente y confirmar la cifra exacta. | §3.5, Ec. `eq:delta_xy` |
| 3.5 | Los ejemplos numéricos ilustrativos usan índices de refracción distintos en cada capítulo ($n_2=2$ en Cap. 2; $n_1=\sqrt2,n_2=2$ en Cap. 4 §4.6; $n_1=2,n_2=3$ en Cap. 5 §5.8; $n_1=1{,}5,n_2=2{,}5$ en los TODOs de Caps. 4 y 6) sin un escenario de referencia único declarado. No es un error, pero definir **un** escenario "canónico" (p. ej. "solo seco típico, $\varepsilon_r=4$") y usarlo consistentemente, marcando las variaciones como estudios paramétricos explícitos, mejoraría la cohesión de lectura. | Varios capítulos |
| 3.6 | §2.5.1 afirma *"A condição ótima de inclinação é $\beta=\psi_0$"* sin prueba ni referencia adelantada — la demostración real aparece recién en el Cap. 3 §3.6. Conviene agregar *"(demonstrado na Seção~\ref{sec:helix_opt})"*. | §2.5.1 |
| 3.7 | Agradecimentos, §1.5 "Publicações" y Apéndice B "Lista de Publicações" están vacíos (`TODO`). Esperable en un borrador, pero son bloqueantes para el depósito final. | Front matter, Apéndice B |
| 3.8 | El Resumo/Abstract y la lista de contribuciones no mencionan el algoritmo multi-arranque (§5.6) ni el estudio de sensibilidad multicapa (§5.8), pese a ser contenido sustancial y original del Capítulo 5. | Resumo/Abstract, §1.3 |

---

## 4. Lo que ya funciona bien (no tocar)

- La derivación de Fresnel TM (§2.2) y la conservación de energía (Apéndice A.1) son correctas — las verifiqué algebraicamente paso a paso.
- La prueba del ángulo de Brewster (Teorema 2.1) y su corolario de ortogonalidad son correctos.
- La formulación variacional de Fermat → ley de Snell (§2.3) y la prueba de unicidad vía $g'(u)>0$ son correctas.
- La fórmula de resolución vertical $\delta_z$ (Ec. 3.13/`eq:Wz`) **está numéricamente validada**: reconstruí a mano el desglose 91 MHz (ancho de banda) + 620 MHz (apertura tomográfica) = 711 MHz de la Fig. 3.1 usando los parámetros de la Tabla 6.1, y coincide con el $\delta_z=21{,}10$ cm reportado — el modelo es internamente consistente y consistente con la simulación.
- La estructura general (fundamentos → resolución → multicamadas → optimización → validación → conclusiones) es lógica y la reordené para reflejar la dependencia real (Cap. 4 antes que Cap. 5) — el compilado con `pdflatex` corre limpio, dos pasadas, sin errores ni referencias indefinidas.
- Las "Limitações do Modelo Atual" y "Perspectivas Futuras" (Cap. 7) están bien pensadas y apropiadamente delimitadas para una tesis doctoral — no son genéricas, apuntan a extensiones concretas (medios con pérdidas, Cramér-Rao, rugosidad de Kirchhoff, validación de campo, multi-alvo).

---

## 5. Checklist priorizado

1. [ ] Corregir $Z_i$ en Ecs. (4.5)/(4.10)-(4.11) del Cap. 4 y reverificar el Lema 4.1, $\Gamma_{\text{TM}}$, $T_{\text{TM}}$ y el Apéndice A.3.
2. [ ] Reconciliar Fig. 5.3 (`pv_fig07_pareto.pdf`, marca $\alpha^*\approx0{,}87$) con el texto ($\alpha^*\in[0{,}2,0{,}4]$) en §5.4 y en la recomendación de §5.8.
3. [ ] Reconciliar Tabla 6.2 vs. Figura 6.2 (errores de validación) y ajustar la cifra "<12%" del Resumo/Abstract al valor correcto y consistente.
4. [ ] Unificar la lista de contribuciones entre §1.3 y §7.1; decidir si el algoritmo §5.6 y el estudio §5.8 se acreditan como contribución adicional.
5. [ ] Completar "Estado da arte" en Caps. 2, 3, 4, 5 con las referencias ya listadas.
6. [ ] Completar §2.4.2 (modelo de señal), §2.4.4 (algoritmo BP), §4.6 ($T_{\text{TM}}$ vs. $\theta_0$), y los Cenários B/C del Cap. 6 — con prioridad en el Cenário C, que habría detectado el hallazgo 1.1 antes.
7. [ ] Revisar las 12 figuras generadas pero no usadas; insertar o descartar.
8. [ ] Pulido menor: VANT en acrónimos, migrar citas a `\cite`/BibTeX, verificar constante $-3$dB de $J_0$, aclarar Proposição de §5.4, unificar ejemplos numéricos.

---

*Este documento es una herramienta de trabajo para el autor, no parte de la tesis. Todas las verificaciones numéricas (Z_i, α de Pareto, δz) se hicieron independientemente a partir del texto y el código fuente del repositorio; los scripts de verificación no se conservan en el repositorio.*
