# analise_sincronismo_3110: sincronismo BiSAR por camino directo

Viabilidad de sincronizar un receptor pasivo (Radar 2) con un Tx/Rx
monostático en dron (Radar 1) **sin GPS, PPS ni reloj común**, usando
solo el camino directo Tx→Rx, y su impacto cuantitativo en una imagen
BiSAR formada por Backprojection.

* **Informe con figuras:** [`informe_sincronismo.tex`](informe_sincronismo.tex)
  (PDF compilado en `build/informe_sincronismo.pdf`; `latexmk -pdf -outdir=build informe_sincronismo.tex`).
* El análisis teórico detallado, con todas las tablas, está en [`analise_teorica.md`](analise_teorica.md).
* La simulación es MATLAB puro (R2021a o posterior, sin toolboxes).

## Reproducción

```matlab
cd simulations/analise_sincronismo_3110
run_all            % ~25 min; o cada exp*.m por separado
tests/test_cadena  % prueba de humo (~10 s)
```

Los resultados se escriben en `resultados/` (CSV, MAT, logs) y las figuras
en `resultados/figs/` (PNG + PDF vectorial).

## Experimentos

| Script | Pregunta que responde |
|---|---|
| `exp1_tolerancias.m` | ¿Qué precisión **hace falta**? Errores inyectados (offset, frecuencia, deriva cuadrática, fase blanca, jitter) frente a la teoría |
| `exp2_cadena_completa.m` | Demostración de extremo a extremo: reloj, asociación, imágenes E0–E3, monostática y fusión mono+bi |
| `exp3_snr_oscilador.m` | ¿Qué precisión **se alcanza**? Monte Carlo frente a SNR del directo (vs CRB) y clase de oscilador |
| `exp4_pulsos_perdidos.m` | Escucha discontinua: separa la pérdida de sincronismo de la de muestreo; puenteo de fase entre huecos |
| `exp5_asociacion.m` | Observabilidad del índice de pulso: stagger y número de ventanas |
| `exp6_navegacion.m` | Acoplamiento entre el error de navegación del Tx y el sincronismo por camino directo |
| `fig_informe.m` | Figuras adicionales del informe (escenario, señal del Rx2, reloj, cortes, resumen de conclusiones); lee los CSV |

## Estrategias de sincronismo comparadas (`lib/sincroniza.m`)

| | Qué usa del camino directo | Fase aplicada al eco |
|---|---|---|
| E0 offset | retardo y fase de la 1.ª ventana | constante |
| E1 modelo | retardos de todas las ventanas → ajuste `δ0 + εt + γt²/2` | `2π fc δ̂(t) + φ̂0` |
| E2 tiempo | retardo de cada ventana | `2π fc δ̂_j + φ̂0` (derivada del retardo) |
| E3 fase | retardo **y fase portadora** de cada ventana | `−(arg r_j + 2π fc τ_d)` (referencia por pulso) |
| Ideal | verdad terreno (solo para evaluar) | exacta |

## Estructura

| Ruta | Contenido |
|---|---|
| `lib/params_sinc.m` | Todos los parámetros (banda, PRF y alturas de `analise_1610/params/common.json`) |
| `lib/osc_clase.m` | Modelos de oscilador OCXO / TCXO / XO / XO_barato |
| `lib/simula_adquisicion.m` | Log del Tx, reloj del Rx2, síntesis de ventanas, detección + ML del directo, cancelación, compresión |
| `lib/estima_directo.m` | Detección del primer camino y estimador ML (Newton) de retardo y fase |
| `lib/asocia_pulsos.m` | Asociación ventana↔pulso por LS sobre hipótesis de índice |
| `lib/sincroniza.m` | Estrategias E0–E3 |
| `lib/bp_bistatico.m` | Backprojection con base de tiempos y fase corregidas |
| `lib/metricas_imagen.m` | Desplazamiento, pérdida de pico, PSLR, ISLR y coherencia frente al reloj ideal |
| `tests/test_cadena.m` | Prueba de humo |
