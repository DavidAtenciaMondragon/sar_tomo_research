# Respuesta al Prof. Hugo sobre el relatorio_1610

Prof. Hugo, muchas gracias por la lectura y por los comentarios. Respondo punto por punto, con lo que ya tengo y con lo que propongo hacer.

## 1. Validación de la trayectoria propuesta

Tiene razón. La selección del Exp. F usó la PSF del punto (ISLR) como métrica de imagen. El Monte Carlo SO(3) del Exp. E se hizo con separaciones constantes (45°, 90°, 135°, 180°), no con el barrido 60°→150°.

Propongo reevaluar no sólo la trayectoria elegida, sino los 10 a 15 mejores candidatos de la frontera de Pareto, con la misma batería para todos:

- biblioteca de objetos con orientación aleatoria SO(3);
- Rx de doble polarización alineado;
- superficie rugosa con clutter.

La métrica de imagen sería la mediana del contraste y la probabilidad de mejora, en lugar del ISLR del punto. Un resultado nuevo, que explico en el punto 2, sugiere que la elección puede cambiar: con objetos reales, la región de 135° a 140° en TM es más valiosa de lo que indicaba el punto.

## 2. El nulo TE a 90°

Tengo que corregir el relatorio en este punto. Después de entregarlo ejecuté el Exp. C2, que compara el punto escalar con una esfera y un cubo metálicos. La ventaja de TE a 90° es un artefacto del punto escalar, que no tiene polarización. El mismo desacople geométrico que anula el clutter HH a 90° también anula el eco de un objeto real.

| SCR en TE [dB] | Mono | Bi 90° |
|---|---|---|
| Punto escalar | 3.9 | 21.7 |
| Esfera | 2.3 | −4.8 |
| Cubo, orientación aleatoria | −8.8 | −9.0 |

Lo que sí sobrevive con objetos reales es el nulo VV "Brewster biestático" en TM, con Δφ entre 135° y 140°. Ahí el bi aislado mejora el SCR del cubo con orientación aleatoria en 18 a 20 dB, y el de la esfera en 8 a 11 dB, frente a la mono.

| SCR en TM [dB] | Mono | Bi 135° | Bi 140° | Fusión 140° |
|---|---|---|---|---|
| Esfera | −0.9 | 7.5 | 9.8 | −0.9 |
| Cubo nominal | −9.5 | 14.7 | 20.1 | 4.7 |
| Cubo, orientación aleatoria | −12.0 | 6.2 | 8.3 | −4.2 |

Este nulo también viene del SPM de primer orden, así que su pregunta sigue en pie. Propongo evaluarlo con:

- un modelo biestático de orden superior (IEM/AIEM o SPM de segundo orden);
- polarización cruzada de la superficie y del blanco, que hoy el modelo omite;
- superficies más realistas: correlación exponencial, rugosidad multiescala y una variación de s y ℓ.

El nulo es estrecho, de 4 a 8° en Δφ, así que también mediría cuánto se rellena y cuánta precisión de Δφ exige.

## 3. Forma de combinar los canales

Coincido en que la suma coherente con pesos iguales no es la mejor. El Exp. C2 lo muestra con claridad. Cuando el canal bi tiene mucho menos clutter, la fusión de igual peso queda dominada por el clutter de la mono y no puede mejorar el SCR más de unos 6 dB. El bi aislado a 140° supera a la fusión por más de 10 dB.

La fusión incoherente ya la evaluamos: no necesita fase, pero rinde 2 dB menos de ISLR y 4.4 dB menos de contraste que la coherente calibrada.

Propongo comparar, sobre el mismo banco de pruebas (punto, objetos SO(3) y clutter), cuatro métodos:

- fusión coherente con pesos iguales, como referencia;
- fusión ponderada, con pesos por canal elegidos para maximizar el SCR o minimizar el ISLR;
- reconstrucción conjunta por inversión regularizada, que además puede estimar la fase entre canales como incógnita;
- fusión incoherente.

## 4. Calibración entre canales

De acuerdo: una única fase constante es el caso ideal. El modelo de error realista sería

```
φ(f, t) = φ₀ + 2π f · Δτ(t) + H(f)
```

Cada término es una fuente distinta:

- Δτ(t) es el retardo de sincronismo variable en la apertura, incluyendo la deriva de reloj y el error de posición del Rx;
- H(f) es la respuesta en frecuencia distinta de las dos cadenas de recepción.

Las simulaciones ya indican la escala de los términos variables. Una deriva lineal de 1 ns en la apertura degrada el PSLR a −6.5 dB, y una deriva de posición de 20 cm rms destruye la ganancia.

Propongo dos pasos:

- simular la sensibilidad a cada término por separado;
- estimarlos pulso a pulso con la señal directa Tx→Rx, cuya geometría se conoce por RTK, y H(f) con un reflector de referencia, evaluando el error residual con una SNR realista.

## 5. Aperturas parciales y suelos húmedos

Todos los resultados usan las dos vueltas completas de la hélice. Con eso el soporte en k es anular y la PSF casi invariante en la escena. Con un arco parcial esto deja de valer, y el Δφ óptimo puede cambiar. Propongo repetir los Exp. B y F con arcos de 45°, 90°, 180° y 360°.

Sobre los suelos húmedos, el resultado es preocupante pero no depende de la trayectoria.

| Suelo | Atenuación [dB/m] | Pérdida de ida y vuelta a 2 m [dB] |
|---|---|---|
| S1 | 7.7 | ≈31 |
| S2 | 13.6 | ≈54 |
| S3 | 27.7 | ≈111 |

Esta pérdida afecta igual a la mono y al bi. Es un límite del sistema, no de la geometría. Propongo dos análisis:

- calcular la profundidad máxima con SCR > 0 para cada suelo;
- ver cómo cambia esa profundidad al bajar la frecuencia, ya que la atenuación crece aproximadamente con f.

Así podremos acotar el dominio de validez de la propuesta en lugar de presentarla como general.

## Prioridad

Si le parece bien, priorizaría el punto 1 junto con el 3, porque definen la trayectoria y el método de fusión, y luego el 2, que es lo que justificaría la novedad del artículo.

Saludos,
David

---

## Notas internas (borrar antes de enviar)

- Las cifras del punto 2 vienen del Exp. C2 (`simulations/io/analise_1610/expC2/expC2_scr_vs_dphi.csv`) y no están en el PDF entregado. Conviene agregarlas al relatorio y corregir las secciones que todavía presentan TE a 90° como el mejor resultado del estudio.
- Las pérdidas de ida y vuelta del punto 5 se calcularon con la atenuación del relatorio y dos trayectos de 2 m. Confirmarlas antes de enviar.
- La prioridad del cierre es una sugerencia. Quitarla si se prefiere que la decida el asesor.
