# Método de cálculo

## Geometría y convenciones

- Deformación plana, por metro de ancho de talud. Unidades: m, kPa, kN/m³ y grados.
- El terreno es una polilínea con x creciente. El talud desciende hacia la izquierda: el pie
  queda a la izquierda, la corona a la derecha y la masa se desliza en el sentido −x.
- Un talud simple se define por altura H y ángulo β: plataforma al pie en y = 0, cara inclinada
  y corona horizontal en y = H.
- El nivel freático es otra polilínea. En un talud simple puede darse en forma relativa: aflora a
  `altura_en_pie` sobre la plataforma del pie y sube en línea recta hasta quedar a
  `profundidad_bajo_corona` bajo la corona, a `distancia_tras_corona` detrás del borde.

## Superficie de falla circular

Cada círculo se define por la abscisa de **salida** x₁ (cerca del pie), la de **entrada** x₂
(cerca de la corona) y el **semiángulo central** ψ del arco. Con la cuerda P₁P₂ de largo L:

    R = (L/2) / sen ψ        d = (L/2) / tan ψ        centro = punto medio + d · n

donde n es la normal unitaria a la cuerda hacia arriba. Un ψ pequeño da un arco casi plano y
uno grande, un arco profundo. Una superficie es admisible si:

1. los dos extremos quedan por debajo del centro, de modo que el arco es una función de x;
2. el arco queda bajo el terreno entre x₁ y x₂ (una sola masa deslizante);
3. no pasa por debajo del estrato firme, si se definió uno;
4. en Bishop, la solución converge y mα ≥ 0.2 en todas las dovelas (criterio práctico habitual
   para descartar soluciones con dovelas de salida muy empinadas).

## Dovelas

La masa se divide en n dovelas verticales del mismo ancho b (100 por defecto). Para cada dovela,
con la base en el centro de la dovela:

- espesor h = cota del terreno − cota de la base;
- altura saturada h_w = cota del nivel freático (recortado al terreno) − cota de la base, entre 0 y h;
- peso W = b · [γ · (h − h_w) + γ_sat · h_w];
- inclinación de la base: sen α = (x − x_c)/R, positiva donde la base sube hacia la corona;
- presión de poros hidrostática u = γ_w · h_w. Se usa la altura vertical del nivel freático sin
  corrección por su inclinación, lo que queda del lado conservador.

## Factor de seguridad

**Fellenius (método ordinario)**: equilibrio de momentos respecto del centro, sin fuerzas entre dovelas.

    FS = Σ [c'·l + (W·cos α − u·l)·tan φ'] / Σ W·sen α,     l = b / cos α

**Bishop simplificado**: agrega el equilibrio vertical de cada dovela y desprecia el corte entre
dovelas (cumple el equilibrio de momentos, no el horizontal).

    FS = Σ { [c'·b + (W − u·b)·tan φ'] / mα } / Σ W·sen α,     mα = cos α + sen α · tan φ' / FS

Como FS aparece en ambos lados, se resuelve por iteración de punto fijo a partir del valor de
Fellenius, hasta un cambio relativo menor que 10⁻⁹. En ambos métodos la fuerza normal efectiva se
limita a valores no negativos.

## Búsqueda del círculo crítico

1. **Malla**: 25 salidas × 25 entradas × 20 valores de ψ (3° a 80°), evaluados con 50 dovelas.
   Por defecto, las salidas van desde 1.5·H antes del pie hasta 3/4 de la cara y las entradas
   desde 1/4 de la cara hasta 2·H detrás de la corona.
2. **Refinamiento**: desde los 4 mejores puntos de la malla, separados entre sí, se aplica una
   búsqueda por patrones con 100 dovelas. Se prueba ±paso en cada parámetro y el paso se reduce a
   la mitad cuando ninguna dirección mejora, hasta 0.001·H en posición y 0.02° en ψ.

Todo el cálculo está vectorizado con numpy: la malla completa (unos 12 000 círculos) se evalúa en
décimas de segundo.

## Análisis probabilístico

**Variables aleatorias.** c' y φ', con distribución normal o lognormal (media y COV). La
correlación ρ(c', φ') se impone en el espacio normal estándar (cópula gaussiana): se generan
normales estándar independientes, se correlacionan con la descomposición de Cholesky y se
transforman a cada distribución. Para una variable lognormal, la correlación indicada es la de
ln(c'). Las normales se recortan en 0 y el número de valores recortados se informa.

**Monte Carlo.** Para cada muestra se calcula el FS de Bishop en varias superficies y se toma el
mínimo, porque la superficie crítica cambia con c' y φ' (con c' baja los círculos críticos son más
someros). Las superficies son:

- la crítica determinista;
- las críticas de cuatro búsquedas completas en los puntos μ ± σ de c' y φ', con tres
  candidatas de cada una;
- diez círculos de bajo FS y distintos entre sí de la malla determinista.

También se informa el resultado sobre la superficie crítica fija, que es el método «global
minimum» de los programas comerciales, como comparación.

**Probabilidad de falla.** P(FS < 1) = número de muestras con FS < 1 / N, con intervalo de
confianza de Wilson al 95 % (válido incluso sin fallas) y error estándar √(p(1−p)/N).

**Índices de confiabilidad.**

    β normal    = (μ_FS − 1) / σ_FS
    β lognormal = ln( μ_FS / √(1 + V²) ) / √( ln(1 + V²) ),    V = σ_FS / μ_FS

**Control con Rosenblueth (1975).** Se evalúa el FS en los 2ᵏ puntos μ ± σ con pesos
P = [1 + Σ_{i<j} sᵢ·sⱼ·ρᵢⱼ] / 2ᵏ, y se comparan la media y la desviación con las de Monte Carlo
sobre la misma superficie fija.

## Ángulo de diseño

El comando `diseno` repite el análisis completo para cada ángulo de la cara y recomienda el mayor
ángulo tal que ese ángulo y todos los más tendidos cumplen FS ≥ FS mínimo y P(FS < 1) ≤ PF máxima.
La resolución está dada por el paso del barrido (2.5° por defecto).
