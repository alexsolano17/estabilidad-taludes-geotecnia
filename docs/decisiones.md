# Decisiones de diseño

Registro breve de las decisiones de arquitectura (ADR): contexto, decisión y consecuencias.

## ADR-001: Bishop simplificado como método principal y Fellenius como control

- **Contexto:** el proyecto necesita un método de equilibrio límite preciso para superficies
  circulares, que pueda verificarse con casos publicados y sea rápido para Monte Carlo.
- **Decisión:** usar Bishop simplificado. Fellenius se calcula siempre como referencia
  conservadora y para explicar la diferencia entre ambos.
- **Consecuencias:** para círculos, Bishop da resultados cercanos a métodos rigurosos (ACADS 1(a):
  0.985 frente a 0.986 de Spencer en Slide2) con una iteración barata. No cumple el equilibrio
  horizontal ni sirve para superficies no circulares.

## ADR-002: Círculos definidos por salida, entrada y semiángulo

- **Contexto:** una malla de centros y radios genera muchos círculos que no cortan el talud o lo
  cortan de forma inadecuada.
- **Decisión:** parametrizar cada círculo por sus puntos de salida y entrada sobre el terreno y
  el semiángulo del arco (búsqueda «entrada-salida»).
- **Consecuencias:** todos los círculos generados pasan por el terreno y los rangos tienen sentido
  físico (pie y corona). Los círculos con entrada por encima del centro se descartan.

## ADR-003: Malla más refinamiento por patrones, sin bibliotecas de optimización

- **Contexto:** la función FS(x₁, x₂, ψ) tiene zonas no admisibles (FS = ∞) y puede tener
  mínimos locales.
- **Decisión:** recorrer una malla gruesa y refinar los 4 mejores puntos distintos con una
  búsqueda por patrones de paso decreciente.
- **Consecuencias:** es robusta frente a zonas no admisibles y no depende de scipy. Una malla más
  fina da el mismo mínimo (diferencia < 0.002).

## ADR-004: c' y φ' como variables aleatorias, correlacionadas mediante una cópula gaussiana

- **Contexto:** la resistencia es la principal fuente de incertidumbre en taludes de suelo, y c'
  y φ' suelen tener correlación negativa al ajustarse juntas desde ensayos.
- **Decisión:** admitir distribuciones normal y lognormal, imponer la correlación en el espacio
  normal estándar con Cholesky y mantener γ determinista.
- **Consecuencias:** el muestreo es simple y verificable. La incertidumbre del peso unitario y del
  nivel freático no se modela; el nivel freático se estudia por escenarios.

## ADR-005: Probabilidad de falla con el mínimo entre superficies candidatas

- **Contexto:** calcular la PF sólo sobre la superficie crítica determinista subestima la
  probabilidad, porque otra superficie puede ser más crítica para ciertas combinaciones de c' y
  φ'. Buscar el círculo crítico en cada una de las 20 000 muestras es demasiado lento.
- **Decisión:** evaluar cada muestra en un conjunto de unas 20 superficies: la crítica, las
  críticas de búsquedas en μ ± σ y círculos de bajo FS de la malla. Se toma el mínimo y también se
  informa el resultado en la superficie fija.
- **Consecuencias:** en el caso de estudio la PF pasa de 10.53 % (superficie fija) a 11.24 %. Frente
  a una búsqueda completa, el error en casos extremos es ≤ 0.17 %.

## ADR-006: Cálculo vectorizado con numpy y dependencias mínimas

- **Contexto:** la búsqueda evalúa miles de círculos y Monte Carlo miles de muestras por
  superficie.
- **Decisión:** escribir Fellenius y Bishop sobre arreglos de forma (..., n_dovelas), combinables
  por broadcasting, y depender sólo de numpy y matplotlib.
- **Consecuencias:** el caso completo corre en unos 5 s. Una implementación independiente dovela
  por dovela verifica que la vectorización no altera el resultado.
