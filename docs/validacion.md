# Validación

Cada comparación se comprueba en una prueba automática (columna «Prueba»). Si un cambio en el
código altera alguno de estos valores, la prueba falla. Los valores obtenidos corresponden a la
versión 1.0.0 con la configuración de búsqueda por defecto.

## 1. Casos de referencia publicados

### ACADS 1(a): talud homogéneo sin agua

Problema de verificación del estudio ACADS (Giam y Donald, 1989), usado por los programas
comerciales para verificar sus resultados. Datos según el manual de verificación de Rocscience
(problema 1): pie en (30, 25), corona en (50, 35), talud 1V:2H de 10 m, base del modelo en
y = 20 m, c' = 3 kPa, φ' = 19.6° y γ = 20 kN/m³. Caso: [`casos/acads_1a.json`](../casos/acads_1a.json).

| Método | Referencia | Este programa | Diferencia |
|---|---|---|---|
| Bishop simplificado | 1.00 (FS de referencia del estudio) | 0.985 | −1.5 % |
| Bishop simplificado | Slide2 0.987 · SVSLOPE 0.989 | 0.985 | −0.2 % / −0.4 % |
| Ordinario (Fellenius) | Slide 0.947 · SVSLOPE 0.945 | 0.942 | −0.5 % / −0.3 % |

El valor 1.00 del estudio ACADS es un promedio de referencia: los programas comerciales
obtienen 0.986–0.990 con distintos métodos de equilibrio límite. Prueba: `tests/test_validacion.py::TestACADS1a`.

Fuentes consultadas:
- [Rocscience, RS2 Slope Stability Verification, parte 1](https://www.rocscience.com/assets/verification-and-theory/RS2/RS2_SlopeStabilityVerification_Part1.pdf) (geometría, FS de referencia 1.00 y resultados de Slide2).
- [OptumG2, Slope Stability Verification](https://optumce.com/Download/Verifications/g2/G2_Slope_Stability_Verification.pdf) (resultados de Slide y SVSLOPE con Bishop y el método ordinario).
- [Registro de Giam y Donald (1989) en TRID](https://trid.trb.org/view/313937).

### Corte vertical en suelo puramente cohesivo (φ = 0)

Para un corte vertical con φ = 0 y superficie circular, la altura crítica es Hc = 3.83·c/γ
(número de estabilidad de Taylor, 1937; ver por ejemplo [este resumen](https://engineeringhulk.com/civil/geotechnical/slope-stability/)).
Con H = 10 m, γ = 20 kN/m³ y c = 50 kPa, el programa obtiene FS·γ·H/c = **3.831**.
Prueba: `tests/test_validacion.py::TestTaylorCorteVertical`.

## 2. Soluciones analíticas

| Comprobación | Valor esperado | Obtenido | Prueba |
|---|---|---|---|
| Talud infinito con c' = 0 (β = 30°, φ' = 25°): un arco casi plano paralelo a la cara tiende a FS = tanφ'/tanβ | 0.80767 | Bishop 0.80770 · Fellenius 0.80768 | `test_limite_talud_infinito_sin_cohesion` |
| Con φ = 0, Bishop y Fellenius coinciden | igualdad | diferencia < 1e-10 | `test_phi_cero_bishop_igual_a_fellenius` |
| Probabilidad de falla con φ = 0: el FS de una superficie fija es k·c, así que P(FS < 1) = Φ((1/k − μ)/σ) con c normal(40, 8) kPa | 40.68 % | 40.75 % (200 000 simulaciones, error estándar 0.11 %) | `test_pf_analitica_con_phi_cero` |
| Si el FS es lognormal, Φ(−β lognormal) iguala la fracción de muestras con FS < 1 | igualdad | dentro de 4 errores estándar | `test_indice_lognormal` |

## 3. Controles internos

| Comprobación | Resultado | Prueba |
|---|---|---|
| Bishop vectorizado frente a una implementación independiente dovela por dovela | igualdad a 9 decimales | `test_bishop_coincide_con_implementacion_independiente` |
| Convergencia con el número de dovelas (círculo de prueba del caso de estudio) | 100 dovelas: 1.25253 · 1600 dovelas: 1.25255 | `test_convergencia_con_el_numero_de_dovelas` |
| Búsqueda con una malla más fina (40 × 40 × 30) | mismo mínimo (diferencia < 0.002) | `test_busqueda_estable_ante_una_malla_mas_fina` |
| El FS baja al subir el nivel freático; un nivel bajo la superficie de falla no lo cambia | se cumple | `test_subir_el_nivel_freatico_reduce_el_fs`, `test_nivel_freatico_profundo_no_cambia_el_fs` |
| Rosenblueth (4 puntos) frente a Monte Carlo (40 000 muestras) en la superficie crítica del caso de estudio | media 1.1231 frente a 1.1232; desviación 0.1020 frente a 0.1034 | `test_rosenblueth_coincide_con_monte_carlo` |
| Mínimo entre superficies candidatas frente a una búsqueda completa, en 6 combinaciones extremas de c' y φ' | diferencia máxima 0.17 % (nunca menor que la búsqueda completa) | `test_superficies_candidatas_frente_a_busqueda_completa` |
| Momentos y correlación de las muestras | media y desviación con error < 1 %; ρ = −0.50 ± 0.01 | `test_momentos_y_correlacion` |

## 4. Criterios de aceptación

Los límites del caso de estudio provienen de la tabla de criterios de aceptación de
Read y Stacey (2009), reproducida en [esta presentación de SRK (Hormazábal, 2018)](https://www.kz.srk.com/sites/default/files/file/EHormazabal_Some_geotechnical_considerations_for_probabilistic_analysis_2018.pdf):

| Escala del talud | Consecuencia | FS estático mínimo | P(FS ≤ 1) máxima |
|---|---|---|---|
| Banco | baja a alta | 1.1 | 25–50 % |
| Interrampa | baja | 1.15–1.2 | 25 % |
| Interrampa | media | 1.2 | 20 % |
| Interrampa | alta | 1.2–1.3 | 10 % |
| Global | baja | 1.2–1.3 | 15–20 % |
| Global | media | 1.3 | 10 % |
| Global | alta | 1.3–1.5 | 5 % |

El caso de estudio usa la fila «interrampa, consecuencia media». Para otro talud basta con
cambiar `criterios` en su archivo JSON.
