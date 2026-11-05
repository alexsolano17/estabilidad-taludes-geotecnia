# Requisitos

Especificación breve del programa, con el criterio de aceptación de cada requisito y la prueba
automática que lo verifica.

## Objetivo

Calcular el factor de seguridad y la probabilidad de falla de un talud en suelo o roca muy
fracturada, con superficies de falla circulares, para apoyar la elección del ángulo de diseño
según criterios de aceptación de minería a cielo abierto.

## Requisitos funcionales

| ID | Requisito | Criterio de aceptación | Pruebas |
|---|---|---|---|
| RF-01 | Definir el talud como talud simple (altura y ángulo) o como perfil por puntos, con estrato firme opcional. | Ambos formatos se leen desde JSON y se analizan. | `test_cli.py` |
| RF-02 | Considerar un nivel freático con presión de poros hidrostática. | El FS baja al subir el nivel; un nivel bajo la superficie no cambia el FS. | `test_subir_el_nivel_freatico_reduce_el_fs`, `test_nivel_freatico_profundo_no_cambia_el_fs` |
| RF-03 | Calcular el FS por Bishop simplificado y por Fellenius. | Igualdad con una implementación independiente; límite de talud infinito; Bishop = Fellenius con φ = 0. | `test_metodos.py` |
| RF-04 | Encontrar el círculo crítico sin que el usuario dé la superficie. | ACADS 1(a) dentro de ±0.01 de Slide2/SVSLOPE; corte vertical de Taylor con error ≤ 1 %. | `test_validacion.py`, `test_busqueda_estable_ante_una_malla_mas_fina` |
| RF-05 | Tratar c' y φ' como variables aleatorias normales o lognormales, con correlación. | Medias y desviaciones de las muestras con error < 1 %; correlación ±0.01; correlación no válida rechazada. | `test_momentos_y_correlacion`, `test_correlacion_no_valida` |
| RF-06 | Estimar P(FS < 1) por Monte Carlo con intervalo de confianza e índices β. | Coincide con la solución exacta para φ = 0 dentro de 4 errores estándar. | `test_pf_analitica_con_phi_cero`, `test_indice_lognormal`, `test_intervalo_de_wilson` |
| RF-07 | Considerar que la superficie crítica cambia entre muestras. | El mínimo entre superficies candidatas difiere < 1 % de una búsqueda completa. | `test_superficies_candidatas_frente_a_busqueda_completa` |
| RF-08 | Controlar Monte Carlo con un método independiente (Rosenblueth). | Media con diferencia < 1 % y desviación < 5 %. | `test_rosenblueth_coincide_con_monte_carlo` |
| RF-09 | Comparar con criterios de aceptación (FS mínimo y PF máxima) y recomendar un ángulo de diseño. | El barrido informa el cumplimiento por ángulo y el ángulo recomendado. | `test_diseno` |
| RF-10 | Generar figuras y resúmenes (JSON, Markdown y CSV de muestras). | Se crean todos los archivos de salida con el contenido esperado. | `test_analizar_caso_determinista`, `test_analizar_caso_probabilistico` |

## Requisitos no funcionales

| ID | Requisito | Cómo se cumple |
|---|---|---|
| RNF-01 | Reproducibilidad | Semilla fija en cada caso; las salidas no incluyen fechas ni tiempos (`test_misma_semilla_mismos_resultados`). |
| RNF-02 | Rendimiento | Análisis completo del caso de estudio en unos 5 s; barrido de 11 ángulos en menos de 1 min; pruebas en unos 10 s. |
| RNF-03 | Dependencias mínimas | Sólo numpy y matplotlib. Las pruebas usan `unittest` de la biblioteca estándar; pytest es opcional. |
| RNF-04 | Idioma | Código, mensajes, figuras y documentación en español. |
| RNF-05 | Trazabilidad | Cada requisito tiene una prueba y cada validación cita su fuente (`docs/validacion.md`). |

## Fuera del alcance

Superficies no circulares, métodos rigurosos (Spencer, Morgenstern-Price), varios materiales,
sismo o cargas externas, grietas de tracción, análisis 3D, flujo transitorio y fallas
estructuralmente controladas en roca.
