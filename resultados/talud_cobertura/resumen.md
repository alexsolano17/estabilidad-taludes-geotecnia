# Resultados: Talud en material de cobertura de un tajo abierto

Caso de estudio ilustrativo: talud de corte de 20 m en suelo residual / roca muy meteorizada con nivel freático. Los parámetros son supuestos razonables, no ensayos: para un diseño real deben reemplazarse por resultados de laboratorio y campo.

## Datos

| Parámetro | Valor |
|---|---|
| Altura | 20 m |
| Ángulo de la cara | 45° |
| c' | 20 kPa |
| φ' | 30° |
| γ / γsat | 19 / 20 kN/m³ |
| Nivel freático | sí |
| c' aleatoria | lognormal, media 20 kPa, COV 30% |
| φ' aleatoria | normal, media 30°, COV 10% |
| Correlación ρ(c', φ') | -0.5 |

## Análisis determinista

| Resultado | Valor |
|---|---|
| FS Bishop simplificado | **1.120** |
| FS Fellenius (ordinario) | 1.010 |
| FS Bishop sin nivel freático | 1.226 |
| Círculo crítico | centro (40.01, 25.03) m, radio 25.03 m |
| Salida / entrada | x = 40.00 m / x = 64.53 m |
| Círculos válidos en la malla | 6 697 de 12 160 |

## Análisis probabilístico

20 000 simulaciones de Monte Carlo con semilla 2026; 19 superficies evaluadas por muestra.

| Resultado | Mínimo entre superficies | Superficie crítica fija |
|---|---|---|
| Media del FS | 1.119 | 1.123 |
| Desviación del FS | 0.102 | 0.103 |
| Percentil 5 % del FS | 0.963 | 0.966 |
| **P(FS < 1)** | **11.24 %** | 10.53 % |
| IC 95 % de P(FS < 1) | 10.81 % – 11.69 % | 10.11 % – 10.96 % |
| β (FS normal) | 1.17 | 1.19 |
| β (FS lognormal) | 1.19 | 1.22 |

Control con Rosenblueth (superficie fija, 4 puntos): media 1.123, desviación 0.102.

## Criterios de aceptación

Talud interrampa con consecuencia de falla media (Read y Stacey, 2009)

| Criterio | Límite | Resultado | Estado |
|---|---|---|---|
| FS mínimo | ≥ 1.2 | 1.120 | no cumple |
| P(FS < 1) máxima | ≤ 20.00 % | 11.24 % | cumple |
