# Guía para agentes y revisores

- Instalar dependencias: `pip install -r requirements.txt` (numpy y matplotlib).
- Pruebas: `python -m unittest` o `pytest`. Deben pasar todas antes de cada commit.
- Analizar un caso: `python -m estabilidad_taludes analizar casos/talud_cobertura.json`.
- Ángulo de diseño: `python -m estabilidad_taludes diseno casos/talud_cobertura.json`.

## Convenciones

- Español en código, mensajes, figuras y documentación.
- Unidades: m, kPa, kN/m³ y grados. El talud desciende hacia la izquierda (pie a la izquierda,
  corona a la derecha).
- Las funciones de cálculo están vectorizadas: los arreglos de dovelas tienen forma (..., n) y
  los parámetros se combinan por broadcasting.
- Requisitos en `docs/requisitos.md` (con su prueba asociada) y decisiones en `docs/decisiones.md`.
  Si un cambio altera un valor de validación, actualizar `docs/validacion.md`.
- `resultados/` se regenera con los comandos de arriba; volver a generarlo si cambia el cálculo.
- No subir credenciales ni datos privados.
