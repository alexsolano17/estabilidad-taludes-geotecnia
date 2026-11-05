"""Validación contra resultados publicados (ver docs/validacion.md).

- ACADS 1(a), Giam y Donald (1989): FS de referencia 1.00. Manual de verificación de
  Rocscience: Slide2 Bishop 0.987; verificación de OptumG2: SVSLOPE Bishop 0.989 y
  método ordinario 0.947 (Slide) / 0.945 (SVSLOPE).
- Corte vertical en suelo con φ = 0 (Taylor, 1937): altura crítica Hc = 3.83·c/γ.
"""

import unittest
from pathlib import Path

from estabilidad_taludes import Material, buscar_circulo_critico, talud_simple
from estabilidad_taludes.casos import cargar_caso

CASOS = Path(__file__).resolve().parents[1] / "casos"


class TestACADS1a(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        caso = cargar_caso(CASOS / "acads_1a.json")
        cls.talud, cls.material = caso.talud(), caso.material

    def test_bishop(self):
        fs = buscar_circulo_critico(self.talud, self.material, "bishop").critica.fs
        self.assertAlmostEqual(fs, 0.988, delta=0.010)  # Slide2 0.987 / SVSLOPE 0.989
        self.assertAlmostEqual(fs, 1.00, delta=0.025)  # FS de referencia del estudio ACADS

    def test_fellenius(self):
        fs = buscar_circulo_critico(self.talud, self.material, "fellenius").critica.fs
        self.assertAlmostEqual(fs, 0.946, delta=0.010)  # Slide 0.947 / SVSLOPE 0.945


class TestTaylorCorteVertical(unittest.TestCase):
    def test_numero_de_estabilidad(self):
        altura, gamma, c = 10.0, 20.0, 50.0
        talud = talud_simple(altura, 90.0)
        fs = buscar_circulo_critico(talud, Material(c=c, phi=0.0, gamma=gamma)).critica.fs
        numero_estabilidad = fs * gamma * altura / c  # Hc = Ns·c/γ
        self.assertAlmostEqual(numero_estabilidad, 3.83, delta=0.04)


if __name__ == "__main__":
    unittest.main()
