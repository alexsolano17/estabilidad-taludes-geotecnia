"""Pruebas de extremo a extremo de la línea de comandos."""

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from estabilidad_taludes.cli import main

CASOS = Path(__file__).resolve().parents[1] / "casos"


class TestLineaDeComandos(unittest.TestCase):
    def _ejecutar(self, *args):
        salida = io.StringIO()
        with contextlib.redirect_stdout(salida):
            codigo = main([str(a) for a in args])
        self.assertEqual(codigo, 0)
        return salida.getvalue()

    def test_analizar_caso_determinista(self):
        with tempfile.TemporaryDirectory() as tmp:
            self._ejecutar("analizar", CASOS / "acads_1a.json", "--salida", tmp)
            datos = json.loads((Path(tmp) / "resumen.json").read_text(encoding="utf-8"))
            self.assertAlmostEqual(datos["determinista"]["fs_bishop"], 0.985, delta=0.005)
            self.assertTrue((Path(tmp) / "talud_critico.png").stat().st_size > 10_000)
            self.assertNotIn("probabilistico", datos)

    def test_analizar_caso_probabilistico(self):
        with tempfile.TemporaryDirectory() as tmp:
            texto = self._ejecutar("analizar", CASOS / "talud_cobertura.json", "--salida", tmp, "--simulaciones", "2000")
            self.assertIn("P(FS < 1)", texto)
            carpeta = Path(tmp)
            for nombre in ("talud_critico.png", "histograma_fs.png", "convergencia_pf.png", "resumen.md"):
                self.assertTrue((carpeta / nombre).exists(), nombre)
            filas = (carpeta / "muestras.csv").read_text(encoding="utf-8").strip().splitlines()
            self.assertEqual(len(filas), 2001)
            datos = json.loads((carpeta / "resumen.json").read_text(encoding="utf-8"))
            self.assertEqual(datos["probabilistico"]["simulaciones"], 2000)
            self.assertIn("cumple", datos)

    def test_diseno(self):
        with tempfile.TemporaryDirectory() as tmp:
            texto = self._ejecutar("diseno", CASOS / "talud_cobertura.json", "--desde", "40", "--hasta", "45",
                                   "--paso", "5", "--salida", tmp, "--simulaciones", "1000")
            self.assertIn("Ángulo recomendado: 40°", texto)
            self.assertTrue((Path(tmp) / "diseno_angulo.png").exists())
            self.assertEqual(len((Path(tmp) / "diseno.csv").read_text(encoding="utf-8").splitlines()), 3)


if __name__ == "__main__":
    unittest.main()
