"""Pruebas del análisis probabilístico."""

import math
import unittest
from pathlib import Path
from statistics import NormalDist

import numpy as np

from estabilidad_taludes import Material, buscar_circulo_critico, talud_simple
from estabilidad_taludes.casos import cargar_caso
from estabilidad_taludes.probabilistico import (
    Variable,
    estadisticas,
    fs_en_superficies,
    intervalo_wilson,
    matriz_correlacion,
    monte_carlo,
    muestrear,
    rosenblueth,
    superficies_para_monte_carlo,
)

CASOS = Path(__file__).resolve().parents[1] / "casos"
NORMAL = NormalDist()


class TestMuestreo(unittest.TestCase):
    def test_momentos_y_correlacion(self):
        variables = [Variable("c", "lognormal", 20.0, 6.0), Variable("phi", "normal", 30.0, 3.0)]
        corr = matriz_correlacion(variables, -0.5)
        muestras, recortadas = muestrear(variables, corr, 400_000, semilla=1)
        c, phi = muestras["c"], muestras["phi"]
        self.assertEqual(recortadas, 0)
        self.assertAlmostEqual(c.mean() / 20.0, 1.0, delta=0.005)
        self.assertAlmostEqual(c.std() / 6.0, 1.0, delta=0.01)
        self.assertAlmostEqual(phi.mean() / 30.0, 1.0, delta=0.002)
        self.assertAlmostEqual(phi.std() / 3.0, 1.0, delta=0.01)
        # La correlación se impone entre ln(c) y φ (espacio normal estándar).
        self.assertAlmostEqual(np.corrcoef(np.log(c), phi)[0, 1], -0.5, delta=0.01)

    def test_misma_semilla_mismos_resultados(self):
        variables = [Variable("c", "normal", 20.0, 4.0)]
        a, _ = muestrear(variables, np.eye(1), 1000, semilla=7)
        b, _ = muestrear(variables, np.eye(1), 1000, semilla=7)
        np.testing.assert_array_equal(a["c"], b["c"])

    def test_correlacion_no_valida(self):
        with self.assertRaises(np.linalg.LinAlgError):
            muestrear([Variable("c", "normal", 20, 4), Variable("phi", "normal", 30, 3)],
                      np.array([[1.0, 1.5], [1.5, 1.0]]), 10, semilla=1)


class TestProbabilidadDeFalla(unittest.TestCase):
    def test_pf_analitica_con_phi_cero(self):
        """Con φ = 0 el FS de una superficie fija es proporcional a c: FS = k·c.

        Si c es normal(μ, σ), entonces P(FS < 1) = Φ((1/k − μ)/σ) de forma exacta.
        """
        talud = talud_simple(10.0, 60.0)
        material = Material(c=40.0, phi=0.0, gamma=20.0)
        critica = buscar_circulo_critico(talud, material).critica
        k = critica.fs / material.c
        variables = [Variable("c", "normal", 40.0, 8.0)]
        mc = monte_carlo(talud, material, [critica], variables, np.eye(1), 200_000, semilla=11)
        est = estadisticas(mc.fs_fija)
        exacta = NORMAL.cdf((1.0 / k - 40.0) / 8.0)
        self.assertGreater(exacta, 0.05)
        self.assertLess(abs(est.pf - exacta), 4.0 * est.error_pf)

    def test_indice_lognormal(self):
        """Para un FS lognormal, Φ(−β_lognormal) coincide con la fracción FS < 1."""
        rng = np.random.default_rng(3)
        media, cov = 1.25, 0.15
        s = math.sqrt(math.log1p(cov**2))
        fs = np.exp(math.log(media) - 0.5 * s * s + s * rng.standard_normal(400_000))
        est = estadisticas(fs)
        self.assertLess(abs(est.pf - est.pf_lognormal), 4.0 * est.error_pf)

    def test_intervalo_de_wilson(self):
        bajo, alto = intervalo_wilson(0, 1000)
        self.assertEqual(bajo, 0.0)
        self.assertAlmostEqual(alto, 3.8415 / 1003.8415, places=4)
        bajo, alto = intervalo_wilson(50, 1000)
        self.assertLess(bajo, 0.05)
        self.assertGreater(alto, 0.05)


class TestCasoDeEstudio(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.caso = cargar_caso(CASOS / "talud_cobertura.json")
        cls.talud = cls.caso.talud()
        cls.busqueda = buscar_circulo_critico(cls.talud, cls.caso.material)
        cls.superficies = superficies_para_monte_carlo(
            cls.talud, cls.caso.material, cls.caso.variables, cls.busqueda, cls.caso.busqueda
        )

    def test_rosenblueth_coincide_con_monte_carlo(self):
        caso, crit = self.caso, self.busqueda.critica
        mc = monte_carlo(self.talud, caso.material, [crit], caso.variables, caso.correlacion, 40_000, semilla=5)
        est = estadisticas(mc.fs_fija)

        def funcion(c, phi):
            return fs_en_superficies(self.talud, caso.material, [crit], c, phi)[0][:, 0]

        ros = rosenblueth(funcion, caso.variables, caso.correlacion)
        self.assertAlmostEqual(ros.media / est.media, 1.0, delta=0.01)
        self.assertAlmostEqual(ros.desviacion / est.desviacion, 1.0, delta=0.05)

    def test_superficies_candidatas_frente_a_busqueda_completa(self):
        """El mínimo entre candidatas debe acercarse al de una búsqueda completa por muestra."""
        caso = self.caso
        pares = [(10.0, 33.0), (14.0, 27.0), (20.0, 30.0), (28.0, 26.0), (34.0, 33.0), (12.0, 36.0)]
        c = np.array([p[0] for p in pares])
        phi = np.array([p[1] for p in pares])
        fs, ok = fs_en_superficies(self.talud, caso.material, self.superficies, c, phi)
        fs_candidatas = np.where(ok, fs, np.inf).min(axis=1)
        for i, (ci, phii) in enumerate(pares):
            completa = buscar_circulo_critico(self.talud, caso.material.con(c=ci, phi=phii)).critica.fs
            self.assertGreaterEqual(fs_candidatas[i], completa - 1e-6)
            self.assertLess(fs_candidatas[i] / completa - 1.0, 0.01, msg=f"c={ci}, φ={phii}")


if __name__ == "__main__":
    unittest.main()
