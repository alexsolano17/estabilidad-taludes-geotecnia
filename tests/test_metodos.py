"""Pruebas del cálculo determinista: dovelas, Fellenius, Bishop y búsqueda."""

import itertools
import math
import unittest

import numpy as np

from estabilidad_taludes import (
    ConfigBusqueda,
    Material,
    Talud,
    buscar_circulo_critico,
    circulo_por_cuerda,
    dovelas,
    evaluar_circulos,
    fs_bishop,
    fs_fellenius,
    talud_simple,
)
from estabilidad_taludes.metodos import factor_seguridad

AGUA = {"altura_en_pie": 0.0, "profundidad_bajo_corona": 5.0, "distancia_tras_corona": 10.0}


def bishop_dovela_por_dovela(talud, material, xc, yc, radio, x1, x2, n):
    """Implementación independiente en Python puro, sin numpy ni vectorización."""
    b = (x2 - x1) / n
    filas = []
    for i in range(n):
        x = x1 + (i + 0.5) * b
        y_base = yc - math.sqrt(radio**2 - (x - xc) ** 2)
        y_sup = float(talud.y_terreno(x))
        h = max(y_sup - y_base, 0.0)
        hw = 0.0
        if talud.nivel_freatico is not None:
            y_agua = min(float(talud.y_freatico(x)), y_sup)
            hw = min(max(y_agua - y_base, 0.0), h)
        peso = b * (material.gamma * (h - hw) + material.gamma_saturado * hw)
        sen = (x - xc) / radio
        filas.append((peso, talud.gamma_agua * hw, sen, math.sqrt(1.0 - sen * sen)))
    t = material.tan_phi
    motor = sum(w * s for w, _, s, _ in filas)
    fs = 1.0
    for _ in range(500):
        nuevo = sum((material.c * b + max(w - u * b, 0.0) * t) / (co + s * t / fs) for w, u, s, co in filas)
        nuevo /= motor
        if abs(nuevo - fs) < 1e-13:
            break
        fs = nuevo
    return nuevo


class TestMetodos(unittest.TestCase):
    def setUp(self):
        self.talud = talud_simple(20.0, 45.0, 40.0, 60.0, nivel_freatico=AGUA)
        self.material = Material(c=20.0, phi=30.0, gamma=19.0, gamma_sat=20.0)
        self.circulo = (42.0, 66.0, 32.0)  # x_salida, x_entrada, psi

    def _dovelas(self, talud=None, material=None, n=100):
        talud, material = talud or self.talud, material or self.material
        x1, x2, psi = self.circulo
        xc, yc, r = circulo_por_cuerda(talud, x1, x2, psi)
        return dovelas(talud, material, xc, yc, r, x1, x2, n=n), (float(xc), float(yc), float(r))

    def test_bishop_coincide_con_implementacion_independiente(self):
        dov, (xc, yc, r) = self._dovelas()
        fs, _, convergio = fs_bishop(dov, self.material.c, self.material.tan_phi)
        x1, x2, _ = self.circulo
        esperado = bishop_dovela_por_dovela(self.talud, self.material, xc, yc, r, x1, x2, 100)
        self.assertTrue(bool(convergio))
        self.assertAlmostEqual(float(fs), esperado, places=9)

    def test_circulo_pasa_por_los_puntos_de_entrada_y_salida(self):
        x1, x2, psi = self.circulo
        xc, yc, r = circulo_por_cuerda(self.talud, x1, x2, psi)
        for x in (x1, x2):
            y = float(self.talud.y_terreno(x))
            self.assertAlmostEqual(math.hypot(x - xc, y - yc), float(r), places=9)

    def test_phi_cero_bishop_igual_a_fellenius(self):
        material = Material(c=40.0, phi=0.0, gamma=19.0)
        dov, _ = self._dovelas(material=material)
        fs_b, _, _ = fs_bishop(dov, material.c, 0.0)
        fs_f = fs_fellenius(dov, material.c, 0.0)
        self.assertAlmostEqual(float(fs_b), float(fs_f), places=10)

    def test_fellenius_menor_que_bishop_con_friccion(self):
        dov, _ = self._dovelas()
        fs_b, _, _ = fs_bishop(dov, self.material.c, self.material.tan_phi)
        fs_f = fs_fellenius(dov, self.material.c, self.material.tan_phi)
        self.assertLess(float(fs_f), float(fs_b))

    def test_limite_talud_infinito_sin_cohesion(self):
        """Arco casi plano paralelo a la cara con c' = 0: FS → tanφ'/tanβ."""
        talud = talud_simple(20.0, 30.0)
        material = Material(c=0.0, phi=25.0, gamma=19.0)
        x_pie, x_corona, _ = talud.cara()
        x1, x2 = x_pie + 5.0, x_corona - 5.0
        esperado = math.tan(math.radians(25.0)) / math.tan(math.radians(30.0))
        for metodo in ("bishop", "fellenius"):
            fs = float(evaluar_circulos(talud, material, x1, x2, 0.5, metodo, n=400))
            self.assertAlmostEqual(fs / esperado, 1.0, delta=0.002, msg=metodo)

    def test_nivel_freatico_profundo_no_cambia_el_fs(self):
        seco = talud_simple(20.0, 45.0, 40.0, 60.0)
        profundo = Talud(seco.superficie, [(0.0, -60.0), (200.0, -60.0)])
        x1, x2, psi = self.circulo
        fs_seco = evaluar_circulos(seco, self.material, x1, x2, psi)
        fs_prof = evaluar_circulos(profundo, self.material, x1, x2, psi)
        self.assertEqual(float(fs_seco), float(fs_prof))

    def test_subir_el_nivel_freatico_reduce_el_fs(self):
        valores = []
        for profundidad in (15.0, 10.0, 5.0, 2.0):
            agua = {"altura_en_pie": 0.0, "profundidad_bajo_corona": profundidad, "distancia_tras_corona": 10.0}
            talud = talud_simple(20.0, 45.0, 40.0, 60.0, nivel_freatico=agua)
            valores.append(buscar_circulo_critico(talud, self.material).critica.fs)
        self.assertTrue(all(a > b for a, b in itertools.pairwise(valores)), valores)

    def test_convergencia_con_el_numero_de_dovelas(self):
        fs_100 = float(factor_seguridad(self._dovelas(n=100)[0], self.material.c, self.material.tan_phi))
        fs_1600 = float(factor_seguridad(self._dovelas(n=1600)[0], self.material.c, self.material.tan_phi))
        self.assertLess(abs(fs_100 / fs_1600 - 1.0), 1e-3)

    def test_superficies_no_admisibles(self):
        # Entrada por encima del centro: el arco dejaría de ser una función de x.
        self.assertEqual(float(evaluar_circulos(self.talud, self.material, 42.0, 66.0, 89.0)), math.inf)
        # Estrato firme por encima del fondo del círculo.
        firme = Talud(self.talud.superficie, self.talud.nivel_freatico, y_estrato_firme=5.0)
        self.assertEqual(float(evaluar_circulos(firme, self.material, 42.0, 66.0, 32.0)), math.inf)
        # Arco que corta la cara: sale por debajo de la cara y vuelve a entrar.
        self.assertEqual(float(evaluar_circulos(self.talud, self.material, 30.0, 75.0, 3.0)), math.inf)

    def test_vectorizacion_por_muestras(self):
        dov, _ = self._dovelas()
        c = np.array([10.0, 20.0, 35.0])
        phi = np.array([25.0, 30.0, 33.0])
        t = np.tan(np.radians(phi))
        fs_vec, _, _ = fs_bishop(dov, c[:, None], t[:, None])
        for i in range(3):
            fs_i, _, _ = fs_bishop(dov, c[i], t[i])
            # Igualdad hasta la tolerancia de la iteración (1e-9 relativa).
            self.assertAlmostEqual(float(fs_vec[i, 0]), float(fs_i), places=8)

    def test_busqueda_estable_ante_una_malla_mas_fina(self):
        base = buscar_circulo_critico(self.talud, self.material).critica.fs
        fina = ConfigBusqueda(n_salida=40, n_entrada=40, n_psi=30, n_refinar=8)
        fs_fina = buscar_circulo_critico(self.talud, self.material, config=fina).critica.fs
        self.assertLess(abs(base - fs_fina), 2e-3)


if __name__ == "__main__":
    unittest.main()
