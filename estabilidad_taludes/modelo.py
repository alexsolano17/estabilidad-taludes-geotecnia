"""Modelo del talud: perfil del terreno, nivel freático y material (deformación plana).

Convención de ejes: x horizontal y y vertical hacia arriba. El talud desciende hacia
la izquierda: el pie queda a la izquierda, la corona a la derecha y la masa deslizante
se mueve en el sentido -x.

Unidades: m, kPa, kN/m³ y grados sexagesimales.
"""

from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass

import numpy as np

GAMMA_AGUA = 9.81  # peso unitario del agua [kN/m³]


@dataclass(frozen=True)
class Material:
    """Suelo homogéneo con resistencia de Mohr-Coulomb en esfuerzos efectivos."""

    c: float  # cohesión efectiva c' [kPa]
    phi: float  # ángulo de fricción efectivo φ' [°]
    gamma: float  # peso unitario sobre el nivel freático [kN/m³]
    gamma_sat: float | None = None  # peso unitario bajo el nivel freático; None = gamma

    def __post_init__(self) -> None:
        if self.c < 0:
            raise ValueError("La cohesión c' no puede ser negativa.")
        if not 0.0 <= self.phi < 90.0:
            raise ValueError("El ángulo de fricción φ' debe estar entre 0° y 90°.")
        if self.gamma <= 0 or (self.gamma_sat is not None and self.gamma_sat <= 0):
            raise ValueError("Los pesos unitarios deben ser positivos.")

    @property
    def gamma_saturado(self) -> float:
        return self.gamma if self.gamma_sat is None else self.gamma_sat

    @property
    def tan_phi(self) -> float:
        return math.tan(math.radians(self.phi))

    def con(self, **cambios: float) -> Material:
        """Copia del material con algunos parámetros cambiados."""
        return dataclasses.replace(self, **cambios)


def _polilinea(puntos, nombre: str) -> np.ndarray:
    arr = np.asarray(puntos, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2 or len(arr) < 2:
        raise ValueError(f"{nombre}: se esperan al menos dos puntos (x, y).")
    if np.any(np.diff(arr[:, 0]) <= 0):
        raise ValueError(f"{nombre}: las abscisas deben ser estrictamente crecientes.")
    return arr


@dataclass(eq=False)
class Talud:
    """Perfil del terreno (polilínea), nivel freático opcional y estrato firme opcional.

    superficie: puntos (x, y) del terreno con x creciente.
    nivel_freatico: puntos (x, y) del nivel freático, o None si el talud está seco.
    y_estrato_firme: cota bajo la cual no pueden pasar las superficies de falla.
    """

    superficie: np.ndarray
    nivel_freatico: np.ndarray | None = None
    y_estrato_firme: float | None = None
    gamma_agua: float = GAMMA_AGUA

    def __post_init__(self) -> None:
        self.superficie = _polilinea(self.superficie, "Superficie del terreno")
        if self.nivel_freatico is not None:
            self.nivel_freatico = _polilinea(self.nivel_freatico, "Nivel freático")

    @property
    def x_min(self) -> float:
        return float(self.superficie[0, 0])

    @property
    def x_max(self) -> float:
        return float(self.superficie[-1, 0])

    def y_terreno(self, x):
        """Cota del terreno en x (interpolación lineal; horizontal fuera del perfil)."""
        return np.interp(x, self.superficie[:, 0], self.superficie[:, 1])

    def y_freatico(self, x):
        """Cota del nivel freático en x sin recortar al terreno; None si no hay agua."""
        if self.nivel_freatico is None:
            return None
        return np.interp(x, self.nivel_freatico[:, 0], self.nivel_freatico[:, 1])

    def cara(self) -> tuple[float, float, float]:
        """(x del pie, x de la corona, altura) a partir de los tramos ascendentes del perfil."""
        dy = np.diff(self.superficie[:, 1])
        sube = np.flatnonzero(dy > 1e-9)
        if sube.size == 0:
            raise ValueError("El perfil no tiene tramos ascendentes: no hay talud.")
        x_pie, y_pie = self.superficie[sube[0]]
        x_corona, y_corona = self.superficie[sube[-1] + 1]
        return float(x_pie), float(x_corona), float(y_corona - y_pie)


def nivel_freatico_simple(
    x_pie: float,
    x_corona: float,
    altura: float,
    x_fin: float,
    altura_en_pie: float = 0.0,
    profundidad_bajo_corona: float = 0.0,
    distancia_tras_corona: float = 0.0,
) -> np.ndarray:
    """Nivel freático lineal para un talud simple.

    Aflora a `altura_en_pie` sobre la plataforma del pie, sube en línea recta hasta la cota
    `altura - profundidad_bajo_corona` a `distancia_tras_corona` detrás de la corona y luego
    sigue horizontal.
    """
    x_int = x_corona + distancia_tras_corona
    y_int = altura - profundidad_bajo_corona
    puntos = [(0.0, altura_en_pie), (x_pie, altura_en_pie), (x_int, y_int)]
    puntos.append((max(x_fin, x_int + 1.0), y_int))
    return np.array(puntos)


def talud_simple(
    altura: float,
    angulo: float,
    largo_pie: float | None = None,
    largo_corona: float | None = None,
    nivel_freatico: dict | None = None,
    y_estrato_firme: float | None = None,
    gamma_agua: float = GAMMA_AGUA,
) -> Talud:
    """Talud simple: plataforma al pie (y = 0), cara a `angulo` grados y corona (y = altura).

    nivel_freatico: argumentos de `nivel_freatico_simple` (altura_en_pie,
    profundidad_bajo_corona, distancia_tras_corona) o None para un talud seco.
    """
    if altura <= 0:
        raise ValueError("La altura del talud debe ser positiva.")
    if not 0.0 < angulo <= 90.0:
        raise ValueError("El ángulo del talud debe estar entre 0° y 90°.")
    largo_pie = 2.0 * altura if largo_pie is None else float(largo_pie)
    largo_corona = 3.0 * altura if largo_corona is None else float(largo_corona)
    x_pie = largo_pie
    # Una cara vertical se representa con un desplazamiento mínimo para mantener x creciente.
    x_corona = x_pie + max(altura / math.tan(math.radians(angulo)), 1e-6)
    x_fin = x_corona + largo_corona
    superficie = [(0.0, 0.0), (x_pie, 0.0), (x_corona, altura), (x_fin, altura)]
    freatico = None
    if nivel_freatico is not None:
        freatico = nivel_freatico_simple(x_pie, x_corona, altura, x_fin, **nivel_freatico)
    return Talud(superficie, freatico, y_estrato_firme, gamma_agua)
