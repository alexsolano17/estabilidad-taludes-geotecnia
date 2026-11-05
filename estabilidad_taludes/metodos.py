"""Superficies de falla circulares, dovelas y factor de seguridad por equilibrio límite.

Métodos implementados:
- Fellenius u ordinario (Fellenius, 1936): equilibrio de momentos sin fuerzas entre dovelas.
- Bishop simplificado (Bishop, 1955): equilibrio de momentos y equilibrio vertical de cada
  dovela, despreciando el corte entre dovelas. Se resuelve por iteración de punto fijo.

Las funciones están vectorizadas: los arreglos de las dovelas tienen forma (..., n) y los
parámetros c y tan(φ) se combinan con la forma inicial (...) mediante broadcasting de numpy.
Así se evalúan miles de círculos (búsqueda) o miles de muestras (Monte Carlo) a la vez.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .modelo import Material, Talud

MA_MIN = 0.2  # valor mínimo admitido de m_alfa en Bishop (criterio práctico habitual)


def circulo_por_cuerda(talud: Talud, x_salida, x_entrada, psi):
    """Centro y radio del círculo que corta el terreno en x_salida (abajo) y x_entrada (arriba).

    psi es el semiángulo central del arco en grados (0 < psi < 90): un psi pequeño da un arco
    casi plano y uno grande, un arco profundo. El centro queda por encima de la cuerda.
    Devuelve (xc, yc, R) con la forma de los argumentos tras el broadcasting.
    """
    x1, x2, psi_rad = np.broadcast_arrays(
        np.asarray(x_salida, float), np.asarray(x_entrada, float), np.radians(np.asarray(psi, float))
    )
    y1, y2 = talud.y_terreno(x1), talud.y_terreno(x2)
    dx, dy = x2 - x1, y2 - y1
    largo = np.hypot(dx, dy)
    with np.errstate(divide="ignore", invalid="ignore"):
        tx, ty = dx / largo, dy / largo
        radio = 0.5 * largo / np.sin(psi_rad)
        d = 0.5 * largo / np.tan(psi_rad)
    xc = 0.5 * (x1 + x2) - d * ty
    yc = 0.5 * (y1 + y2) + d * tx
    return xc, yc, radio


@dataclass
class Dovelas:
    """Dovelas de una o varias superficies circulares (arreglos de forma (..., n))."""

    x: np.ndarray  # abscisa del centro de cada dovela [m]
    ancho: np.ndarray  # ancho b de las dovelas, forma (..., 1) [m]
    y_base: np.ndarray  # cota de la base (sobre el círculo) [m]
    y_sup: np.ndarray  # cota del terreno [m]
    peso: np.ndarray  # peso W por metro de talud [kN/m]
    sen_a: np.ndarray  # seno de la inclinación α de la base
    cos_a: np.ndarray  # coseno de α
    u: np.ndarray  # presión de poros en la base [kPa]
    valida: np.ndarray  # (...): la superficie es geométricamente admisible


def dovelas(
    talud: Talud,
    material: Material,
    xc,
    yc,
    radio,
    x_salida,
    x_entrada,
    n: int = 100,
    prof_min: float = 0.0,
    tol: float = 1e-6,
) -> Dovelas:
    """Divide la masa deslizante en n dovelas verticales del mismo ancho.

    La superficie es válida si el arco entra y sale del terreno sólo en sus extremos, queda
    por encima del estrato firme y su espesor máximo alcanza `prof_min`.
    El peso usa la altura en el centro de cada dovela y la presión de poros es hidrostática:
    u = γw · (cota del nivel freático − cota de la base), recortando el nivel al terreno.
    """
    xc, yc, radio, x1, x2 = np.broadcast_arrays(
        *(np.asarray(a, float) for a in (xc, yc, radio, x_salida, x_entrada))
    )
    fraccion = (np.arange(n) + 0.5) / n
    ancho = ((x2 - x1) / n)[..., None]
    x = x1[..., None] + fraccion * (x2 - x1)[..., None]
    dx = x - xc[..., None]
    r = radio[..., None]
    with np.errstate(invalid="ignore"):
        disc = r**2 - dx**2
        raiz = np.sqrt(np.clip(disc, 0.0, None))
    y_base = yc[..., None] - raiz
    y_sup = talud.y_terreno(x)
    espesor = y_sup - y_base

    y1, y2 = talud.y_terreno(x1), talud.y_terreno(x2)
    valida = np.isfinite(radio) & (radio > 0) & (x2 > x1) & (y1 < yc) & (y2 < yc)
    valida &= np.all(disc > 0, axis=-1) & np.all(espesor >= -tol, axis=-1)
    valida &= np.max(espesor, axis=-1) >= max(prof_min, tol)
    if talud.y_estrato_firme is not None:
        valida &= np.min(y_base, axis=-1) >= talud.y_estrato_firme - tol

    espesor = np.clip(espesor, 0.0, None)
    if talud.nivel_freatico is not None:
        y_agua = np.minimum(talud.y_freatico(x), y_sup)
        saturado = np.clip(y_agua - y_base, 0.0, espesor)
    else:
        saturado = np.zeros_like(espesor)
    peso = ancho * (material.gamma * (espesor - saturado) + material.gamma_saturado * saturado)
    with np.errstate(invalid="ignore", divide="ignore"):
        sen_a = dx / r
        cos_a = raiz / r
    u = talud.gamma_agua * saturado
    return Dovelas(x, ancho, y_base, y_sup, peso, sen_a, cos_a, u, valida)


def _por_dovela(valor) -> np.ndarray:
    """Agrega un eje final para combinar un parámetro con los arreglos de dovelas."""
    return np.asarray(valor, float)[..., None]


def fs_fellenius(dov: Dovelas, c, tan_phi) -> np.ndarray:
    """FS = Σ[c'·l + (W·cosα − u·l)·tanφ'] / Σ W·senα  (fuerza normal efectiva ≥ 0)."""
    c_, t_ = _por_dovela(c), _por_dovela(tan_phi)
    with np.errstate(invalid="ignore", divide="ignore"):
        largo = dov.ancho / dov.cos_a
        normal_ef = np.maximum(dov.peso * dov.cos_a - dov.u * largo, 0.0)
        resistente = np.sum(c_ * largo + normal_ef * t_, axis=-1)
        motor = np.sum(dov.peso * dov.sen_a, axis=-1)
        fs = resistente / motor
    return np.where(motor > 0, fs, np.nan)


def fs_bishop(dov: Dovelas, c, tan_phi, tol: float = 1e-9, max_iter: int = 200):
    """Bishop simplificado por iteración de punto fijo.

    FS = Σ{[c'·b + (W − u·b)·tanφ'] / mα} / Σ W·senα,  con mα = cosα + senα·tanφ'/FS.
    Devuelve (fs, ma_min, convergio): ma_min es el menor mα de la superficie, útil para
    descartar soluciones con dovelas de salida muy empinadas.
    """
    c_, t_ = _por_dovela(c), _por_dovela(tan_phi)
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        motor = np.sum(dov.peso * dov.sen_a, axis=-1)
        numerador = c_ * dov.ancho + np.maximum(dov.peso - dov.u * dov.ancho, 0.0) * t_
        fs = fs_fellenius(dov, c, tan_phi)
        fs = np.where(np.isfinite(fs) & (fs > 0), fs, 1.0)
        convergio = np.zeros(fs.shape, dtype=bool)
        for _ in range(max_iter):
            ma = dov.cos_a + dov.sen_a * t_ / fs[..., None]
            nuevo = np.sum(numerador / ma, axis=-1) / motor
            convergio = np.abs(nuevo - fs) <= tol * np.maximum(1.0, np.abs(nuevo))
            fs = nuevo
            if np.all(convergio | ~np.isfinite(fs)):
                break
        ma_min = np.min(dov.cos_a + dov.sen_a * t_ / fs[..., None], axis=-1)
    fs = np.where(motor > 0, fs, np.nan)
    return fs, ma_min, convergio


def factor_seguridad(dov: Dovelas, c, tan_phi, metodo: str = "bishop", ma_min: float = MA_MIN):
    """FS de cada superficie; +inf donde la superficie o la solución no son admisibles."""
    if metodo == "fellenius":
        fs = fs_fellenius(dov, c, tan_phi)
        ok = dov.valida & np.isfinite(fs) & (fs > 0)
    elif metodo == "bishop":
        fs, ma, convergio = fs_bishop(dov, c, tan_phi)
        ok = dov.valida & convergio & np.isfinite(fs) & (fs > 0) & (ma >= ma_min)
    else:
        raise ValueError(f"Método desconocido: {metodo!r} (use 'bishop' o 'fellenius').")
    return np.where(ok, fs, np.inf)


def evaluar_circulos(
    talud: Talud,
    material: Material,
    x_salida,
    x_entrada,
    psi,
    metodo: str = "bishop",
    n: int = 100,
    ma_min: float = MA_MIN,
    prof_min: float = 0.0,
) -> np.ndarray:
    """FS de los círculos definidos por (x_salida, x_entrada, psi); +inf si no son válidos."""
    xc, yc, radio = circulo_por_cuerda(talud, x_salida, x_entrada, psi)
    dov = dovelas(talud, material, xc, yc, radio, x_salida, x_entrada, n=n, prof_min=prof_min)
    return factor_seguridad(dov, material.c, material.tan_phi, metodo, ma_min)
