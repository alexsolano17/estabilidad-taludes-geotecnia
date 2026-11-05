"""Búsqueda del círculo de falla crítico (mínimo factor de seguridad).

Cada círculo se define por tres parámetros: abscisa de salida (cerca del pie), abscisa de
entrada (cerca de la corona) y semiángulo central del arco psi. La búsqueda recorre primero
una malla de esas tres variables y luego refina los mejores puntos con una búsqueda por
patrones (pasos en cada eje que se reducen a la mitad cuando ya no hay mejora).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .metodos import MA_MIN, circulo_por_cuerda, evaluar_circulos
from .modelo import Material, Talud


@dataclass(frozen=True)
class Superficie:
    """Círculo de falla y su factor de seguridad."""

    x_salida: float
    x_entrada: float
    psi: float
    xc: float
    yc: float
    radio: float
    fs: float


@dataclass
class ConfigBusqueda:
    """Rangos y resolución de la búsqueda. Los rangos None se calculan a partir del perfil."""

    x_salida: tuple[float, float] | None = None
    x_entrada: tuple[float, float] | None = None
    psi: tuple[float, float] = (3.0, 80.0)
    n_salida: int = 25
    n_entrada: int = 25
    n_psi: int = 20
    n_dovelas_malla: int = 50
    n_dovelas: int = 100
    n_refinar: int = 4
    ma_min: float = MA_MIN
    prof_min: float = 0.0


@dataclass
class ResultadoBusqueda:
    critica: Superficie
    metodo: str
    n_dovelas: int
    evaluadas: int
    validas: int
    malla: np.ndarray = field(repr=False)  # filas (x_salida, x_entrada, psi, fs) válidas


def rangos_por_defecto(talud: Talud) -> tuple[tuple[float, float], tuple[float, float]]:
    """Salidas desde 1,5·H antes del pie hasta 3/4 de la cara; entradas desde 1/4 de la cara
    hasta 2·H detrás de la corona (recortados a los límites del perfil)."""
    x_pie, x_corona, altura = talud.cara()
    cara = x_corona - x_pie
    salida = (max(talud.x_min, x_pie - 1.5 * altura), x_pie + 0.75 * cara)
    entrada = (x_pie + 0.25 * cara, min(talud.x_max, x_corona + 2.0 * altura))
    return salida, entrada


def _superficie(talud: Talud, p: np.ndarray, fs: float) -> Superficie:
    xc, yc, radio = circulo_por_cuerda(talud, p[0], p[1], p[2])
    return Superficie(float(p[0]), float(p[1]), float(p[2]), float(xc), float(yc), float(radio), float(fs))


def _evaluar_lotes(talud, material, params, metodo, n, cfg, lote=4000) -> np.ndarray:
    fs = np.empty(len(params))
    for i in range(0, len(params), lote):
        p = params[i : i + lote]
        fs[i : i + lote] = evaluar_circulos(
            talud, material, p[:, 0], p[:, 1], p[:, 2], metodo, n, cfg.ma_min, cfg.prof_min
        )
    return fs


def _refinar(funcion, p0: np.ndarray, paso0: np.ndarray, paso_min: np.ndarray):
    """Búsqueda por patrones: prueba ±paso en cada eje; si nada mejora, divide el paso por 2."""
    p = np.array(p0, float)
    fp = funcion(p[None, :])[0]
    paso = np.array(paso0, float)
    for _ in range(2000):
        if np.all(paso < paso_min):
            break
        vecinos = np.concatenate([p + np.diag(paso), p - np.diag(paso)])
        fv = funcion(vecinos)
        j = int(np.argmin(fv))
        if fv[j] < fp - 1e-12:
            p, fp = vecinos[j], fv[j]
        else:
            paso = paso * 0.5
    return p, fp


def _distintos(params: np.ndarray, escala: np.ndarray, k: int, separacion: float) -> list[int]:
    """Índices (en el orden dado) de hasta k puntos separados al menos `separacion` pasos."""
    elegidos: list[int] = []
    for i in range(len(params)):
        if all(np.max(np.abs(params[i] - params[j]) / escala) >= separacion for j in elegidos):
            elegidos.append(i)
            if len(elegidos) == k:
                break
    return elegidos


def buscar_circulo_critico(
    talud: Talud,
    material: Material,
    metodo: str = "bishop",
    config: ConfigBusqueda | None = None,
) -> ResultadoBusqueda:
    """Encuentra el círculo con el menor factor de seguridad."""
    cfg = config or ConfigBusqueda()
    salida_def, entrada_def = rangos_por_defecto(talud)
    salida = cfg.x_salida or salida_def
    entrada = cfg.x_entrada or entrada_def

    xs = np.linspace(*salida, cfg.n_salida)
    xe = np.linspace(*entrada, cfg.n_entrada)
    ps = np.linspace(*cfg.psi, cfg.n_psi)
    malla = np.stack(np.meshgrid(xs, xe, ps, indexing="ij"), axis=-1).reshape(-1, 3)
    malla = malla[malla[:, 1] > malla[:, 0]]
    fs_malla = _evaluar_lotes(talud, material, malla, metodo, cfg.n_dovelas_malla, cfg)
    validos = np.isfinite(fs_malla)
    if not np.any(validos):
        raise RuntimeError("Ningún círculo de la malla es válido: revise rangos y geometría.")

    orden = np.argsort(fs_malla)
    orden = orden[np.isfinite(fs_malla[orden])]
    paso0 = np.array(
        [
            np.ptp(xs) / max(cfg.n_salida - 1, 1),
            np.ptp(xe) / max(cfg.n_entrada - 1, 1),
            np.ptp(ps) / max(cfg.n_psi - 1, 1),
        ]
    )
    paso0 = np.where(paso0 > 0, paso0, 1.0)
    _, _, altura = talud.cara()
    paso_min = np.array([1e-3 * altura, 1e-3 * altura, 0.02])

    def funcion(p: np.ndarray) -> np.ndarray:
        p = np.atleast_2d(p)
        return _evaluar_lotes(talud, material, p, metodo, cfg.n_dovelas, cfg)

    mejor_p, mejor_fs = None, np.inf
    for i in _distintos(malla[orden], paso0, cfg.n_refinar, separacion=2.0):
        p, fs = _refinar(funcion, malla[orden][i], paso0, paso_min)
        if fs < mejor_fs:
            mejor_p, mejor_fs = p, fs

    tabla = np.column_stack([malla[validos], fs_malla[validos]])
    return ResultadoBusqueda(
        critica=_superficie(talud, mejor_p, mejor_fs),
        metodo=metodo,
        n_dovelas=cfg.n_dovelas,
        evaluadas=len(malla),
        validas=int(validos.sum()),
        malla=tabla,
    )


def superficies_candidatas(
    talud: Talud, resultado: ResultadoBusqueda, k: int = 10, separacion: float = 0.05
) -> list[Superficie]:
    """Superficies de la malla con FS bajo y distintas entre sí (separación relativa a H)."""
    _, _, altura = talud.cara()
    tabla = resultado.malla[np.argsort(resultado.malla[:, 3])]
    escala = np.array([altura, altura, 90.0])
    idx = _distintos(tabla[:, :3], escala, k, separacion)
    return [_superficie(talud, tabla[i, :3], tabla[i, 3]) for i in idx]
