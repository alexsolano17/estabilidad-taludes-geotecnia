"""Análisis probabilístico del factor de seguridad.

- Muestreo de c' y φ' con distribución normal o lognormal y correlación impuesta en el
  espacio normal estándar (cópula gaussiana, descomposición de Cholesky).
- Monte Carlo con Bishop simplificado sobre un conjunto de superficies: la crítica
  determinista y otras candidatas. Para cada muestra se reporta el FS en la superficie fija y
  el mínimo entre todas las candidatas, porque la superficie crítica cambia con c' y φ'.
- Estimaciones puntuales de Rosenblueth (1975) como control independiente de media y desviación.
- Probabilidad de falla P(FS < 1) con intervalo de Wilson e índices de confiabilidad β.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from statistics import NormalDist

import numpy as np

from .busqueda import ConfigBusqueda, ResultadoBusqueda, Superficie, buscar_circulo_critico, superficies_candidatas
from .metodos import MA_MIN, dovelas, fs_bishop
from .modelo import Material, Talud

_NORMAL = NormalDist()
VARIABLES_ADMITIDAS = ("c", "phi")


@dataclass(frozen=True)
class Variable:
    """Parámetro aleatorio del material ('c' en kPa o 'phi' en grados)."""

    nombre: str
    distribucion: str  # 'normal' o 'lognormal'
    media: float
    desviacion: float
    minimo: float | None = 0.0  # recorte inferior aplicado sólo a la distribución normal

    def __post_init__(self) -> None:
        if self.nombre not in VARIABLES_ADMITIDAS:
            raise ValueError(f"Variable no admitida: {self.nombre!r} (use 'c' o 'phi').")
        if self.distribucion not in ("normal", "lognormal"):
            raise ValueError(f"Distribución no admitida: {self.distribucion!r}.")
        if self.desviacion < 0:
            raise ValueError("La desviación estándar no puede ser negativa.")
        if self.distribucion == "lognormal" and self.media <= 0:
            raise ValueError("Una variable lognormal necesita media positiva.")

    @property
    def cov(self) -> float:
        return self.desviacion / self.media

    def desde_normal_estandar(self, z: np.ndarray) -> np.ndarray:
        """Transforma valores normales estándar a la distribución de la variable."""
        if self.distribucion == "normal":
            x = self.media + self.desviacion * z
            return x if self.minimo is None else np.maximum(x, self.minimo)
        s = math.sqrt(math.log1p(self.cov**2))
        m = math.log(self.media) - 0.5 * s**2
        return np.exp(m + s * z)


def matriz_correlacion(variables: list[Variable], rho_c_phi: float = 0.0) -> np.ndarray:
    """Matriz de correlación con ρ(c', φ') si ambas variables son aleatorias."""
    k = len(variables)
    corr = np.eye(k)
    nombres = [v.nombre for v in variables]
    if "c" in nombres and "phi" in nombres:
        i, j = nombres.index("c"), nombres.index("phi")
        corr[i, j] = corr[j, i] = rho_c_phi
    return corr


def muestrear(variables: list[Variable], correlacion: np.ndarray, n: int, semilla: int | None):
    """Muestras correlacionadas. Devuelve (dict nombre → arreglo, número de valores recortados)."""
    chol = np.linalg.cholesky(np.asarray(correlacion, float))
    rng = np.random.default_rng(semilla)
    z = rng.standard_normal((n, len(variables))) @ chol.T
    muestras = {v.nombre: v.desde_normal_estandar(z[:, j]) for j, v in enumerate(variables)}
    recortadas = sum(
        int(np.sum(v.media + v.desviacion * z[:, j] < v.minimo))
        for j, v in enumerate(variables)
        if v.distribucion == "normal" and v.minimo is not None
    )
    return muestras, recortadas


def fs_en_superficies(
    talud: Talud,
    material: Material,
    superficies: list[Superficie],
    c: np.ndarray,
    phi: np.ndarray,
    n_dovelas: int = 100,
    ma_min: float = MA_MIN,
    lote: int = 1000,
):
    """FS de Bishop para cada muestra (filas) y superficie (columnas), y su máscara de validez."""
    x1 = np.array([s.x_salida for s in superficies])
    x2 = np.array([s.x_entrada for s in superficies])
    xc = np.array([s.xc for s in superficies])
    yc = np.array([s.yc for s in superficies])
    radio = np.array([s.radio for s in superficies])
    dov = dovelas(talud, material, xc, yc, radio, x1, x2, n=n_dovelas)
    c = np.asarray(c, float)
    tan_phi = np.tan(np.radians(np.asarray(phi, float)))
    fs = np.empty((len(c), len(superficies)))
    ok = np.empty(fs.shape, dtype=bool)
    for i in range(0, len(c), lote):
        f, ma, convergio = fs_bishop(dov, c[i : i + lote, None], tan_phi[i : i + lote, None])
        fs[i : i + lote] = f
        ok[i : i + lote] = convergio & (ma >= ma_min) & dov.valida & np.isfinite(f) & (f > 0)
    return fs, ok


@dataclass
class ResultadoMonteCarlo:
    c: np.ndarray
    phi: np.ndarray
    fs_fija: np.ndarray  # FS en la superficie crítica determinista (columna 0)
    fs_min: np.ndarray  # mínimo entre todas las superficies candidatas
    superficie_min: np.ndarray  # índice de la superficie que da el mínimo
    superficies: list[Superficie]
    avisos_bishop: int  # muestras sin solución admisible de Bishop en la superficie fija
    recortadas: int  # valores normales recortados en el mínimo admitido


def monte_carlo(
    talud: Talud,
    material: Material,
    superficies: list[Superficie],
    variables: list[Variable],
    correlacion: np.ndarray,
    n: int,
    semilla: int | None = None,
    n_dovelas: int = 100,
    ma_min: float = MA_MIN,
) -> ResultadoMonteCarlo:
    """Simulación de Monte Carlo. La primera superficie debe ser la crítica determinista."""
    muestras, recortadas = muestrear(variables, correlacion, n, semilla)
    c = muestras.get("c", np.full(n, material.c))
    phi = muestras.get("phi", np.full(n, material.phi))
    fs, ok = fs_en_superficies(talud, material, superficies, c, phi, n_dovelas, ma_min)
    fs_fija = fs[:, 0].copy()
    validos = np.where(ok, fs, np.inf)
    idx = np.argmin(validos, axis=1)
    fs_min = validos[np.arange(n), idx]
    sin_solucion = ~np.isfinite(fs_min)
    fs_min[sin_solucion] = fs_fija[sin_solucion]
    idx[sin_solucion] = 0
    return ResultadoMonteCarlo(
        c=c,
        phi=phi,
        fs_fija=fs_fija,
        fs_min=fs_min,
        superficie_min=idx,
        superficies=list(superficies),
        avisos_bishop=int(np.sum(~ok[:, 0])),
        recortadas=recortadas,
    )


def _distinta(s: Superficie, lista: list[Superficie], escala: float) -> bool:
    return all(
        max(abs(s.x_salida - o.x_salida) / escala, abs(s.x_entrada - o.x_entrada) / escala, abs(s.psi - o.psi) / 90.0)
        > 0.01
        for o in lista
    )


def superficies_para_monte_carlo(
    talud: Talud,
    material: Material,
    variables: list[Variable],
    busqueda: ResultadoBusqueda,
    config: ConfigBusqueda | None = None,
    k_malla: int = 10,
) -> list[Superficie]:
    """Crítica determinista + críticas en los puntos μ ± σ de Rosenblueth + candidatas de la malla.

    Las búsquedas en μ ± σ cubren las superficies que se vuelven críticas cuando c' y φ'
    se alejan de su media (por ejemplo, círculos más someros con c' baja).
    """
    _, _, altura = talud.cara()
    lista = [busqueda.critica]
    for signos in itertools.product((1, -1), repeat=len(variables)):
        cambios = {}
        for v, s in zip(variables, signos):
            valor = v.media + s * v.desviacion
            cambios[v.nombre] = max(valor, 0.0)
        r = buscar_circulo_critico(talud, material.con(**cambios), "bishop", config)
        for sup in [r.critica, *superficies_candidatas(talud, r, k=3)]:
            if _distinta(sup, lista, altura):
                lista.append(sup)
    for sup in superficies_candidatas(talud, busqueda, k=k_malla):
        if _distinta(sup, lista, altura):
            lista.append(sup)
    return lista


@dataclass(frozen=True)
class Estadisticas:
    n: int
    media: float
    desviacion: float
    minimo: float
    p05: float
    p50: float
    p95: float
    fallas: int
    pf: float  # P(FS < 1) por conteo
    pf_ic95: tuple[float, float]  # intervalo de Wilson al 95 %
    error_pf: float  # error estándar binomial
    beta_normal: float
    beta_lognormal: float
    pf_normal: float  # Φ(−β) suponiendo FS normal
    pf_lognormal: float  # Φ(−β) suponiendo FS lognormal


def intervalo_wilson(fallas: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """Intervalo de confianza de Wilson para una proporción (válido también con 0 fallas)."""
    p = fallas / n
    den = 1.0 + z * z / n
    centro = (p + z * z / (2 * n)) / den
    medio = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    bajo = 0.0 if fallas == 0 else max(0.0, centro - medio)
    alto = 1.0 if fallas == n else min(1.0, centro + medio)
    return bajo, alto


def estadisticas(fs: np.ndarray) -> Estadisticas:
    fs = np.asarray(fs, float)
    n = fs.size
    media = float(np.mean(fs))
    desv = float(np.std(fs, ddof=1))
    fallas = int(np.sum(fs < 1.0))
    pf = fallas / n
    beta_n = (media - 1.0) / desv
    v2 = (desv / media) ** 2
    beta_ln = math.log(media / math.sqrt(1.0 + v2)) / math.sqrt(math.log1p(v2))
    p05, p50, p95 = (float(x) for x in np.percentile(fs, [5, 50, 95]))
    return Estadisticas(
        n=n,
        media=media,
        desviacion=desv,
        minimo=float(np.min(fs)),
        p05=p05,
        p50=p50,
        p95=p95,
        fallas=fallas,
        pf=pf,
        pf_ic95=intervalo_wilson(fallas, n),
        error_pf=math.sqrt(pf * (1 - pf) / n),
        beta_normal=beta_n,
        beta_lognormal=beta_ln,
        pf_normal=_NORMAL.cdf(-beta_n),
        pf_lognormal=_NORMAL.cdf(-beta_ln),
    )


def pf_acumulada(fs: np.ndarray) -> np.ndarray:
    """Estimación de P(FS < 1) a medida que se agregan muestras (para ver la convergencia)."""
    fs = np.asarray(fs, float)
    return np.cumsum(fs < 1.0) / np.arange(1, fs.size + 1)


@dataclass(frozen=True)
class ResultadoRosenblueth:
    valores: list[dict]  # valores de las variables en cada punto
    pesos: np.ndarray
    fs: np.ndarray
    media: float
    desviacion: float


def rosenblueth(funcion, variables: list[Variable], correlacion: np.ndarray) -> ResultadoRosenblueth:
    """Estimaciones puntuales de Rosenblueth (1975) con 2^k puntos μ ± σ.

    Pesos para variables correlacionadas: P = [1 + Σ_{i<j} s_i·s_j·ρ_ij] / 2^k.
    `funcion(c, phi)` recibe arreglos y devuelve el FS de cada punto.
    """
    corr = np.asarray(correlacion, float)
    k = len(variables)
    valores, pesos = [], []
    for signos in itertools.product((1, -1), repeat=k):
        suma = sum(signos[i] * signos[j] * corr[i, j] for i in range(k) for j in range(i + 1, k))
        pesos.append((1.0 + suma) / 2**k)
        valores.append({v.nombre: v.media + s * v.desviacion for v, s in zip(variables, signos)})
    c = np.array([v.get("c", np.nan) for v in valores])
    phi = np.array([v.get("phi", np.nan) for v in valores])
    fs = np.asarray(funcion(c, phi), float)
    pesos = np.array(pesos)
    media = float(np.sum(pesos * fs))
    desviacion = float(math.sqrt(max(np.sum(pesos * fs**2) - media**2, 0.0)))
    return ResultadoRosenblueth(valores, pesos, fs, media, desviacion)
