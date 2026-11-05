"""Flujo completo de un caso: búsqueda determinista, Monte Carlo, Rosenblueth y diseño."""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np

from .busqueda import ResultadoBusqueda, buscar_circulo_critico
from .casos import Caso
from .modelo import Talud
from .probabilistico import (
    Estadisticas,
    ResultadoMonteCarlo,
    ResultadoRosenblueth,
    estadisticas,
    fs_en_superficies,
    monte_carlo,
    rosenblueth,
    superficies_para_monte_carlo,
)


@dataclass
class ResultadoAnalisis:
    caso: Caso
    angulo: float | None
    talud: Talud
    bishop: ResultadoBusqueda
    fellenius: ResultadoBusqueda
    fs_seco: float | None  # Bishop sin nivel freático (sólo si el caso tiene agua)
    monte_carlo: ResultadoMonteCarlo | None
    est_fija: Estadisticas | None  # FS en la superficie crítica determinista
    est_min: Estadisticas | None  # FS mínimo entre superficies candidatas
    rosenblueth: ResultadoRosenblueth | None
    segundos: float

    @property
    def fs(self) -> float:
        return self.bishop.critica.fs

    @property
    def pf(self) -> float | None:
        return None if self.est_min is None else self.est_min.pf

    def cumple(self) -> dict | None:
        """Comparación con los criterios de aceptación del caso (None si no hay criterios)."""
        crit = self.caso.criterios
        if crit is None:
            return None
        cumple_fs = self.fs >= crit.fs_min
        cumple_pf = None if self.pf is None else self.pf <= crit.pf_max
        return {
            "fs": cumple_fs,
            "pf": cumple_pf,
            "global": bool(cumple_fs and (cumple_pf is None or cumple_pf)),
        }


def analizar(
    caso: Caso,
    angulo: float | None = None,
    simulaciones: int | None = None,
    probabilistico: bool = True,
) -> ResultadoAnalisis:
    """Ejecuta el análisis determinista y, si hay variables aleatorias, el probabilístico."""
    inicio = time.perf_counter()
    talud = caso.talud(angulo)
    bishop = buscar_circulo_critico(talud, caso.material, "bishop", caso.busqueda)
    fellenius = buscar_circulo_critico(talud, caso.material, "fellenius", caso.busqueda)
    fs_seco = None
    if talud.nivel_freatico is not None:
        seco = caso.talud(angulo, con_agua=False)
        fs_seco = buscar_circulo_critico(seco, caso.material, "bishop", caso.busqueda).critica.fs

    mc = est_fija = est_min = ros = None
    if probabilistico and caso.variables:
        superficies = superficies_para_monte_carlo(
            talud, caso.material, caso.variables, bishop, caso.busqueda
        )
        mc = monte_carlo(
            talud,
            caso.material,
            superficies,
            caso.variables,
            caso.correlacion,
            simulaciones or caso.simulaciones,
            caso.semilla,
            caso.busqueda.n_dovelas,
            caso.busqueda.ma_min,
        )
        est_fija = estadisticas(mc.fs_fija)
        est_min = estadisticas(mc.fs_min)

        def fs_superficie_critica(c, phi):
            c = np.where(np.isnan(c), caso.material.c, c)
            phi = np.where(np.isnan(phi), caso.material.phi, phi)
            fs, _ = fs_en_superficies(
                talud, caso.material, [bishop.critica], c, phi, caso.busqueda.n_dovelas
            )
            return fs[:, 0]

        ros = rosenblueth(fs_superficie_critica, caso.variables, caso.correlacion)

    return ResultadoAnalisis(
        caso=caso,
        angulo=angulo if angulo is not None else caso.geometria.get("angulo"),
        talud=talud,
        bishop=bishop,
        fellenius=fellenius,
        fs_seco=fs_seco,
        monte_carlo=mc,
        est_fija=est_fija,
        est_min=est_min,
        rosenblueth=ros,
        segundos=time.perf_counter() - inicio,
    )


@dataclass(frozen=True)
class FilaDiseno:
    angulo: float
    fs: float
    fs_seco: float | None
    pf: float | None
    pf_ic95: tuple[float, float] | None
    cumple_fs: bool | None
    cumple_pf: bool | None


def barrido_angulos(caso: Caso, angulos, simulaciones: int | None = None) -> list[FilaDiseno]:
    """Analiza el talud para cada ángulo de cara (requiere geometría de talud simple)."""
    filas = []
    for angulo in angulos:
        r = analizar(caso, float(angulo), simulaciones)
        cumple = r.cumple() or {}
        filas.append(
            FilaDiseno(
                angulo=float(angulo),
                fs=r.fs,
                fs_seco=r.fs_seco,
                pf=r.pf,
                pf_ic95=None if r.est_min is None else r.est_min.pf_ic95,
                cumple_fs=cumple.get("fs"),
                cumple_pf=cumple.get("pf"),
            )
        )
    return filas


def angulo_recomendado(filas: list[FilaDiseno]) -> float | None:
    """Mayor ángulo tal que él y todos los ángulos menores cumplen FS y PF."""
    recomendado = None
    for fila in sorted(filas, key=lambda f: f.angulo):
        if fila.cumple_fs and fila.cumple_pf is not False:
            recomendado = fila.angulo
        else:
            break
    return recomendado
