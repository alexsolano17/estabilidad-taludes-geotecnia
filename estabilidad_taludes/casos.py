"""Lectura de casos de análisis desde archivos JSON (ver la carpeta casos/)."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .busqueda import ConfigBusqueda
from .modelo import GAMMA_AGUA, Material, Talud, talud_simple
from .probabilistico import Variable, matriz_correlacion


@dataclass(frozen=True)
class Criterios:
    """Criterios de aceptación: FS mínimo y probabilidad de falla máxima."""

    fs_min: float
    pf_max: float
    descripcion: str = ""


@dataclass
class Caso:
    nombre: str
    descripcion: str
    geometria: dict
    material: Material
    nivel_freatico: dict | None = None
    y_estrato_firme: float | None = None
    gamma_agua: float = GAMMA_AGUA
    busqueda: ConfigBusqueda = field(default_factory=ConfigBusqueda)
    variables: list[Variable] = field(default_factory=list)
    rho_c_phi: float = 0.0
    simulaciones: int = 20000
    semilla: int = 2026
    criterios: Criterios | None = None

    @property
    def es_simple(self) -> bool:
        """True si la geometría es un talud simple (altura y ángulo) y no una polilínea."""
        return "puntos" not in self.geometria

    @property
    def correlacion(self) -> np.ndarray:
        return matriz_correlacion(self.variables, self.rho_c_phi)

    def talud(self, angulo: float | None = None, con_agua: bool = True) -> Talud:
        """Construye el talud; `angulo` reemplaza el ángulo de un talud simple."""
        agua = self.nivel_freatico if con_agua else None
        puntos_agua = agua.get("puntos") if agua else None
        if self.es_simple:
            g = self.geometria
            relativo = None if (agua is None or puntos_agua is not None) else agua
            t = talud_simple(
                g["altura"],
                g["angulo"] if angulo is None else angulo,
                g.get("largo_pie"),
                g.get("largo_corona"),
                relativo,
                self.y_estrato_firme,
                self.gamma_agua,
            )
            if puntos_agua is not None:
                t = Talud(t.superficie, puntos_agua, self.y_estrato_firme, self.gamma_agua)
            return t
        if angulo is not None:
            raise ValueError("Cambiar el ángulo sólo es posible con una geometría de talud simple.")
        return Talud(self.geometria["puntos"], puntos_agua, self.y_estrato_firme, self.gamma_agua)


def _variables(datos: dict, material: Material) -> list[Variable]:
    variables = []
    for nombre in ("c", "phi"):
        if nombre not in datos:
            continue
        d = datos[nombre]
        media = float(d.get("media", getattr(material, nombre)))
        desviacion = float(d["desviacion"]) if "desviacion" in d else float(d["cov"]) * media
        variables.append(
            Variable(nombre, d.get("distribucion", "normal"), media, desviacion, d.get("minimo", 0.0))
        )
    return variables


def caso_desde_dict(datos: dict) -> Caso:
    material = Material(**datos["material"])
    busqueda = ConfigBusqueda(
        **{k: (tuple(v) if isinstance(v, list) else v) for k, v in datos.get("busqueda", {}).items()}
    )
    prob = datos.get("probabilistico", {})
    crit = datos.get("criterios")
    return Caso(
        nombre=datos["nombre"],
        descripcion=datos.get("descripcion", ""),
        geometria=datos["geometria"],
        material=material,
        nivel_freatico=datos.get("nivel_freatico"),
        y_estrato_firme=datos.get("y_estrato_firme"),
        gamma_agua=datos.get("gamma_agua", GAMMA_AGUA),
        busqueda=busqueda,
        variables=_variables(prob.get("variables", {}), material),
        rho_c_phi=float(prob.get("correlacion_c_phi", 0.0)),
        simulaciones=int(prob.get("simulaciones", 20000)),
        semilla=int(prob.get("semilla", 2026)),
        criterios=Criterios(**crit) if crit else None,
    )


def cargar_caso(ruta: str | Path) -> Caso:
    return caso_desde_dict(json.loads(Path(ruta).read_text(encoding="utf-8")))
