"""Estabilidad de taludes por equilibrio límite con análisis probabilístico."""

from .busqueda import ConfigBusqueda, ResultadoBusqueda, Superficie, buscar_circulo_critico
from .metodos import circulo_por_cuerda, dovelas, evaluar_circulos, fs_bishop, fs_fellenius
from .modelo import Material, Talud, talud_simple

__all__ = [
    "ConfigBusqueda",
    "Material",
    "ResultadoBusqueda",
    "Superficie",
    "Talud",
    "buscar_circulo_critico",
    "circulo_por_cuerda",
    "dovelas",
    "evaluar_circulos",
    "fs_bishop",
    "fs_fellenius",
    "talud_simple",
]

__version__ = "1.0.0"
